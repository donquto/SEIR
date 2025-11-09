# ============================================================
# SEIR 拟合脚本（融合 Excel + 稳健修正版）
# - 自动解析 Excel 参数（R0, sigma, gamma, population）
# - beta_shape 归一化（mean=1）
# - 自适应人口规模，避免 μ 远超观测
# - L-BFGS-B 带边界优化，必要时回退 Nelder–Mead
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(deSolve); library(readxl)
})

# -------------------------------
# 可配置项
# -------------------------------
if (!exists("output_dir")) output_dir <- file.path(getwd(), "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

preproc_rds_path      <- file.path(output_dir, "preproc.rds")
excel_path_guess      <- file.path(getwd(), "流感参数.xlsx")
force_excel_population <- FALSE  # 若 TRUE 则强制使用 Excel 的人口；否则按观测自适应

eta0_init <- 0.05  # 观测比例的初值（更易收敛；如偏高可调小到 0.02）

# -------------------------------
# 解析工具
# -------------------------------
num_only <- function(x) as.numeric(gsub("[^0-9.+-eE]", "", x))

parse_range <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  if (grepl("~", x)) {
    parts <- strsplit(x, "~")[[1]]
    vals  <- suppressWarnings(as.numeric(gsub("[^0-9.+-eE]", "", parts)))
    if (length(vals) == 2 && all(is.finite(vals))) return(vals)
  }
  v <- suppressWarnings(as.numeric(gsub("[^0-9.+-eE]", "", x)))
  if (is.finite(v)) return(v)
  NA_real_
}

parse_fraction_center <- function(x) {
  x <- as.character(x)
  if (!grepl("/", x)) {
    v <- suppressWarnings(as.numeric(gsub("[^0-9.+-eE]", "", x)))
    return(ifelse(is.finite(v), v, NA_real_))
  }
  rhs <- sub("^.*?/", "", x)
  if (grepl("~", rhs)) {
    rr <- parse_range(rhs)
    if (length(rr) == 2 && all(is.finite(rr))) return(1 / mean(rr))
  }
  rhs_center <- sub("±.*$", "", rhs)
  center <- suppressWarnings(as.numeric(gsub("[^0-9.+-eE]", "", rhs_center)))
  if (is.finite(center) && center != 0) return(1 / center)
  NA_real_
}

parse_population <- function(x) {
  v <- suppressWarnings(as.numeric(gsub("[^0-9]", "", as.character(x))))
  if (is.finite(v)) return(v)
  NA_real_
}

get_excel_value <- function(df, key_pattern) {
  cand <- df %>% filter(grepl(key_pattern, 参数))
  if (nrow(cand) == 0) return(NA_character_)
  as.character(cand$数值[1])
}

# -------------------------------
# 读取预处理结果
# -------------------------------
if (!file.exists(preproc_rds_path)) {
  stop(paste0("未找到预处理文件: ", preproc_rds_path,
              "\n请先生成 preproc.rds 或修改 preproc_rds_path。"))
}
obj <- readRDS(preproc_rds_path)
surv <- as.data.frame(obj$surv)
beta_shape <- as.numeric(obj$model_total$beta_shape)

# -------------------------------
# 读取 Excel 参数（带回退）
# -------------------------------
R0_default    <- c(1.2, 1.6)
sigma_default <- 7 / 1.55  # ≈0.645
gamma_default <- 7 / 4.0   # 2~6 天中值
N_default     <- 1e6


if (file.exists(excel_path_guess)) {
  excel_df <- read_excel(excel_path_guess, sheet = 1)
  names(excel_df) <- c("参数", "数值")
  
  # R0
  R0_raw <- get_excel_value(excel_df, "再生数|R0|R₀")
  R0_parsed <- parse_range(R0_raw)
  if (all(is.finite(R0_parsed))) {
    R0_range <- if (length(R0_parsed) == 2) sort(R0_parsed) else rep(R0_parsed, 2)
  } else R0_range <- R0_default
  
  # sigma
  sigma_raw <- get_excel_value(excel_df, "潜伏率|潜伏期|σ")
  sigma_val <- parse_fraction_center(sigma_raw)
  if (!is.finite(sigma_val)) sigma_val <- suppressWarnings(num_only(sigma_raw))
  if (!is.finite(sigma_val)) sigma_val <- sigma_default
  
  # gamma
  gamma_raw <- get_excel_value(excel_df, "恢复率|恢复期|γ")
  gamma_val <- parse_fraction_center(gamma_raw)
  if (!is.finite(gamma_val)) gamma_val <- suppressWarnings(num_only(gamma_raw))
  if (!is.finite(gamma_val)) gamma_val <- gamma_default
  
  # population
  N_raw   <- get_excel_value(excel_df, "人口|总人口")
  N_excel <- parse_population(N_raw)
  if (!is.finite(N_excel) || N_excel <= 0) N_excel <- N_default
  
  message("✅ 已从 Excel 读取参数表：")
  print(excel_df)
} else {
  warning(paste0("⚠️ 未找到 Excel 文件：", excel_path_guess, "，使用默认参数。"))
  R0_range <- R0_default
  sigma_val <- sigma_default
  gamma_val <- gamma_default
  N_excel <- N_default
}

R0_mid   <- mean(R0_range)
beta0_val <- R0_mid * gamma_val

# -------------------------------
# 时间对齐 & beta_shape 归一化
# -------------------------------
# 对齐：当 beta_shape 与 weeks 不等长时，按较小长度截取并同步 surv
weeks_all <- sort(unique(surv$week_id))
T_all <- length(weeks_all)
T_b   <- length(beta_shape)
T     <- min(T_all, T_b)

if (T < T_all) {
  weeks <- weeks_all[seq_len(T)]
  surv  <- surv[surv$week_id %in% weeks, ]
} else {
  weeks <- weeks_all
}
time_idx <- seq_len(T)
id_map   <- data.frame(week_id = weeks, t_idx = time_idx)
beta_shape <- beta_shape[seq_len(T)]

# 归一化到 mean=1
bs_mean_before <- mean(beta_shape, na.rm = TRUE)
beta_shape <- beta_shape / bs_mean_before
message(sprintf("beta_shape: 原均值=%.4f, 归一化后均值=%.4f", bs_mean_before, mean(beta_shape)))

# -------------------------------
# 自适应人口规模：避免 μ 远超观测
# -------------------------------
y_peak <- max(surv$y, na.rm=TRUE)
N_suggest <- max(1e6, 20 * y_peak / eta0_init)

if (force_excel_population) {
  N_total <- N_excel
  message(sprintf("使用 Excel 人口 N=%.0f（已强制）", N_total))
} else {
  # 若 Excel 人口比建议人口大很多（>50倍），默认采用建议值以免 μ 过大
  if (is.finite(N_excel) && N_excel / N_suggest > 50) {
    N_total <- N_suggest
    message(sprintf("Excel 人口(%.0f) 远大于观测规模，改用自适应人口 N=%.0f", N_excel, N_total))
  } else {
    N_total <- if (is.finite(N_excel)) N_excel else N_suggest
    message(sprintf("使用人口 N=%.0f", N_total))
  }
}

# -------------------------------
# 傅里叶一阶季节基
# -------------------------------
w  <- 2 * pi / 52
F1 <- cbind(sin = sin(w * time_idx), cos = cos(w * time_idx))

# -------------------------------
# SEIR 动力学与模拟
# -------------------------------
seir_ode_series <- function(t, state, pars) {
  S <- state[["S"]]; E <- state[["E"]]; I <- state[["I"]]; R <- state[["R"]]
  beta_series <- pars$beta_series; sigma <- pars$sigma; gamma <- pars$gamma
  N <- S + E + I + R
  idx <- max(1L, min(length(beta_series), floor(t) + 1L))
  beta_t <- beta_series[idx]
  dS <- - beta_t * S * I / N
  dE <-   beta_t * S * I / N - sigma * E
  dI <-   sigma * E - gamma * I
  dR <-   gamma * I
  list(c(dS, dE, dI, dR), I_inc = sigma * E, N = N)
}

simulate_weekly_series <- function(init_state, beta_series, sigma, gamma){
  out <- deSolve::ode(y = init_state, times = time_idx, func = seir_ode_series,
                      parms = list(beta_series = beta_series, sigma = sigma, gamma = gamma),
                      method = "lsoda")
  df <- as.data.frame(out); names(df)[names(df)=="time"] <- "t_idx"
  if (!"I_inc" %in% names(df)) df$I_inc <- sigma * df$E
  merge(id_map, df[, c("t_idx","I_inc")], by="t_idx", all.x=TRUE)[, c("week_id","I_inc")]
}

# -------------------------------
# 参数边界投影
# -------------------------------
project_par <- function(p){
  p["log_beta0"]   <- pmax(log(0.05), pmin(log(3.0),  p["log_beta0"]))
  p["log_sigma"]   <- pmax(log(1/14), pmin(log(2),    p["log_sigma"]))  # 潜伏期 0.5~14 天
  p["log_gamma"]   <- pmax(log(1/21), pmin(log(2),    p["log_gamma"]))  # 恢复期 0.5~21 天
  p["logit_eta0"]  <- pmax(qlogis(1e-4), pmin(qlogis(0.8), p["logit_eta0"]))
  p["log_k"]       <- pmax(log(1),   pmin(log(500),  p["log_k"]))
  p["b1"]          <- pmax(-2.0, pmin(2.0, p["b1"]))
  p["b2"]          <- pmax(-2.0, pmin(2.0, p["b2"]))
  p
}

# -------------------------------
# 负二项似然
# -------------------------------
nbinom_loglik <- function(par_vec, data_df, init_state){
  p <- project_par(par_vec)
  sigma <- exp(p["log_sigma"]); gamma <- exp(p["log_gamma"])
  k     <- exp(p["log_k"])
  beta0 <- exp(p["log_beta0"])
  logit_eta0 <- p["logit_eta0"]
  b1 <- p["b1"]; b2 <- p["b2"]
  
  eta_t <- plogis(as.numeric(logit_eta0 + b1 * F1[,1] + b2 * F1[,2]))
  sim <- simulate_weekly_series(init_state, beta_series = beta0 * beta_shape,
                                sigma = sigma, gamma = gamma)
  dfm <- merge(data_df, sim, by="week_id", all.x=TRUE)
  t_idx_vec <- match(dfm$week_id, weeks)
  if (any(!is.finite(dfm$I_inc)) || max(dfm$I_inc, na.rm=TRUE) > 1e9) return(1e12)
  mu <- pmax(1e-12, eta_t[t_idx_vec] * dfm$I_inc)
  -sum(dfm$omega * dnbinom(dfm$y, size = k, mu = mu, log = TRUE), na.rm = TRUE)
}

# -------------------------------
# 初值与拟合
# -------------------------------
I0 <- max(10, round(mean(surv$y, na.rm=TRUE)))
E0 <- I0; 
R0 <- I0; 
S0 <- max(1, N_total - I0 - E0 - R0)
init_state <- c(S=S0, E=E0, I=I0, R=R0)

par_init <- c(
  log_beta0  = log(beta0_val),
  log_sigma  = log(sigma_val),
  log_gamma  = log(gamma_val),
  logit_eta0 = qlogis(eta0_init),
  log_k      = log(30),
  b1 = 0.35,
  b2 = 0.15
)

test_ll <- nbinom_loglik(par_init, data_df = as.data.frame(surv), init_state = init_state)
message("自检负对数似然：", round(test_ll, 2))

# 优先 L-BFGS-B（带边界）；失败则回退 Nelder–Mead
lower <- c(log(0.05), log(1/14), log(1/21), qlogis(1e-4), log(1),  -2,  -2)
upper <- c(log(3.0),  log(2),    log(2),    qlogis(0.8),  log(500), 2,   2)

opt <- tryCatch(
  optim(par = par_init, fn = nbinom_loglik,
        data_df = as.data.frame(surv), init_state = init_state,
        method = "L-BFGS-B", lower = lower, upper = upper,
        control = list(maxit = 3000, factr = 1e7)),
  error = function(e) e
)

if (inherits(opt, "error")) {
  warning("L-BFGS-B 失败，回退到 Nelder–Mead：", opt$message)
  opt <- optim(par = par_init, fn = nbinom_loglik,
               data_df = as.data.frame(surv), init_state = init_state,
               method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-8))
}
best <- project_par(opt$par)

# -------------------------------
# 结果与导出
# -------------------------------
pars_best <- list(
  beta0 = exp(best["log_beta0"]),
  sigma = exp(best["log_sigma"]),
  gamma = exp(best["log_gamma"]),
  k     = exp(best["log_k"]),
  eta_baseline = plogis(best["logit_eta0"]),
  b1 = best["b1"],
  b2 = best["b2"]
)
R0_est <- pars_best$beta0 / pars_best$gamma

param_tbl <- data.table::data.table(
  param = c("beta0","sigma_week","gamma_week","k","R0_approx",
            "eta_baseline","b1_sin","b2_cos","beta_source",
            "R0_input_low","R0_input_high","N_total","beta_shape_mean_before"),
  value = c(pars_best$beta0, pars_best$sigma, pars_best$gamma, pars_best$k,
            R0_est,
            pars_best$eta_baseline, pars_best$b1, pars_best$b2,
            "preproc_beta_shape(归一化) + Excel",
            R0_range[1], R0_range[2], N_total, bs_mean_before)
)
fwrite(param_tbl, file.path(output_dir, "fitted_params.csv"))
saveRDS(pars_best, file = file.path(output_dir, "pars_best.rds"))

# 预测与图
sim <- simulate_weekly_series(init_state, beta_series = pars_best$beta0 * beta_shape,
                              sigma = pars_best$sigma, gamma = pars_best$gamma)
plot_df <- merge(as.data.frame(surv), sim, by="week_id", all.x=TRUE)
eta_t <- plogis(as.numeric(best["logit_eta0"] + best["b1"]*F1[,1] + best["b2"]*F1[,2]))
t_idx_vec <- match(plot_df$week_id, weeks)
plot_df$mu <- pmax(1e-12, eta_t[t_idx_vec] * plot_df$I_inc)

p1 <- ggplot(plot_df, aes(x = week_id)) +
  geom_col(aes(y = y), alpha = 0.55, fill = "steelblue") +
  geom_line(aes(y = mu), linewidth = 1, color = "orange") +
  labs(x = "连续周序号", y = "观测阳性数 / 拟合期望",
       title = "SEIR 拟合：观测 vs 期望（beta_shape已归一化；人口自适应；Excel参数）") +
  theme_bw()

ggsave(file.path(output_dir, "fit_observed_vs_mu.png"), p1, width = 11, height = 5.2, dpi = 300)

message("✅ 拟合完成，输出目录：", normalizePath(output_dir))

# ============ 追加：更多可视化（保持你的变量命名） ============

# 1) β 形状与 β(t)（两轴一致的相对值，可直观看季节强迫）
beta_series <- as.numeric(pars_best$beta0 * beta_shape)
df_beta <- data.frame(
  week_id = weeks,
  beta_shape = as.numeric(beta_shape),
  beta_t = beta_series
)
p_beta <- ggplot(df_beta, aes(x = week_id)) +
  geom_line(aes(y = beta_shape), linewidth = 0.9, color = "#1f77b4") +
  geom_line(aes(y = scales::rescale(beta_t, to = range(beta_shape, na.rm = TRUE))),
            linewidth = 1, color = "#d62728", alpha = .9) +
  labs(title = "β 形状与 β(t)（红：β(t)经同轴缩放）",
       x = "连续周序号", y = "相对幅度（单位化）") +
  theme_bw()
ggsave(file.path(output_dir, "beta_shape_and_series.png"),
       p_beta, width = 10, height = 3.6, dpi = 300)

# 2) η_t 季节序列
eta_df <- data.frame(week_id = weeks, eta_t = eta_t)
p_eta <- ggplot(eta_df, aes(week_id, eta_t)) +
  geom_line(color = "#ff7f0e", linewidth = 1) +
  labs(title = "η_t（可见比例）的季节变化",
       x = "连续周序号", y = "η_t") +
  theme_bw()
ggsave(file.path(output_dir, "eta_series.png"),
       p_eta, width = 10, height = 3.2, dpi = 300)

# ---- 修复：保证 init_state 有 S/E/I/R 的名字 ----
init_state <- setNames(as.numeric(init_state), c("S","E","I","R"))

# ---- 安全版方程：既支持有名也支持无名 ----
seir_ode_series_safe <- function(time, state, pars){
  # 兼容：若没名字，用位置 1:4 取；若有名字，用名字取
  if (is.null(names(state)) || !all(c("S","E","I","R") %in% names(state))) {
    S <- state[1]; E <- state[2]; I <- state[3]; R <- state[4]
  } else {
    S <- state[["S"]]; E <- state[["E"]]; I <- state[["I"]]; R <- state[["R"]]
  }
  beta_series <- pars$beta_series
  sigma <- pars$sigma; gamma <- pars$gamma
  
  N <- S + E + I + R
  idx <- max(1L, min(length(beta_series), floor(time) + 1L))
  beta_t <- beta_series[idx]
  
  dS <- - beta_t * S * I / N
  dE <-   beta_t * S * I / N - sigma * E
  dI <-   sigma * E - gamma * I
  dR <-   gamma * I
  list(c(dS, dE, dI, dR), I_inc = sigma * E, N = N)
}

# ---- 直接把命名好的 init_state 传入 ----
beta_series <- as.numeric(pars_best$beta0 * beta_shape)
sim_states <- deSolve::ode(
  y = init_state,
  times = time_idx,
  func  = seir_ode_series_safe,
  parms = list(beta_series = beta_series,
               sigma = pars_best$sigma,
               gamma = pars_best$gamma),
  method = "lsoda"
) |> as.data.frame()

names(sim_states)[names(sim_states)=="time"] <- "t_idx"
sim_states <- merge(id_map, sim_states, by = "t_idx", all.x = TRUE)

# 下面照旧：做四仓室比例曲线
N_series <- with(sim_states, S + E + I + R)
plot_states <- tidyr::pivot_longer(sim_states,
                                   cols = c("E","I","R"),
                                   names_to = "comp", values_to = "value") |>
  dplyr::mutate(prop = value / rep(N_series, 3L))

p_state <- ggplot(plot_states, aes(week_id, prop, color = comp)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(E="#ff7f0e",I="#2ca02c",R="#d62728"),
                     name = "Compartment") +
  labs(title = "SEIR 四仓室相对规模（按周）",
       x = "连续周序号", y = "比例") +
  theme_bw()
ggsave(file.path(output_dir, "seir_states_3.png"),
       p_state, width = 10, height = 4.0, dpi = 300)
# 下面照旧：做四仓室比例曲线
N_series <- with(sim_states, S + E + I + R)
plot_states <- tidyr::pivot_longer(sim_states,
                                   cols = c("S","E","I","R"),
                                   names_to = "comp", values_to = "value") |>
  dplyr::mutate(prop = value / rep(N_series, 4L))

p_state <- ggplot(plot_states, aes(week_id, prop, color = comp)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(S="#1f77b4",E="#ff7f0e",I="#2ca02c",R="#d62728"),
                     name = "Compartment") +
  labs(title = "SEIR 四仓室相对规模（按周）",
       x = "连续周序号", y = "比例") +
  theme_bw()
ggsave(file.path(output_dir, "seir_states_4.png"),
       p_state, width = 10, height = 4.0, dpi = 300)

# 4) 残差与拟合优度
plot_df$resid <- with(plot_df, y - mu)
p_resid <- ggplot(plot_df, aes(week_id, resid)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_col(fill = "#6baed6", alpha = .8) +
  labs(title = "周残差（y - μ）", x = "连续周序号", y = "残差") +
  theme_bw()
ggsave(file.path(output_dir, "residuals_by_week.png"),
       p_resid, width = 10, height = 3.6, dpi = 300)

p_fit <- ggplot(plot_df, aes(mu, y)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0, color = "orange", linewidth = 1) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "拟合优度散点：观测 vs 期望",
       x = "期望 μ", y = "观测 y") +
  theme_bw()
ggsave(file.path(output_dir, "fit_scatter_mu_vs_y.png"),
       p_fit, width = 4.6, height = 4.6, dpi = 300)

# 5) 可选：QQ 图（粗查负二项残差是否近似正态）
if (requireNamespace("ggqqplot", quietly = TRUE)) {
  # 若没有 ggqqplot 就跳过；保持脚本健壮
}


# ============================
# 敏感性分析：LHS + PRCC
# ============================
req_pkgs <- c("lhs","ppcor")
for (pk in req_pkgs) if (!requireNamespace(pk, quietly = TRUE)) install.packages(pk)
library(lhs); library(ppcor)

set.seed(2025)

# ----- 1) 采样范围（可按需微调） -----
# β0 从 Excel 给的 R0 区间推导；若不可用，就以拟合值 ±20% 为范围
if (exists("R0_range") && length(R0_range)==2 && all(is.finite(R0_range))) {
  beta0_lo <- R0_range[1] * pars_best$gamma
  beta0_hi <- R0_range[2] * pars_best$gamma
  if (!is.finite(beta0_lo) || !is.finite(beta0_hi) || beta0_lo<=0) {
    beta0_lo <- 0.8 * pars_best$beta0; beta0_hi <- 1.2 * pars_best$beta0
  }
} else {
  beta0_lo <- 0.8 * pars_best$beta0
  beta0_hi <- 1.2 * pars_best$beta0
}

sigma_lo <- 1/7   ; sigma_hi <- 2.0      # 潜伏期 0.5–7 天（周率）：[~1/7, 2]
gamma_lo <- 1/21  ; gamma_hi <- 2.0      # 感染期 0.5–21 天（周率）
I0_lo    <- max(1, 0.2 * init_state["I"])
I0_hi    <- max(I0_lo+1, 5.0 * init_state["I"])

# ----- 2) LHS 抽样 -----
n_mc <- 1000  # 抽样规模
U <- randomLHS(n = n_mc, k = 4)
LHS <- data.frame(
  beta0 = beta0_lo + U[,1] * (beta0_hi - beta0_lo),
  sigma = sigma_lo + U[,2] * (sigma_hi - sigma_lo),
  gamma = gamma_lo + U[,3] * (gamma_hi - gamma_lo),
  I0    = I0_lo + U[,4] * (I0_hi - I0_lo)
)

# η_t 固定为本次拟合得到的季节序列（只对 β0/σ/γ/I0 做敏感性）
eta_t <- plogis(as.numeric(best["logit_eta0"] + best["b1"]*F1[,1] + best["b2"]*F1[,2]))
t_idx_vec <- match(weeks, weeks)  # 1:T

# ----- 3) 评估函数 -----
eval_once <- function(b0, sg, gm, I0){
  # 初始状态：只替换 I0/E0；其余用你脚本中的 N_total
  E0 <- I0; R0 <- 0
  S0 <- max(1, N_total - I0 - E0 - R0)
  init_tmp <- c(S=S0, E=E0, I=I0, R=R0)
  
  # 模拟 I_inc
  sim <- simulate_weekly_series(init_tmp, beta_series = b0 * beta_shape,
                                sigma = sg, gamma = gm)
  # 对齐周序
  sim <- merge(data.frame(week_id = weeks), sim, by="week_id", all.x=TRUE)
  I_inc <- pmax(0, sim$I_inc)
  
  # 期望 μ
  mu <- pmax(1e-12, eta_t[t_idx_vec] * I_inc)
  
  # 指标
  final_frac   <- sum(I_inc, na.rm=TRUE) / N_total
  peak_I       <- max(I_inc, na.rm=TRUE)
  t_peak       <- which.max(I_inc)
  rmse         <- sqrt(mean((surv$y - mu)^2, na.rm=TRUE))
  total_mu     <- sum(mu, na.rm=TRUE)
  
  c(final_frac=final_frac, peak_I=peak_I, t_peak=t_peak, rmse=rmse, total_mu=total_mu)
}

# ----- 4) 批量计算 -----
mc_out <- t(mapply(eval_once, LHS$beta0, LHS$sigma, LHS$gamma, LHS$I0))
mc_df  <- cbind(LHS, as.data.frame(mc_out))
mc_df  <- as.data.table(mc_df)

# 保存结果
fwrite(mc_df, file.path(output_dir, "mc_sensitivity_results.csv"))

# ----- 5) PRCC（偏斯皮尔曼相关） -----
# 对每个输出 Y，计算控制其他自变量的偏相关（秩相关）
rank_df <- as.data.frame(apply(mc_df, 2, rank, ties.method = "average"))
Xvars <- c("beta0","sigma","gamma","I0")

calc_prcc <- function(yname){
  pp <- ppcor::pcor(rank_df[, c(Xvars, yname)], method = "spearman")
  data.frame(
    output = yname,
    var    = rownames(pp$estimate)[1:length(Xvars)],
    prcc   = as.numeric(pp$estimate[1:length(Xvars), ncol(pp$estimate)]),
    pval   = as.numeric(pp$p.value[1:length(Xvars), ncol(pp$p.value)])
  )
}

prcc_tbl <- do.call(rbind, lapply(c("final_frac","peak_I","t_peak","rmse","total_mu"), calc_prcc))
fwrite(as.data.table(prcc_tbl), file.path(output_dir, "mc_prcc_table.csv"))

# ----- 6) 可视化 -----

# 6.1 β×γ 对最终感染比例的热图（把 σ 固化在采样自身，不额外分层）
cut_beta  <- cut(mc_df$beta0, breaks = quantile(mc_df$beta0, probs = seq(0,1,0.2), na.rm=TRUE), include.lowest = TRUE)
cut_gamma <- cut(mc_df$gamma, breaks = quantile(mc_df$gamma, probs = seq(0,1,0.2), na.rm=TRUE), include.lowest = TRUE)
heat_df   <- mc_df[, .(final_frac = median(final_frac, na.rm=TRUE)), by=.(cut_beta, cut_gamma)]
p_heat <- ggplot(heat_df, aes(cut_beta, cut_gamma, fill = final_frac)) +
  geom_tile() +
  scale_fill_viridis_c(name = "最终感染比例", option = "C") +
  labs(x = "β0分组（分位数）", y = "γ分组（分位数）",
       title = "β0 与 γ 对最终感染比例的影响（分位数热图，σ/I0 已随机）") +
  theme_bw()
ggsave(file.path(output_dir, "mc_beta_gamma_heatmap.png"), p_heat, width = 7.6, height = 6, dpi = 300)

# 6.2 β 分组箱线图（最终感染比例）
beta_grp <- cut(mc_df$beta0, breaks = quantile(mc_df$beta0, probs = c(0,1/3,2/3,1), na.rm=TRUE), include.lowest = TRUE,
                labels = c("低β","中β","高β"))
p_box <- ggplot(data.frame(beta_grp, final_frac = mc_df$final_frac),
                aes(beta_grp, final_frac, fill = beta_grp)) +
  geom_boxplot(alpha = .8, outlier_alpha = .4) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "β0分组", y = "最终感染比例", title = "不同 β0 水平的最终感染比例") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(output_dir, "mc_beta_groups_finalfrac.png"), p_box, width = 6.8, height = 5.4, dpi = 300)

# 6.3 β 与最终规模散点 + LOESS
p_sc <- ggplot(mc_df, aes(beta0, final_frac)) +
  geom_point(alpha = .35, size = 1) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1, color = "steelblue") +
  labs(x = expression(beta[0]), y = "最终感染比例", title = expression(beta[0]*" 与最终感染比例")) +
  theme_bw()
ggsave(file.path(output_dir, "mc_beta_vs_finalfrac.png"), p_sc, width = 7.2, height = 5.0, dpi = 300)

message("✅ 敏感性分析完成：",
        "\n- mc_sensitivity_results.csv",
        "\n- mc_prcc_table.csv",
        "\n- mc_beta_gamma_heatmap.png",
        "\n- mc_beta_groups_finalfrac.png",
        "\n- mc_beta_vs_finalfrac.png")
