# ============================
# 敏感性分析：LHS + PRCC
# ============================
req_pkgs <- c("lhs","ppcor")
for (pk in req_pkgs) if (!requireNamespace(pk, quietly = TRUE)) install.packages(pk)
library(lhs); library(ppcor)

set.seed(2025)

# ----- 1) 采样范围（可按需微调） -----
# β0 从 Excel 给的 R0 区间推导；若不可用，就以拟合值 ±20% 为范围
# ======【周尺度一致性设置】======
# 若你在 Excel/文献里用“天”为单位给了先验，就在这里转成“每周速率”
latent_days_range     <- c(2, 7)     # 潜伏期天数的合理区间（示例）
infectious_days_range <- c(2, 10)    # 传染期天数的合理区间（示例）
R0_range_days         <- if (exists("R0_range")) R0_range else c(1.2, 1.8)

# 把天数转为每周速率
sigma_lo <- 7 / max(latent_days_range)   # latent_days 大 → sigma 小
sigma_hi <- 7 / min(latent_days_range)
gamma_lo <- 7 / max(infectious_days_range)
gamma_hi <- 7 / min(infectious_days_range)

# 用“周速率的 γ”把 R0 区间换成 β0 区间
# 这里用 γ 的中值代表：也可用区间保守化处理（见下方注释）
gamma_mid <- 7 / mean(infectious_days_range)
beta0_from_R0 <- R0_range_days * gamma_mid   # β0 区间（每周）

# 兜底：若你希望以拟合值为中心 ±20%
beta0_lo_fallback <- 0.8 * pars_best$beta0
beta0_hi_fallback <- 1.2 * pars_best$beta0

# 最终 β0 区间（优先 R0 推导，失败用拟合±20%）
if (all(is.finite(beta0_from_R0))) {
  beta0_lo <- min(beta0_from_R0); beta0_hi <- max(beta0_from_R0)
} else {
  beta0_lo <- beta0_lo_fallback;  beta0_hi <- beta0_hi_fallback
}

# I0 区间：以当前拟合 I0 为中心按比例给范围（周粒度的观测）
I0_base <- max(1, init_state["I"])
I0_lo   <- max(1, 0.2 * I0_base)
I0_hi   <- max(I0_lo + 1, 5.0 * I0_base)

message(sprintf("[敏感性·单位检查] sigma_week ∈ [%.3f, %.3f]  (latent_days %d–%d)",
                sigma_lo, sigma_hi, min(latent_days_range), max(latent_days_range)))
message(sprintf("[敏感性·单位检查] gamma_week ∈ [%.3f, %.3f]  (infectious_days %d–%d)",
                gamma_lo, gamma_hi, min(infectious_days_range), max(infectious_days_range)))
message(sprintf("[敏感性·单位检查] beta0_week ∈ [%.3f, %.3f]  (from R0 %g–%g, gamma_mid=%.3f)",
                beta0_lo, beta0_hi, min(R0_range_days), max(R0_range_days), gamma_mid))
message(sprintf("[敏感性·单位检查] I0 ∈ [%.0f, %.0f] (基于拟合 I0=%.0f)", I0_lo, I0_hi, I0_base))

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
  E0 <- I0; 
  R0 <- I0
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
