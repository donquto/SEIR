# 数据读取 + 预处理 + 生成 beta_shape（STL优先，失败回退傅里叶），并保存到 preproc.rds

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(stringr)
  library(ggplot2); library(zoo)
})

if (!exists("output_dir")) output_dir <- file.path(getwd(), "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!exists("data_path")) data_path <- file.path(getwd(), "FluNet_格式化_104周.csv")
if (!exists("meta_path")) meta_path <- file.path(getwd(), "region_meta.csv")

parse_npct <- function(x){
  n_val  <- suppressWarnings(as.numeric(str_extract(x, "^[0-9]+")))
  pct_tx <- str_replace(str_extract(x, "(?<=\\().+?(?=%\\))"), ",", ".")
  pct_val<- suppressWarnings(as.numeric(pct_tx))
  list(n = n_val, pct = pct_val/100)
}
cn_week_to_int <- function(x) as.integer(str_extract(x, "[0-9]+"))

# ---- 读取 ----
raw_bytes <- read_file_raw(data_path)
enc <- readr::guess_encoding(raw_bytes)$encoding[1]
df <- read_csv(data_path, locale = locale(encoding = enc), show_col_types = FALSE) |> as.data.frame()
names(df) <- sub("^\ufeff", "", names(df), useBytes = TRUE)
names(df) <- gsub("\\s+", "", names(df), perl = TRUE)
names(df) <- chartr("（）", "()", names(df))

col_year   <- grep("^年", names(df), value = TRUE)[1]
col_week   <- grep("^周次", names(df), value = TRUE)[1]
col_region <- grep("^地区", names(df), value = TRUE)[1]
col_tests  <- grep("^检测数", names(df), value = TRUE)[1]
col_pos    <- grep("^阳性数\\(%\\)", names(df), value = TRUE)[1]
col_A      <- grep("^A型\\(%\\)", names(df), value = TRUE)[1]
col_H1     <- grep("^A\\(H1N1\\)pdm09\\(%\\)", names(df), value = TRUE)[1]
col_H3     <- grep("^A\\(H3N2\\)\\(%\\)", names(df), value = TRUE)[1]
if (any(is.na(c(col_year, col_week, col_region, col_tests, col_pos, col_A, col_H1, col_H3)))) {
  stop("监测数据缺少必要列，请检查文件格式。")
}

pos_np <- parse_npct(df[[col_pos]])
A_np   <- parse_npct(df[[col_A]])
H1_np  <- parse_npct(df[[col_H1]])
H3_np  <- parse_npct(df[[col_H3]])

df_processed <- df |>
  mutate(
    year   = as.integer(.data[[col_year]]),
    week   = cn_week_to_int(.data[[col_week]]),
    region = .data[[col_region]],
    tests  = suppressWarnings(as.numeric(.data[[col_tests]])),
    pos_n  = pos_np$n,  pos_rate = pos_np$pct,
    A_n    = A_np$n,    A_rate   = A_np$pct,
    H1_n   = H1_np$n,   H1_rate  = H1_np$pct,
    H3_n   = H3_np$n,   H3_rate  = H3_np$pct
  ) |>
  select(year, week, region, tests,
         pos_n, pos_rate, A_n, A_rate, H1_n, H1_rate, H3_n, H3_rate) |>
  as.data.frame()

df_processed <- df_processed[order(df_processed$year, df_processed$week, df_processed$region), ]
weeks_unique <- unique(df_processed[, c("year","week")])
weeks_unique <- weeks_unique[order(weeks_unique$year, weeks_unique$week), ]
weeks_unique$week_id <- seq_len(nrow(weeks_unique))
surv <- merge(df_processed, weeks_unique, by = c("year","week"), all.x = TRUE)
surv <- as.data.table(surv)
setorder(surv, week_id, region)

# 元数据（可选）
if (file.exists(meta_path)) {
  meta_raw <- read_file_raw(meta_path)
  enc_m <- guess_encoding(meta_raw)$encoding[1]
  meta_df <- read_csv(meta_path, locale = locale(encoding = enc_m), show_col_types = FALSE) |> as.data.frame()
  names(meta_df) <- sub("^\ufeff", "", names(meta_df), useBytes = TRUE)
  names(meta_df) <- gsub("\\s+", "", names(meta_df), perl = TRUE)
  names(meta_df) <- chartr("（）", "()", names(meta_df))
  col_reg_m <- grep("^region|^地区", names(meta_df), value = TRUE)[1]
  col_pop_m <- grep("^pop|人口", names(meta_df), value = TRUE)[1]
  col_sent_m<- grep("^sent|哨点", names(meta_df), value = TRUE)[1]
  if (!any(is.na(c(col_reg_m, col_pop_m, col_sent_m)))) {
    meta_df <- meta_df[, c(col_reg_m, col_pop_m, col_sent_m)]
    names(meta_df) <- c("region","pop","sent")
    meta_df$pop <- suppressWarnings(as.numeric(meta_df$pop))
    meta_df$sent<- suppressWarnings(as.numeric(meta_df$sent))
    surv <- merge(surv, meta_df, by = "region", all.x = TRUE)
  } else {
    surv[, `:=`(pop = NA_real_, sent = NA_real_)]
  }
} else {
  surv[, `:=`(pop = NA_real_, sent = NA_real_)]
}

# 权重（仅用于样本/方差）
surv[, test_rate := ifelse(is.finite(pop), tests / pmax(pop, 1), tests)]
surv[, E := test_rate / median(test_rate[is.finite(test_rate)], na.rm = TRUE), by = region]
surv[, P := pmin(1, pmax(0, H3_rate))]
surv[, w := pmax(1e-3, E * P)]
surv[, omega := pmax(1, tests)]
surv[, y := pmax(0, H3_n)]

# ===== 用 STL/傅里叶 提取季节形状（周频=52），并做振幅控制 =====
x_raw <- as.numeric(surv$pos_rate)
x_fix <- na.approx(x_raw, maxgap = 8, na.rm = FALSE)
x_fix <- na.locf(x_fix, na.rm = FALSE); x_fix <- na.locf(x_fix, fromLast = TRUE)
med   <- median(x_fix[is.finite(x_fix)], na.rm = TRUE)
x_fix[!is.finite(x_fix)] <- med
x_fix <- pmax(1e-6, x_fix)

if (length(x_fix) > 104) x_fix <- tail(x_fix, 104)
if (length(x_fix) < 104) stop("可用周数不足 104，当前=", length(x_fix))
ts_rate <- ts(x_fix, start = c(1,1), frequency = 52)

stl_season_or_fourier <- function(ts_rate){
  season <- tryCatch({
    stl_fit <- stl(ts_rate, s.window = "periodic", robust = TRUE)
    as.numeric(stl_fit$time.series[, "seasonal"])
  }, error = function(e) rep(NA_real_, length(ts_rate)))

  if (all(!is.finite(season))) {
    t <- seq_along(ts_rate); w <- 2*pi/52
    X <- cbind(s1 = sin(w*t),  c1 = cos(w*t),
               s2 = sin(2*w*t), c2 = cos(2*w*t))
    fit <- lm(as.numeric(ts_rate) ~ X)
    season <- as.numeric(fitted(fit)) - mean(fitted(fit))
  }
  season_s <- zoo::rollmedian(season, k = 5, fill = "extend", align = "center")
  shape_raw <- season_s / mean(season_s, na.rm = TRUE)
  shape_clip <- pmin(2.0, pmax(0.5, shape_raw))  # 原来 0.6–1.6 → 放宽到 0.5–2.0
  lambda <- 0.70                                 # 原来 0.45 → 0.70（提高峰顶）
  beta_shape <- 1 + lambda * (shape_clip - 1)
  beta_shape
}

beta_shape <- stl_season_or_fourier(ts_rate)
cat("beta_shape 范围=", range(beta_shape), " 平均=", mean(beta_shape), "\n")

# 保存
model_total <- list(beta_shape = beta_shape)
saveRDS(list(surv = surv, model_total = model_total), file = file.path(output_dir, "preproc.rds"))
message("预处理完成：surv / beta_shape 已保存到 preproc.rds")
