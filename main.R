# 主脚本（V3）：beta_shape 在预处理阶段生成；拟合阶段直接使用
# 使用：setwd 到本目录 -> source("main.R", encoding="UTF-8")
setwd("E:/a第七学期/ai比赛/调试")
output_dir <- file.path(getwd(), "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

data_path  <- file.path(getwd(), "FluNet_格式化_104周.csv")
meta_path  <- file.path(getwd(), "region_meta.csv")   # 可选

source(file.path(getwd(), "data_preprocessing_holtwinters.R"), encoding = "UTF-8")
source(file.path(getwd(), "seasonal_seir_model_fit.R"),         encoding = "UTF-8")

message("✔ 全流程完成。输出：", normalizePath(output_dir))
