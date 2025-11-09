# Seasonal SEIR V3（beta_shape 在预处理阶段生成）

**改动点：**
- 在 `data_preprocessing_holtwinters.R` 中，用 **STL（失败回退傅里叶）** 从阳性率提取季节分量；
- 做 **5周中位数平滑 + Winsor(0.6,1.6) + 线性缩幅(λ=0.45)**；
- 将 `beta_shape` 连同 `surv` 一起写入 `preproc.rds`；
- `seasonal_seir_model_fit.R` **直接读取 `beta_shape`**，不再重复计算；
- 继续使用 **相对时间索引**、**lsoda**、**参数护栏** 和 **mu 不乘 w**。

**用法：**
1. 将 `FluNet_格式化_104周.csv` 放在本目录；
2. 在 R 中执行：`source("main.R", encoding="UTF-8")`；
3. 输出：`output/fit_observed_vs_mu.png`, `output/fitted_params.csv` 等。

**观测量 I：** CSV 中 `A(H3N2)(%)` 的括号前计数（每周 H3N2 确证数）。
