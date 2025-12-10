# subset_realdata.R
set.seed(2025)

## ---- 1. 读入原始数据 ----
cov_file <- "input/covar_stool.txt"
abd_file <- "input/pheno_stool_norrf_count.txt"

cov <- read.table(cov_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
abd <- read.table(abd_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# 确定样本 ID 列名
id_col_cov <- colnames(cov)[1]
id_col_abd <- colnames(abd)[1]

## ---- 2. 随机选取 100 个样本 ----
common_ids <- intersect(cov[[id_col_cov]], abd[[id_col_abd]])
if (length(common_ids) < 100) stop("样本不足 100 个！")
ids_keep <- sample(common_ids, size = 100, replace = FALSE)

## ---- 3. 提取 cov 中需要的列 ----
# 尝试匹配 PC1–PC5 / PC01–PC05 的列名
pc_cols <- grep("^PC0?[1-5]$", colnames(cov), value = TRUE)
sel_cov_cols <- c(id_col_cov, "SAMPLE_COLLECTION_AGE_MONTHS_DEV", "SEX", pc_cols)
sel_cov_cols <- sel_cov_cols[sel_cov_cols %in% colnames(cov)]

cov_sub <- cov[cov[[id_col_cov]] %in% ids_keep, sel_cov_cols]

## ---- 4. 提取 abd 中前 3 个 feature ----
abd_features <- setdiff(colnames(abd), id_col_abd)
sel_abd_cols <- c(id_col_abd, abd_features[1:3])

abd_sub <- abd[abd[[id_col_abd]] %in% ids_keep, sel_abd_cols]

## ---- 5. 按样本顺序对齐 ----
cov_sub <- cov_sub[match(ids_keep, cov_sub[[id_col_cov]]), ]
abd_sub <- abd_sub[match(ids_keep, abd_sub[[id_col_abd]]), ]

## ---- 6. 输出 ----
write.table(cov_sub, "cov_small.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(abd_sub, "abd_small.txt", quote = FALSE, sep = "\t", row.names = FALSE)

## ---- 7. 生成 keep 文件给 PLINK ----
keep <- data.frame(FID = 0, IID = ids_keep)
write.table(keep, "keep_samples.txt", quote = FALSE, sep = " ",
            row.names = FALSE, col.names = FALSE)

cat("✅ 已生成文件：cov_small.tsv, abd_small.tsv, keep_samples.txt\n")
