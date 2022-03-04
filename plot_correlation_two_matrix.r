#!/home/genesky/software/r/4.0.5/bin/Rscript

library(docopt)

"Usage: two_matrix_correlation.r  --matrix1 <file> --matrix2 <file> -s <file> -o <dir> [--cor_method <string> --width <numeric> --height <numeric> --pvalue_cutoff <numeric> --title <string> --tl_cex <numeric> --number_cex <numeric> --pch_cex <numeric> --cl_cex <numeric> --shape <string> --order_item <string> --item_color <string> --show_coef_value <color> --rlib <dir>]

Options:
	--matrix1 <file>                输入矩阵1，每行是一个特征，每列是一个样本
	--matrix2 <file>                输入矩阵2，每行是一个特征，每列是一个样本
	-s, --samplefile <file>         需要分析的样本名文件，一列样本名，无表头
	-o, --output <dir>              输出文件目录
	--cor_method <string>           相关性计算方法，pearson/spearman [default: spearman]
	--shape <string>                每一个cell的形状。circle, square, ellipse, number, shade, color, pie [default: circle]
	--show_coef_value <color>       在矩阵中显示相关性值，同时设置颜色 [default: white]
	--order_item <string>           对特征排序，排序方法 original AOE FPC hclust alphabet [default: original]
	--item_color <string>           特征名称颜色 [default: black]
	--pvalue_cutoff <numeric>       pvalue阈值，大于该阈值的cell，会标注特殊字符 [default: 0.05]
	--tl_cex <numeric>              特征名称字体大小 [default: 0.5]
	--number_cex <numeric>          相关性值字体大小 [default: 0.5]
	--pch_cex <numeric>             特殊字符大小 [default: 0.8]
	--cl_cex <numeric>              图例字体大小 [default: 0.5]
	--width <numeric>               pdf宽度 [default: 10]
	--height <numeric>              pdf高度 [default: 20]
	--title <string>                标题名称 [default: NA]
	--rlib <dir>                    R包路径 [default: /home/genesky/software/r/4.0.5/lib64/R/library]" -> doc

opts   <- docopt(doc, version='lhj 272，两个矩阵之间计算相关性\n')

matrix1         <- opts$matrix1
matrix2         <- opts$matrix2
samplefile      <- opts$samplefile
output          <- opts$output
cor_method      <- opts$cor_method
shape           <- opts$shape
show_coef_value <- opts$show_coef_value
order_item      <- opts$order_item
item_color      <- opts$item_color
pvalue_cutoff   <- as.numeric(opts$pvalue_cutoff)
tl_cex          <- as.numeric(opts$tl_cex)
number_cex      <- as.numeric(opts$number_cex)
pch_cex         <- as.numeric(opts$pch_cex)
cl_cex          <- as.numeric(opts$cl_cex)
width           <- as.numeric(opts$width)
height          <- as.numeric(opts$height)
title           <- opts$title
rlib            <- opts$rlib


### 测试数据
# matrix1     <- "/home/lhj/research/r_script/correlationship/species.taxon.Abundance.txt"
# matrix2     <- "/home/lhj/research/r_script/correlationship/env.txt"
# samplefile  <- "/home/lhj/research/r_script/correlationship/sample.txt"
# cor_method  <- "spearman"
# output      <- "/home/lhj/research/r_script/correlationship/out/"
# shape       <- "circle"
# show_coef_value  <- "white"
# order_item       <- "original"
# pvalue_cutoff    <- 0.05
# tl_cex           <- 0.5
# number_cex       <- 0.5
# pch_cex          <- 0.8
# cl_cex           <- 0.5
# title            <- NA
# width       <- 10
# height      <- 50
# rlib        <- "/home/genesky/software/r/4.0.5/lib64/R/library"


.libPaths(rlib)
library(corrplot)

data1   = read.table(matrix1, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "", fill = T)
data2   = read.table(matrix2, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "", fill = T)
samples = read.table(samplefile, header = F, sep = "\t", check.name = F, comment.char = "", quote= "", fill = T)

lost_samples1 = samples[!samples[, 1] %in% colnames(data1), 1]
if (length(lost_samples1) > 0){
    message("[Error] 部分需要分析的样本在 matrix1 文件中没有找到 : ", lost_samples1)
    q()
}
lost_samples2 = samples[!samples[, 1] %in% colnames(data2), 1]
if (length(lost_samples2) > 0){
    message("[Error] 部分需要分析的样本在 matrix2 文件中没有找到 : ", lost_samples2)
    q()
}

data_need1 = data1[, samples[, 1], drop = F]
data_need2 = data2[, samples[, 1], drop = F]

res = t(sapply(1:nrow(data_need1), function(x) {sapply(1:nrow(data_need2), function(y) {res = cor.test(as.matrix(data_need1)[x,], as.matrix(data_need2)[y,], method = cor_method);c(res[[3]],res[[4]])})}))
cols_p = seq(1, ncol(res)-1, 2)
cols_r = seq(2, ncol(res), 2)
cor_p = res[, cols_p, drop = F]
cor_r = res[, cols_r, drop = F]

rownames(cor_r) = rownames(data1)
colnames(cor_r) = rownames(data2)
rownames(cor_p) = rownames(data1)
colnames(cor_p) = rownames(data2)
write.csv(cor_r, paste0(output,"/Correlation_r.csv"))
write.csv(cor_p, paste0(output,"/Correlation_p.csv"))
result = data.frame(Matrix1 = rep(rownames(cor_r), ncol(cor_r)), Matrix2 = rep(colnames(cor_r), each = nrow(cor_r)), Cor = as.vector(cor_r), Pvalue = as.vector(cor_p))
write.csv(result, paste0(output,"/Correlation_All.csv"), row.names = F)

corr0 = signif(cor_r, 2)

# color <- colorRampPalette(c("blue", "white", "red"))
color <- colorRampPalette(c("blue", "red"))

pdf(paste(output, "Correlation.pdf", sep = "/"), width = width, height = height)
corrplot(corr0, 
        method = shape, 
        title = title, 
        mar = c(0, 0, 1, 0), 
        addCoef.col = show_coef_value, 
        order = order_item, 
        tl.cex = tl_cex, 
        tl.col = item_color, 
        number.cex = number_cex,
        pch.cex = pch_cex,
        pch.col = "black",
        # cl.lim = c(-0.5, 0.5),
        cl.cex = cl_cex,
        p.mat = cor_p, 
        sig.level = pvalue_cutoff, 
        col = color(100)
        )
dev.off()
