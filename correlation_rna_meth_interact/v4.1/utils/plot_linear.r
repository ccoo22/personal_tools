#/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)
"Usage: plot_linear.r -i <file> -o <file> --rna <file> --meth <file> --sample <file>  [--top <integer>]

Options:
  -i <file>                   输入文件，应该输入文件 filter.correlation.result.txt。需要存在几个关键的表头： RNA_ID METH_ID Pvalue
  -o <pdf>                    输出pdf路径，例如： ./linear_plot.pdf
  --rna <file>                RNA表达量数据矩阵，第一列是RNA_ID, 后面的列可以是样本或注释
  --meth <file>               甲基化数据矩阵，第一列是METH_ID, 后面的列可以是样本或注释
  --sample <file>             样本文件，一列，务必与表达量、甲基化一致,没有表头
  --top <integer>             输入文件中，对前n个结果绘图 [default: 300]"-> doc

opts          <- docopt(doc)
input         <- opts$i
output        <- opts$o
rna           <- opts$rna
meth          <- opts$meth
sample        <- opts$sample
top           <- as.integer(opts$top)

# rna ='gene.xls'
# meth = 'methyl.xls'
# sample = 'samples.list'
# input = 'result/filter.correlation.result.txt'
message("load data")
data_rna      <- read.table(rna, sep='\t', header = TRUE, check.names=FALSE, stringsAsFactors=FALSE)
data_meth     <- read.table(meth, sep='\t', header = TRUE, check.names=FALSE, stringsAsFactors=FALSE)
data_sample   <- read.table(sample, sep='\t', header = FALSE, check.names=FALSE, stringsAsFactors=FALSE)
data_input    <- read.table(input, sep='\t', header = TRUE, check.names=FALSE, stringsAsFactors=FALSE)
rownames(data_rna) = data_rna[, 1]
rownames(data_meth) = data_meth[, 1]

data_input = data_input[order(data_input$Pvalue), ]  # Pvalue 排序
# 样本
samples <- as.character(data_sample[, 1])


# (2) 检查样本是否有缺失
message("check sample")

lost_sample_rna <- samples[which(!samples %in% colnames(data_rna))]
lost_sample_meth <- samples[which(!samples %in% colnames(data_meth))]
if(length(lost_sample_rna) > 0) stop('检测到部分样本在rna数据中缺失 : ', lost_sample_rna)
if(length(lost_sample_meth) > 0) stop('检测到部分样本在meth数据中缺失 : ', lost_sample_meth)

if(top > nrow(data_input)) top = nrow(data_input)

message("start plot")

pdf(output)
for(row in 1:top)
{
    rna_id = data_input[row, 'RNA_ID']
    meth_id = data_input[row, 'METH_ID']
    pvalue  = format(data_input[row, 'Pvalue'], scientific=TRUE, digit=4) # 科学计数，更美观，字符型数据
    
    data_plot = data.frame(rna = as.numeric(data_rna[rna_id, samples]), meth = as.numeric(data_meth[meth_id, samples]))
    data_plot = data_plot[complete.cases(data_plot), , drop = FALSE] # 去掉NA
    data_plot = data_plot[apply(data_plot == "", 1, sum) == 0, ]  # 去掉空值
    if(is.null(nrow(data_plot))) next #只有一个样本或一个位点，无法进行分析

    # 拟合
    x = data_plot$meth
    y = data_plot$rna
    lm_fit <- lm(y ~ x)
    plot(x = x, y = y, xlab = meth_id, ylab = rna_id , pch = 16, main = paste0(rna_id, " ", meth_id) ) # 黑色散点
    lines(x, fitted(lm_fit), col="red") # 拟合曲线

    text_formula <- paste("y = ", round(lm_fit$coefficients['x'], 5), "x + ", round(lm_fit$coefficients['(Intercept)'], 5) , sep="")
    text_rsquared <- paste("R-squared = ", round(summary(lm_fit)$r.squared, 4), sep = "")

    legend("topright", legend = c(text_formula, text_rsquared))
}

dev.off()
