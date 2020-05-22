#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: metabolite_normalized.r -i <file>  -o <file> 
Options:
   -i, --input <file>         代谢物输入文件矩阵
   -o, --output <file>        输出文件" -> doc

opts                     <- docopt(doc, version = 'Program : 代谢物样本水平归一化 v1.0 \n          甘斌 129\n')
input       <- opts$input
group_file  <- opts$group_file
output      <- opts$output

message("start metabolite_normalized.r")
data_raw <- read.table(input, head = T, check.names = F, stringsAsFactors =F, sep = "\t")

norm_sum <- function(x)
{
   1000000 * x/sum(x)
}

# 矫正
data_norm <- apply(data_raw[, 3:ncol(data_raw)], 2, norm_sum)

# 补充名称
data_norm <- data.frame(data_raw[,1:2], data_norm)

# 输出
write.table(data_norm, output, sep = "\t", quote = F, row.names = F)

message("finish metabolite_normalized.r")
