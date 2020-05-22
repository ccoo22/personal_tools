#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: metabolite_replace_missing.r -i <file>  -o <file> 
Options:
   -i, --input <file>         代谢物输入文件矩阵
   -o, --output <file>        输出文件" -> doc

opts                     <- docopt(doc, version = 'Program : 代谢物缺失值填充 v1.0 \n          甘斌 129\n')
input       <- opts$input
group_file  <- opts$group_file
output      <- opts$output
 
message("start metabolite_replace_missing.r")

data_raw <- read.table(input, head = T, check.names = F, stringsAsFactors =F, sep = "\t")

# 获取非0的最小值的一半
min_value <- min( data_raw[, 3:ncol(data_raw)][data_raw[, 3:ncol(data_raw)] > 0 ], na.rm = TRUE ) / 2 

# 填充NA
data_raw[, 3:ncol(data_raw)][is.na(data_raw[, 3:ncol(data_raw)])] <- min_value

# 速出
write.table(data_raw, output, sep = "\t", quote = F, row.names = F)

message("finish metabolite_replace_missing.r")