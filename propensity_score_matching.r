#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: propensity_score_matching.r  -i <file> -o <prefix>  [--rlib <string>]

Options:
    -i, --input <file>              输入文件，一样一个样本，第一列样本ID, 第二列分组信息，必须是1/0，分别代表case/control，第三列及之后的数据为临床表型数据。
                                    不要有缺失值，当然，流程会自动存在缺失值的样本。
                                    第一行表头，不要存在空格、“-”等特殊符号，尽可能是 字母、数字、下划线构成。
    -o, --output <prefix>           输出文件前缀,示例：./a.pdf
    --rlib <string>                 r包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，PSM分析，挑选合适的配对样本集，减少其他环境因素的影响，再说的白话一点，就是挑选样本，使当前case/control条件下，其他的临床数据不要有差异\n')
input              <- opts$input
output             <- opts$output
rlib               <- opts$rlib
 
.libPaths(rlib)
library(MatchIt)

message("读入数据")
data = read.table(input, head = T, row.names = 1, sep = "\t", stringsAsFactors=F)
if(sum(grepl(" ", colnames(data)) > 0 or sum(grepl("-", colnames(data)))
{
    message("输入文件的表头含有空格或-字符，请避免使用！\n");
    q()
}
data_clean = data[complete.cases(data), ]


# 构建模型
yname = colnames(data_clean)[1]
envs = colnames(data_clean)[2:ncol(data_clean)]
formula <- as.formula(paste(yname,' ~ ', paste(envs, collapse = ' + ') ) )

m.out = matchit(formula, data = data_clean, method ="nearest", ratio =1)

summary (m.out)

# PSM 得分分布
psm_pdf = paste0(output, ".psm.pdf")
pdf(psm_pdf)
plot(m.out, type = "jitter")
user.prompt()
dev.off()

# match之后，每个特征的点状图，x/y 分别是control样本，case样本。配对的。
matched_variable_pdf = paste0(output, ".matched_variable.pdf")
pdf(matched_variable_pdf)
plot(m.out)
user.prompt()
dev.off()

# 配对情况
match_pair = m.out$match.matrix
match_pair = cbind(case_sample = rownames(match_pair), control_samples = match_pair[,1])

match_pair_txt = paste0(output, ".matched_pair_sample.txt")
write.table(match_pair, file = match_pair_txt, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 配对后的临床矩阵
match_result = match.data(m.out)
match_result = cbind(sample = row.names(match_result), match_result)

result_txt = paste0(output, ".matched_result.txt")
write.table(match_result, file = result_txt, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


message("[result] 配对后的样本矩阵：        ", result_txt)
message("[result] 配对样本对：             ", match_pair_txt)
message("[result] 配对后的ps值绘图：       ", psm_pdf)
message("[result] 配对后的每个特征的点状图：", matched_variable_pdf)
 