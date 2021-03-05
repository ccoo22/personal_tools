#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: select_gene_based_pamr.r --exp <file> --case_sample <string> --control_sample <string> --output_dir <dir> [ --threshold <numeric>]
Options:
    --exp <file>                 表达量矩阵，每一行是一个基因，每一列是一个样本
    --case_sample <string>       case组样本列表，样本之间用逗号分隔
    --control_sample <string>    control组样本列表，样本之间用逗号分隔
    --output_dir <dir>           结果输出目录
    --threshold <numeric>        基因筛选阈值 默认情况下，自动选择
 " -> doc

opts           <- docopt(doc, version = '根据分组，使用pamr包筛选基因 \n        甘斌 129\n')
exp            <- opts$exp
case_sample    <- opts$case_sample
control_sample <- opts$control_sample
threshold      <- as.numeric(opts$threshold)
output_dir     <- opts$output_dir

# exp = './dmp_sig.txt'
# case_sample = '19,26,27,33'
# control_sample = '43,123,119,16'

library(pamr)


data_exp        = read.table(exp, header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names=F)
case_samples    = unlist(strsplit(case_sample, ',' ))
control_samples = unlist(strsplit(control_sample, ',' ))

# 输入数据准备
x         <- as.matrix(data_exp[, c(case_samples, control_samples)])
y         <- c(rep('case', length(case_samples)), rep('control', length(control_samples)))
geneid    <- rownames(x)
mydata <- list(x = x, y=factor(y, levels=c('case', 'control')), geneid = geneid)


# train
mytrain<- pamr.train(mydata)

# 最多选择10个基因
mythreshold_10 = mytrain$threshold[which(mytrain$nonzero <= 10)[1]]  
mythreshold_100 = mytrain$threshold[which(mytrain$nonzero <= 100)[1]]  
 
# gene list
mylistgenes_10 <- pamr.listgenes(mytrain, mydata, threshold=mythreshold_10)
mylistgenes_100 <- pamr.listgenes(mytrain, mydata, threshold=mythreshold_100)
write.table(mylistgenes_100, paste0(output_dir, "/top100.genelist.txt"), quote=F, row.names=F, col.names=T)

# gene plot
pdf(paste0(output_dir, "/top10.pdf"))
pamr.geneplot(mytrain, mydata, threshold=mythreshold_10)
dev.off()
