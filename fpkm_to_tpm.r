#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: fpkm_to_tpm.r -i <file> -o <file> 

Options:
    -i, --input <file>              样本fpkm矩阵，每一行为一个基因，每一列为一个样本，第一列为基因名，第一行为样本名。
    -o, --output <file>             输出文件
    " -> doc

opts   <- docopt(doc, version='甘斌，FPKM 转 TPM\n')
input               <- opts$input
output              <- opts$output
 
###################################################################### 主程序
fpkmToTpm <- function(fpkm)
{
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


message("read input")
expMatrix<-read.table(input,header = T,row.names = 1, check.names=F)

message("convert to tpm")
tpms <- apply(expMatrix, 2, fpkmToTpm)
tpms = cbind(gene = rownames(tpms), tpms)

message("output")
write.table(tpms, output, quote = FALSE, row.names = FALSE, sep = '\t')

