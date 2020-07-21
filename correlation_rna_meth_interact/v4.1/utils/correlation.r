#!/usr/bin/env Rscript

library(docopt)
"Usage: correlation.r --rna_file <file> --meth_file <file> --output_file <file> [--correct_methon <string>]
Options:
   -r, --rna_file <file>          表达量数据，每一行为基因、转录本，每一列为样本,含表头，前4列分别是id/chr/start/end
   -m, --meth_file <file>         甲基化数据，每一行为甲基化，每一列为样本，含表头，前4列分别是id/chr/start/end
   -o, --output_file <file>       结果输出目录 
   -c, --correct_methon <string>  计算方法，支持pearson/spearman[default: pearson]" -> doc

 
opts           <- docopt(doc)
rna_file       <- opts$rna_file
meth_file      <- opts$meth_file
output_file    <- opts$output_file
correct_methon <- opts$correct_methon

# (1) 数据读入
data_rna      <- read.table(rna_file, sep='\t', header = TRUE, check.names=FALSE, stringsAsFactors=FALSE)
data_meth     <- read.table(meth_file, sep='\t', header = TRUE, check.names=FALSE, stringsAsFactors=FALSE)
 
# (2) 计算
result = matrix(nrow = nrow(data_rna) * nrow(data_meth), ncol = 10);

row = 0
for(row_rna in 1:nrow(data_rna))
{
    for(row_meth in 1:nrow(data_meth))
    {   
        row = row + 1
        data <- cbind(rna=as.numeric(data_rna[row_rna, 5:ncol(data_rna)]), meth=as.numeric(data_meth[row_meth, 4:ncol(data_meth)]))
        data <- data.frame(data[complete.cases(data),])
        data = data[apply(data == "", 1, sum) == 0, ]  # 去掉空值
    
        nmiss = nrow(data)
        pvalue = NA
        estimate = NA
        if(nmiss > 1)
        {
            tmp <- cor.test(data[, 1], data[, 2], alternative = "two.sided", method = correct_methon, conf.level = 0.95)        
            pvalue <- tmp$p.value
            estimate <- tmp$estimate[[1]]
        }

        result[row, ] = c(as.matrix(data_rna[row_rna, 1:4])[1, ], as.matrix(data_meth[row_meth, 1:3])[1,], nmiss, pvalue, estimate )
    }
}
write.table(result, output_file , sep = "\t", quote = F, row.names = F ,col.names = F)
