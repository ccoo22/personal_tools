#!/home/genesky/software/r/3.6.1/bin/Rscript
library(docopt)
"Usage: split_data.r --rna_file <file> --meth_file <file> --sample_file <file> --interact <file> --output_dir <dir>
Options:
   -r, --rna_file <file>         表达量数据，每一行为基因、转录本，每一列为样本，含表头
   -m, --meth_file <file>        甲基化数据，每一行为甲基化，每一列为样本，含表头
   -s, --sample_file <file>      样本文件，一列，务必与表达量、甲基化一致,没有表头
   -i, --interact <file>         经过bedtools interact 取交集后的文件，其中第4列为表达量的第一列信息，第8列为甲基化的第一列信息
   -o, --output_dir <dir>        结果输出目录" -> doc

opts                <- docopt(doc, version = '根据区域交集信息，对数据进行拆分。 v1.0 \n')
rna_file     <- opts$rna_file
meth_file    <- opts$meth_file
sample_file  <- opts$sample_file
interact     <- opts$interact
output_dir   <- opts$output_dir


# (1) 数据读入
data_rna      <- read.table(rna_file, sep='\t', header = TRUE, check.names=FALSE, stringsAsFactors=FALSE)
data_meth     <- read.table(meth_file, sep='\t', header = TRUE, check.names=FALSE, stringsAsFactors=FALSE)
data_sample   <- read.table(sample_file, sep='\t', header = FALSE, check.names=FALSE, stringsAsFactors=FALSE)
data_interact <- read.table(interact, sep='\t', header = FALSE, check.names=FALSE, stringsAsFactors=FALSE)
rownames(data_rna) = data_rna[, 1]
rownames(data_meth) = data_meth[, 1]
# 样本
samples <- as.character(data_sample[, 1])

# (2) 检查样本是否有缺失
lost_sample_rna <- samples[which(!samples %in% colnames(data_rna))]
lost_sample_meth <- samples[which(!samples %in% colnames(data_meth))]
if(length(lost_sample_rna) > 0) stop('检测到部分样本在rna数据中缺失 : ', lost_sample_rna)
if(length(lost_sample_meth) > 0) stop('检测到部分样本在meth数据中缺失 : ', lost_sample_meth)

# (3) 确认存在交集的基因、转录本
for(gene in unique(data_interact[, 4]))
{  
   # 交集甲基化
   meth_site = data_interact[, 8][data_interact[, 4] == gene] # 区域交集的甲基化位点
   
   # 拆分数据
   split_gene = cbind(data_rna[gene, 1:4, drop=FALSE], data_rna[gene, samples])
   split_meth = cbind(data_meth[meth_site, 1:3, drop=FALSE], data_meth[meth_site, samples])
   
   # 输出
   split_gene_file = paste0(output_dir, '/', gene, '.gene.txt') 
   split_meth_file = paste0(output_dir, '/', gene, '.meth.txt') 
   write.table(split_gene, split_gene_file, sep = "\t", quote = F, row.names = F, col.names = TRUE )
   write.table(split_meth, split_meth_file, sep = "\t", quote = F, row.names = F, col.names = TRUE )
}