# 检测 -> 脚本输入
ARGS <- commandArgs(T)
if(length(ARGS) < 2){
  cat("
Program: kallisto_to_reads_count 
Version: v1.0
Contact: 129 甘斌
      
Usage:   Rscript stringdb.r ARGS1 ARGS2 ARGS3
      
Options:
      
      file_list         kallisto定量输出结果的文件列表。目前仅支持人。两列，第一列：样本名，第二列：文件路径，没有表头
      OUT               输出结果
      

      \n");
      q()
}

library(tximport)

inputfile  = normalizePath(ARGS[1])
output     = normalizePath(ARGS[2])

# 输入文件
file_table = read.table(inputfile, head = F, stringsAsFactors = F)
file_list  = file_table[[2]]
names(file_list) = file_table[[1]] # 样本名


# 转录本->genesymbol映射关系
library("EnsDb.Hsapiens.v86")
library("ensembldb")
db <- EnsDb.Hsapiens.v86
tx <- transcripts(db, return.type = "DataFrame")
tx <- tx [, c("tx_id", "gene_id")]

gene      <- genes(db, return.type = "DataFrame")
gene_type <- gene[, c("symbol", "gene_biotype")]
gene_type <- gene_type[!duplicated(gene_type$symbol) ,] # 安基因名，去重复
gene      <- gene[, c("gene_id", "symbol")]

tx2genes <- merge(tx, gene, by = "gene_id")
tx2genes <- tx2genes[, 2:3]

# 读入表达量,并在基因水平汇总
txi <- tximport(files = file_list, type = "kallisto", tx2gene = tx2genes, ignoreTxVersion = T) # 必须加上ignoreTxVersion参数，防止转录本版本号导致无法匹配

# 基因biotype注释
gene_anno = data.frame(symbol = row.names(txi$counts))
gene_anno = merge(gene_anno, gene_type, by = 'symbol', all.x = T, sort = F) # 获取基因的类型

# 输出
count = cbind(Gene = row.names(txi$counts), txi$counts)
count_anno = cbind(count, gene_anno$gene_biotype)
write.table(count, file = output, sep = "\t", quote = F, row.names = FALSE)
write.table(count_anno, file = paste0(output, '.anno_gene_biotype.txt'), sep = "\t", quote = F, row.names = FALSE)
