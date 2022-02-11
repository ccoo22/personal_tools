#!/home/genesky/software/r/4.0.3/bin/Rscript

library(docopt)
"Usage: enricher.r -g <file> -o <file> -d <file> [--rlib <dir>]
Options:
        -g, --gene <file>            差异基因列表,没有表头，一列基因名称
        -o, --output <dir>           输出路径

        -d, --database <file>        个人创建的KEGG数据库， Gene_ID  列是基因名称， koid是对应的KO基因ID，Pathway列是基因所属的通路，多个通路用分号分隔，每个通路的格式为 '通路ID|通路描述'
                                     注意：
                                         （1）输入的基因列表必须能与Gene_ID对应上
                                         （2）必须保证database文件中，Gene_ID没有重复.
        --rlib <dir>                 R包路径 [default: /home/genesky/software/r/4.0.3/lib64/R/library]" -> doc

opts     <- docopt(doc, version = 'Program : 基于自定义的通路数据库，通路富集分析 \n      甘斌 129\n')

gene_file       <- opts$gene
output_dir      <- opts$output
database_file   <- opts$database

library(clusterProfiler)
library(dplyr)
library(tidyr)

data_gene = read.table(gene_file, header=F, sep = '\t', check.names=F, comment.char='', stringsAsFactors=F, quote = '')
kegg_data = read.table(database_file, header=T, sep = '\t', check.names=F, comment.char='', stringsAsFactors=F, quote = '')  

if(! 'Gene_ID' %in% colnames(kegg_data) || !'Pathway' %in% colnames(kegg_data) || !'koid' %in% colnames(kegg_data))
{
    message('[Error] database文件 必须含有 Gene_ID、koid、Pathway 列，且符合格式定义!')
    q()
}

# 数据库格式化
database = kegg_data %>% 
    select(Gene_ID, koid, Pathway) %>% # 选择Pathway列
    separate_rows(Pathway, sep = ';', convert = F) %>%  # 按分号拆掉
    filter(Pathway != '')  %>%  # 却掉空值
    separate(col = Pathway, into = c('ID', 'Description'), sep = '[|]', convert=FALSE)   # 把ID/Description 拆成两列

# koid 与 gene_id映射关系
map_geneid_koid = as.data.frame(database[, c('Gene_ID', 'koid')])
map_geneid_koid = map_geneid_koid[!duplicated(map_geneid_koid$Gene_ID), ]
rownames(map_geneid_koid) = map_geneid_koid$Gene_ID


# 清理输入的基因
gene = data_gene[, 1] 
gene = gene[!duplicated(gene)]  
gene = gene[gene %in% database$Gene_ID]


# 富集分析
enrich <- enricher(gene = gene,
                  TERM2GENE = database[c('ID','Gene_ID')],
                  TERM2NAME = database[c('ID','Description')],
                  pvalueCutoff = 1,
                  pAdjustMethod = 'BH',
                  qvalueCutoff = 1) 

result = as.data.frame(enrich)

# 添加KOID
result$KOID <- unlist(lapply(1:length(result$geneID), function(row){
                   paste(map_geneid_koid[unlist(strsplit(result$geneID[row], "/")), 'koid'], collapse = "/" )
                     }))

# 添加在线查看地址，标记颜色,基于geneid的
result$keggLink <- unlist(lapply(1:length(result$geneID), function(row){
                        koid = unlist(strsplit(result$KOID[row], "/"))
                        koid = koid[!duplicated(koid)]
                        gene_color = paste(paste(koid,    '%09red', sep = ''), collapse = '/')
                        # 最终地址
                        paste0('https://www.kegg.jp/kegg-bin/show_pathway?', result$ID[row], '/', gene_color)
                   }))

# （1）富集通路点图
kegg_pdf <- paste(output_dir, "kegg_dotplot.pdf", sep="/")
pdf(kegg_pdf,height=7,width=12)
p <-dotplot(enrich)
print(p)
dev.off()

#（2）表格输出
kegg_enrichment <- paste(output_dir, "kegg_enrichment.xls", sep="/")
write.table(result, kegg_enrichment, sep="\t", quote=F, row.names=F)

