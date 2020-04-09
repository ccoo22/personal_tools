#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)

"Usage: clusterProfiler_selfdb.R  -g <file> --pathway_desc <file> --pathway_geneid_map <file> --geneid_symbol_map <file> --output_dir <dir>
Options:
   -g, --gene <file>              需要分析的基因列表文件（genesymbol）,没有表头，第一列必须是基因名
   --pathway_desc <file>          通路描述文件，两列数据，含有表头。第一列是通路id, 第二列是通路描述信息
   --pathway_geneid_map <file>    通路与基因id映射关系，或者说是通路里包含哪些基因。两列数据，含有表头。第一列是通路id,第二列是通路下基因id编号，如果一个通路有n个基因，则会有n行。
   --geneid_symbol_map <file>     基因id与genesymbol映射关系。三列数据，含有表头。第一列是基因id, 第二列genesymbol，第三列是该symbol的优先级，值越小越优先，用于应对同一个symbol有不同的基因id的情况，如果symbol相同，优先级也相同，那就只有听天由命了。
                                  因为通路中的基因id要求唯一性，而我们经常使用的genesymbol往往含有别名，故通过该映射文件把输入的genesymbol转换为唯一基因id。
                                  但是，在kegg的注释里，一个基因名可能会有多个不同的id（历史遗留问题，或者别名导致的）, 故需要有一个权重数据协助选择，这就是第三列的作用了， 默认选择编号较小的。
   --output_dir <dir>             结果输出目录" -> doc


opts                <- docopt(doc, version='甘斌，基于自建数据库的通路分析\n')
gene                <- opts$gene
pathway_desc        <- opts$pathway_desc 
pathway_geneid_map  <- opts$pathway_geneid_map 
geneid_symbol_map   <- opts$geneid_symbol_map 
output_dir          <- opts$output_dir

dir.create(output_dir, showWarnings = FALSE)
library(clusterProfiler)

# pathway_desc       = "/home/ganb/work/research/kegg/self_database/kegg_hsa_20200408/hsa.pathway.desc.txt"
# pathway_geneid_map = "/home/ganb/work/research/kegg/self_database/kegg_hsa_20200408/hsa.pathway.geneid.map.txt"
# geneid_symbol_map  = "/home/ganb/work/research/kegg/self_database/kegg_hsa_20200408/hsa.geneid.symbol.map.txt"
# gene = "/home/ganb/work/research/kegg/self_database/example_gene.txt"

# 读入数据
message('loading data')
data_pathway_desc       <- read.table(pathway_desc,        header = T, sep = "\t" ,stringsAsFactors = F, quote = "", comment.char ="")
data_pathway_geneid_map <- read.table(pathway_geneid_map,  header = T, sep = "\t" ,stringsAsFactors = F, quote = "", comment.char ="")
data_geneid_symbol_map  <- read.table(geneid_symbol_map,   header = T, sep = "\t" ,stringsAsFactors = F, quote = "", comment.char ="")
data_gene               <- read.table(gene,                header = F, sep = "\t" ,stringsAsFactors = F, quote = "", comment.char ="")
colnames(data_geneid_symbol_map)     <- c("geneid", 'genesymbol', 'count')

# 取unique,且大写,防止输入重复
data_gene_symbol <- toupper(unique(data_gene[, 1])) 

message('match gene symbol')
# 提取symbol 与 geneid 映射关系，并只保留映射成功的
match_1 <- data_geneid_symbol_map[data_geneid_symbol_map[, 'genesymbol'] %in% data_gene_symbol, ]  # 找出所有symbol匹配的
match_2 <- match_1[order(match_1$count), ]  # 按照count排序
match_final <- match_2[!duplicated(match_2$genesymbol), ]  # 重复的只保留count最小的


if(nrow(match_final) == 0)
{
    message("没有匹配到基因id, 无法做富集分析. 请仔细检查输入的基因名是否与data_geneid_symbol_map文件中的genesymbol列对应")
    q()
}

# 富集分析
message('enrichment analysis')
kk = enricher(match_final[, 'geneid'], 
	pvalueCutoff = 1,
	qvalueCutoff = 1,
	minGSSize = 1,
	TERM2GENE = data_pathway_geneid_map, 
	TERM2NAME = data_pathway_desc)

message('output')
kk_enrich = as.data.frame(kk)

# geneid 重新转换为genesymbol
kk_enrich$geneID <- unlist(
                           lapply(1:nrow(kk_enrich), function(row){  
                                                                paste(match_final[match_final[, 'geneid'] %in% unlist(strsplit(kk_enrich$geneID[row], "/")), 'genesymbol'], collapse = "/")  
                                                                }
                                 )
                           )
# 输出
kegg_pdf <- paste(output_dir, "kegg_dotplot.pdf", sep="/")
kegg_enrichment <- paste(output_dir, "kegg_enrichment.xls", sep="/")

if(nrow(kk_enrich) >= 1)
{
	pdf(kegg_pdf, height=7, width=12)
	p <-dotplot(kk)
	print(p)
	dev.off()
	write.table(kk_enrich, kegg_enrichment, sep="\t", quote=F, row.names=F)
}




