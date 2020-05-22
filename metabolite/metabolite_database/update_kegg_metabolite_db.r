#!/home/genesky/software/r/3.6.1/bin/Rscript
library(docopt)

"Usage: update_kegg_metabolite_db.r  --compound_desc <file> --pathway_desc <file> --pathway_compound_map <file> --kgml_dir <dir>  --output_file <file>
Options:
   --compound_desc <file>         代谢物描述文件，两列数据，含有表头，第一列是代谢物id, 第二列是代谢物描述信息
   --pathway_desc <file>          通路描述文件，两列数据，含有表头。第一列是通路id, 第二列是通路描述信息
   --pathway_compound_map <file>  通路与代谢物id映射关系，或者说是通路里包含哪些代谢物。两列数据，含有表头。第一列是通路id,第二列是通路下代谢物id编号，如果一个通路有n个代谢物，则会有n行。
   --kgml_dir <dir>               通路的kgml文件路径，文件名必须是pathwayid.kgml。这些文件用于生成拓扑结果，计算# relative betweenness centrality # out degree centrality
   --output_file <file>           输出文件" -> doc


opts                  <- docopt(doc, version='甘斌， 生成R格式数据库\n')
compound_desc         <- opts$compound_desc
pathway_desc          <- opts$pathway_desc 
pathway_compound_map  <- opts$pathway_compound_map 
kgml_dir              <- opts$kgml_dir 
output_file           <- opts$output_file

# compound_desc         <- "./test/kegg.hsa.compound.desc.txt"
# pathway_desc          <- "./test/kegg.hsa.pathway.desc.txt"
# pathway_compound_map  <- "./test/kegg.hsa.pathway.compound.map.txt" 
# kgml_dir              <- "./test/tmp/hsa"

library(KEGGgraph)  # kgml文件解析工具
library(RBGL)  # 拓扑计算工具

# 读入数据
message('loading data')
data_pathway_desc         <- read.table(pathway_desc,         header = T, sep = "\t" ,stringsAsFactors = F, quote = "", comment.char ="")
pathway_compound_map      <- read.table(pathway_compound_map, header = T, sep = "\t" ,stringsAsFactors = F, quote = "", comment.char ="")
compound_desc             <- read.table(compound_desc,        header = T, sep = "\t" ,stringsAsFactors = F, quote = "", comment.char ="", row.names = 1)
 
message('building database: this may take some minutes')
# (1) 通路列表数据库
path.ids <- data_pathway_desc[, 1]
names(path.ids) <- data_pathway_desc[, 2]

# (2) 通路包含的代谢物数据库
mset.list <- list()  # 记录通路的代谢物id
rbc.list  <- list()  # 记录通路的代谢物拓扑值
for(row in 1:length(path.ids))
{  
	pathway_id <- path.ids[row]  # 通路id
   
   # 代谢物id获取
	compound <- pathway_compound_map[pathway_compound_map[, 1] == pathway_id, 2]  # 通路包含的代谢物id
	names(compound) <- compound_desc[compound, 1]
	mset.list[[row]] <- compound

   # 代谢物拓扑计算
   kgml <- paste0(kgml_dir, "/", pathway_id, ".kgml")
   map  <- parseKGML(kgml)
   cg   <- KEGGpathway2reactionGraph(map)  # 提取代谢物graph

   if(is.null(cg))  # 可能所有的代谢物彼此都没有联系
   {
      rbc_tmp = rep(0, length(compound))  # 都设为0
      names(rbc_tmp) <- compound
   }else {
      bcc  <- brandes.betweenness.centrality(cg)
      rbc_tmp <- bcc$relative.betweenness.centrality.vertices[1,]
      rbc_tmp <- rbc_tmp / sum(rbc_tmp)  # 使求和为 1
      names(rbc_tmp) <- gsub("cpd:","", x= names(rbc_tmp))  # 去掉多余的前缀
   }
   

   # 初始化
   rbc = rep(0, length(compound))
   names(rbc) <- compound
   rbc[names(rbc_tmp)] <- rbc_tmp  # 填充。上一步解析出的代谢图cg中，部分kgml中的代谢物没有包含在内，因为从通路图结构上来看，他们没有彼此相连。故，这里采用这种形式制作rbc
   rbc.list[[row]] <- rbc

}
names(mset.list) <- path.ids
names(rbc.list) <- path.ids

# (3) 数据整合
metpa <- list(mset.list = mset.list, path.ids = path.ids, rbc.list = rbc.list)

# (4) 保存
save(metpa, file=output_file)





