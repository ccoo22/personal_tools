#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: string.r  -g <file> -o <dir> --case_group_name <string> --control_group_name <string> [--cor_file <file> --pvalue_cutoff <numeric> --log2fc_cutoff <numeric> --Rlib <dir>]
Options:
    -g <file>                       基因列表，一列数据，没有表头
    -o, --output_dir <dir>          结果输出目录
    --species_id <string>           物种ID编号，人的是 9606， 小鼠的是 10090， 大鼠的是 10116。可以在https://stringdb-static.org/download/species.v11.0.txt 查找。
    --ppi_db_dir <dir>              string本地数据库目录，包含三个string分析的关键文件 protein_aliases_tf.tsv.gz protein_links.tsv.gz proteins.tsv.gz
                                    例如： 
                                    人 /home/genesky/database_new/string_ppi/Homo_sapiens_9606
                                    小鼠 /home/genesky/database_new/string_ppi/Mus_musculus_10090
                                    大鼠 /home/genesky/database_new/string_ppi/Rattus_norvegicus_10116
    --Rlib <dir>                    R包路径  [default: /home/genesky/software/r/3.5.1/lib64/R/library] 
    
    注：不建议使用，因为R包只支持10.0版本，https://string-db.org 已经更新到11.0
    " -> doc

opts                     <- docopt(doc, version = 'string ppi 分析 \n          甘斌 129\n')
gene                     <- opts$gene
output_dir               <- opts$output_dir
species_id               <- as.numeric(opts$species_id)
ppi_db_dir               <- opts$ppi_db_dir
Rlib                     <- opts$Rlib
.libPaths(Rlib)


library(STRINGdb)
library(igraph)


input = read.table(gene, header = F, sep = '\t', check.names = F, comment.char = '', quote = '', fill = T)
colnames(input)[1] = 'gene'
 
#加载并初始化数据库
message('load string db')
string_db   = STRINGdb$new(version = '10', species = species_id, score_threshold = 0, input_directory = ppi_db_dir)
# 将gene名和数据库中的蛋白标识符mapping
message('mapping string id')
data_mapped = string_db$map(input,'gene',removeUnmappedRows = TRUE)
rownames(data_mapped) = data_mapped$STRING_id
# 提取映射后的stringid
hits        = data_mapped$STRING_id



# PPI 绘图
message('plot string network ')
if(length(hits) < 400) {
    pdf(paste(output_dir, 'ppi_string.pdf',sep = '/'))
    string_db$plot_network(hits) # ppi分析，绘图，该步骤会从官网下载数据库，大约35MB
    dev.off()
}else{
    message('[warning] 抱歉，无法绘图，因为基因数量超过400个')
}

# PPI 基因间打分
inter    = string_db$get_interactions(hits)
node1    = data_mapped[inter$from, 'gene']
node2    = data_mapped[inter$to, 'gene']

inter_res = cbind(node1, node2, inter)
write.table(inter_res, paste(output_dir, 'ppi_interactions.xls', sep = '/'), sep = '\t', quote = F, row.names = F)


#igraph构建网络
net = graph.data.frame(inter_res, directed = F)

# 网络拓扑属性 
network_attr <- data.frame(vcount(net), ecount(net), edge_density(net), transitivity(net, type='global'))
colnames(network_attr) <- c('vertices', 'edges', 'density', 'clustering coefficient')

out_file <- paste(output_dir, 'ppi_network.attr.xls', sep = '/')
write.table(network_attr, out_file, sep = '\\t', row.names = F, quote = F)

# 网络拓扑 节点属性
node_attr <- data.frame(
    degree(net, mode='all'),
    betweenness(net),
    closeness(net),
    transitivity(net, type='local')
)
colnames(node_attr) <- c('degree', 'betweenness centrality', 'closeness centrality', 'clustering coefficient')
out_file <- paste(output_dir, 'ppi_node.attr.xls', sep = '/')
write.table(node_attr, out_file, sep = '\\t', quote = F, col.names = NA)


