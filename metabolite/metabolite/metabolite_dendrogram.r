#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: MetaboAnalyst_clustering -i <file> --clusterplot <file> --sample_group_file <file> --smplDist <string> --clstDist <string> [--color <string>]
Options:
	-i, --input <file>                 the input file, each row is a feature, each column is sample
	--clusterplot <file>               the output pca pdf file
	--sample_group_file <file>         the input group file, column 1 is sampleinfo, column 2 is groupinfo
	--smplDist <string>                the method to calculate sample distance: euclidean, binary, pearson, kendall, spearman
	--clstDist <string>                the method to calculate clustering distance: ward.D, ward.D2, single, complete, average(UPGMA), mcquitty(WPGMA), median(WPGMC), centroid(UPGMC)
	--color <string>                   format: \'group1,red;group2,blue;group3,green\' [default: NA]" -> doc

opts              <- docopt(doc, version = 'Program : MetaboAnalyst Dendrogram plot v1.0 \n      lhj 272\n')
inputfile         <- opts$input
clusterplot       <- opts$clusterplot
groupfile         <- opts$sample_group_file
color             <- opts$color
smplDist          <- opts$smplDist
clstDist          <- opts$clstDist

# 测试用参数
# inputfile         <- "/home/lhj/research/MetaboAnalyst/Dendrogram/input/normalized_data.txt"
# clusterplot       <- "/home/lhj/research/MetaboAnalyst/Dendrogram/output/Sampletree.pdf"
# groupfile         <- "/home/lhj/research/MetaboAnalyst/Dendrogram/input/sample_group_4.xls"
# color             <- "B1,Red;9,Blue;6,Black;A4,Gray"
# smplDist          <- "pearson"
# clstDist          <- "average"
# # color             <- "NA"

message("\nstart MetaboAnalyst_clustering.r\n")
data <- read.table(inputfile, sep = '\t', header = TRUE, check.names = FALSE) # 读取定量数据，列为样本
sample.groups <- read.table(groupfile, sep = "\t", header = FALSE, quote = "", stringsAsFactors = FALSE) # 读取样本以及分组数据

# 获取指定分析的样本
# 判断指定的样本是否存在于原始数据
sample_need  = sample.groups[,1]
if (length(sample_need[!(sample_need %in% colnames(data))]) != 0)
{
	cat('[Error] please check the sample_name_list !\n')
	q()
}
data_need    = data[, sample_need, drop = FALSE]
sample_group = as.character(sample.groups[,2])

# 判断分组信息是否有缺失值
if (NA %in% sample_group)
{
	cat('[Error] please check the group_name_list !\n')
	q()
}

data_new <- t(as.matrix(data_need))               # 预处理完用于后续分析的数据
colnames(data_new) <- data[,1]                    # 原始数据id号（第一列）作为列名,确保无重复

# set up distance matrix
if(smplDist == 'euclidean')
{
	dist.mat <- dist(data_new, method = smplDist)
}else{
	dist.mat <- dist(1-cor(t(data_new), method = smplDist))
}

hc_tree <- hclust(dist.mat, method = clstDist)

pdf(clusterplot, width = 13, height = 13)
par(cex = 0.8, mar = c(4,2,2,8))
clusDendro <- as.dendrogram(hc_tree)

if (color == 'NA')
{
	message("使用默认颜色设置")
	mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)  
	mycol <- colors()[rep(mycol, 50)]
 
	# 对应分组颜色信息给每个样本填颜色，得到相应颜色向量
	ordergroup = unique(sample_group)
	ordercolor = rep(1, length(sample_group))
	for(i in 1:length(ordergroup))
	{
		ordercolor[sample_group %in% ordergroup[i]] = mycol[i]
	}
}else{
	message("使用指定颜色")
	ordergroup = unique(sample_group)
	group_col = unlist(strsplit(unlist(strsplit(color, ';')), ','))
	if (length(group_col) != 2 * length(ordergroup))
	{
		message("指定颜色数量与分组数不一致")
		q()
	}
	color_define = t(matrix(group_col, nrow = 2))
	mycol        = color_define[, 2]
	names(mycol) = color_define[, 1]  # 组名与颜色一一对应

	# 对应分组颜色信息给每个样本填颜色，得到相应颜色向量
	ordercolor = rep(1, length(sample_group))
	for(i in 1:length(ordergroup))
	{
		ind <- which(names(mycol) == ordergroup[i])
		ordercolor[sample_group %in% ordergroup[i]] = mycol[ind]
	}
}

# 定标签颜色
names(ordercolor) <- sample_need
labelColors <- ordercolor[hc_tree$order]

colLab <- function(n)
{
    if(is.leaf(n))
    {
    	a <- attributes(n)
        labCol <- labelColors[a$label]
        attr(n, "nodePar") <- 
          if(is.list(a$nodePar)) c(a$nodePar, lab.col = labCol, pch = NA) else
            list(lab.col = labCol, pch = NA)
    }
    n
}
clusDendro <- dendrapply(clusDendro, colLab) # 对谱系图所有节点应用函数
plot(clusDendro, horiz = T, axes = T)
par(cex = 1)

legend("topleft", legend = unique(sample_group), pch = 15, col = unique(ordercolor), bty = "n")  # 标注分组信息
dev.off()

message("\nfinish MetaboAnalyst_clustering.r\n")