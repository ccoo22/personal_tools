#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: MetaboAnalyst_pca -i <file> --pcaplot <file> --biplot <file> --pca_coordinate <file> --sample_group_file <file> [--id_file <file> --color <string> --show_confidence]
Options:
	-i, --input <file>                 the input file, each row is a feature, each column is sample
	--pcaplot <file>                   the output pca pdf file
	--biplot <file>                    the output biplot pdf file
	--sample_group_file <file>         the input group file, column 1 is sampleinfo, column 2 is groupinfo
	--pca_coordinate <file>            pca coordinate file
	--color <string>                   format: \'group1,red;group2,blue;group3,green\' [default: NA]
    --id_file <file>                   需要绘制的特征值ID文件[default: all]
	--show_confidence                  if show confidence in picture" -> doc

opts              <- docopt(doc, version = 'Program : MetaboAnalyst PCA plot v1.0 \n      lhj 272\n')
inputfile         <- opts$input
pcaplot           <- opts$pcaplot
biplot            <- opts$biplot
pca_coordinate    <- opts$pca_coordinate
groupfile         <- opts$sample_group_file
color             <- opts$color
confidence        <- opts$show_confidence
id                <- opts$id_file


# 测试用参数
# inputfile         <- '/home/lhj/research/MetaboAnalyst/PCA/input/normalized_data.txt'
# pcaplot           <- '/home/lhj/research/MetaboAnalyst/PCA/output/PCA.pdf'
# biplot            <- "/home/lhj/research/MetaboAnalyst/PCA/output/PCAbiplot.pdf"
# pca_coordinate    <- '/home/lhj/research/MetaboAnalyst/PCA/output/pca.csv'
# groupfile         <- '/home/lhj/research/MetaboAnalyst/PCA/input/sample_group.xls'
# color             <- '1,Red;3,Black;2,Purple;4,Gray;5,skyblue;6,orangered;7,White;8,Yellow;9,Blue'
# confidence        <- TRUE

# options(warn=-1) #忽视任何警告
message("\nstart MetaboAnalyst_pca.r\n")
library(ggplot2, quietly = TRUE) 
data <- read.table(inputfile, sep='\t', header = T, row.names = 1, check.names=FALSE, stringsAsFactors = FALSE, quote = "") # 读取定量数据，列为样本
sample.groups <- read.table(groupfile, sep = "\t", header = F, quote="", stringsAsFactors = FALSE, colClasses = 'character') # 读取样本以及分组数据
feature_id       = rownames(data)
# PCA分析要求
# （1）输入数据，行为特征值，列为样本
# （2）每一列不能存在缺失数据

# 获取指定分析的样本
# 判断指定的样本是否存在于原始数据
sample_need  = sample.groups[,1]
if (length(sample_need[!(sample_need %in% colnames(data))]) != 0)
{
	cat('[Error] 分组文件中，部分样本编号在丰度文件中不存在! 请仔细检查分组文件 !\n')
	q()
}
if(id != 'all')
{
    data_id       <- read.table(id, header = F, quote = "")
    feature_id    <- as.character(unlist(data_id))
    if(sum(feature_id %in% rownames(data)) != length(feature_id))
    {
        lost_ID = feature_id[!feature_id %in% rownames(data)]
        message("[Error] 部分输入的特征ID在输入文件中没有找到 : ", lost_ID)
        q()
    }
}
data_need    = data[feature_id, sample_need, drop = FALSE]
sample_group = sample.groups[,2]

# 判断分组信息是否有缺失值
if (NA %in% sample_group)
{
	cat('[Error] please check the group_name_list !\n')
	q()
}

# pca
pca <- prcomp(t(data_need))
summaryInfo <- summary(pca)
pc1Prop <- 100 * summaryInfo$importance['Proportion of Variance', 'PC1'] # PC1主元解释百分比
pc2Prop <- 100 * summaryInfo$importance['Proportion of Variance', 'PC2'] # PC2主元解释百分比
coordinate <- as.data.frame(pca$x)
coordinate$Group = factor(sample_group, levels = unique(sample_group)) # factor数据，levels与输入顺序一致，以固定pca图上的分组顺序

# 输出坐标
write.csv(coordinate, file=pca_coordinate)

###plot###
message('开始绘图')

if (color != 'NA')
{
	message("使用指定颜色")
	Groups = unique(sample.groups[,2])
	group_col = unlist(strsplit(unlist(strsplit(color, ';')), ','))
	if (length(group_col) != 2 * length(Groups))
	{
		message("指定颜色数量与分组数不一致")
		q()
	}

	# 根据分组映射颜色，指定names实现
	color_define = t(matrix(group_col,nrow=2))
	mycol        = color_define[, 2]
	names(mycol) = color_define[, 1]

}else{
	message("使用默认颜色设置")
	mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)  
	mycol <- colors()[rep(mycol, 50)]
}

# 形状模版
myshape <- rep(c(15,16,17,18,19,20,21,22,23,24,25,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14),4)

# 绘制pca
pdf(file=pcaplot, width=13, height=13)
# 判断是否展示置信区间
if(confidence == TRUE)
{
	message("展示95%置信区间")
	p = ggplot(data=coordinate, aes(x = PC1, y = PC2, colour = Group, shape = Group)) +       # 映射
        geom_hline(yintercept = 0, colour = "gray65") +                                   # 添加水平线
        geom_vline(xintercept = 0, colour = "gray65") +                                   # 添加垂直线
        scale_shape_manual(values = myshape ) +                                           # 指定形状
        scale_color_manual(values = mycol ) +                                             # 指定颜色
        geom_point(size = 3, alpha = 1) +                                                 # 设置点的大小、透明度
        stat_ellipse(level = 0.95, show.legend = F)  +                                    # 添加置信区间
        ggtitle("PCA plot of MetaboAnalyst") + xlab(paste("PC1 ", "(", pc1Prop, "%)", sep="")) + ylab(paste("PC2 ", "(", pc2Prop, "%)", sep="")) + theme(plot.title = element_text(hjust = 0.5))   # 定义标题、横纵坐标
}else{
	message("不展示置信区间")
	p = ggplot(data=coordinate, aes(x = PC1, y = PC2, colour = Group, shape = Group)) +
        geom_hline(yintercept = 0, colour = "gray65") +
        geom_vline(xintercept = 0, colour = "gray65") +
        scale_shape_manual(values = myshape ) +
        scale_color_manual(values = mycol ) +
        geom_point(size = 3, alpha = 1) +  
        ggtitle("PCA plot of MetaboAnalyst") + xlab(paste("PC1 ", "(", pc1Prop, "%)", sep="")) + ylab(paste("PC2 ", "(", pc2Prop, "%)", sep="")) + theme(plot.title = element_text(hjust = 0.5))
}

p
dev.off()

# 绘制biplot
choices = c(1, 2)                                 # PC1、PC2
score   <- pca$x
lam     <- pca$sdev[choices]                      # 提取 $sdev 中PC1、PC2的值
n       <- nrow(score)
lam     <- lam * sqrt(n)
pdf(file=biplot, width=13, height=13)
x       <- t(t(score[, choices]) / lam)
y       <- t(t(pca$rotation[, choices]) * lam)    # 提取 $rotation 中PC1、PC2的值并进行处理
rownames(y) = data[feature_id, 1]                          # 给定行名，用于图中展示
biplot(x, 
	   y, 
	   xpd =T, 
	   cex=0.9, 
	   xlab=paste("PC1 ", "(", pc1Prop, "%)", sep=""), 
	   ylab=paste("PC2 ", "(", pc2Prop, "%)", sep=""))
dev.off()
message("\nfinish MetaboAnalyst_pca.r\n")
