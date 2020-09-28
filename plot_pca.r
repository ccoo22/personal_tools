#!/home/genesky/software/r/3.5.1/bin/Rscript

.libPaths("/home/genesky/software/r/3.5.1/lib64/R/library/")
library(docopt)
"Usage: pca_plot.r --pca_file <file> --group_file <file> --output <pdf>  [--size <int> ] 
Options:
   --pca_file <file>           PCA坐标文件，第一列样本名，第二、三列分别是pc1/pc2坐标。含有表头。
   --group_file <file>         分组文件，含有表头。最多三列，可以只放入两列。第一列样本名，第二列分组1，图中用颜色区分。第三列分组2，图中用形状区分。
   --size <int>                点的大小 [default: 1]
   --output <pdf>              pdf output file" -> doc

opts                     <- docopt(doc, version = 'Program : plot pca v1.0 \n          甘斌 129\n')
pca_file                 <- opts$pca_file
group_file               <- opts$group_file
output                   <- opts$output
size                     <- as.numeric(opts$size)
 
# # 测试用参数
# pca_file                 <- '/home/ganb/work/tmp/17B1116A_18B0511B/18B0511B/output/sample_pca.txt'
# group_file               <- '/home/ganb/work/tmp/17B1116A_18B0511B/18B0511B/group.txt'
# output                   <- './test.pdf'

library(ggplot2) # 正态分布检测要用

pca_data   = read.table(pca_file, head = T, check.names = F, stringsAsFactors=F, sep='\t')
group_data = read.table(group_file, head = T, check.names = F, stringsAsFactors=F, sep='\t', colClasses = 'character')
group_col_names = colnames(group_data)
colnames(pca_data)[1:3] = c('Sample', 'PC1', 'PC2')
colnames(group_data)[1:2] = c('Sample', 'Group')
pca_data$Sample = as.character(pca_data$Sample)
group_data$Sample = as.character(group_data$Sample)

if(ncol(group_data) == 3) colnames(group_data)[3] = 'Group2'
rownames(pca_data)   = pca_data$Sample
rownames(group_data) = group_data$Sample


# 取出分组名，没有的话，定义为NO_GROUP

pca_data$Group = group_data[pca_data$Sample, 'Group']
pca_data$Group[is.na(pca_data$Group)] = 'NO_GROUP'
pca_data$Group = factor(pca_data$Group, levels=unique(c(group_data$Group, 'NO_GROUP')))
if(ncol(group_data) == 3)
{
    pca_data$Group2 = group_data[pca_data$Sample, 'Group2']
    pca_data$Group2[is.na(pca_data$Group2)] = 'NO_GROUP'
    pca_data$Group2 = factor(pca_data$Group2, levels=unique(c(group_data$Group2, 'NO_GROUP')))
}



# 颜色模版
mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)  
mycol <- colors()[rep(mycol, 50)] 

# 形状模版
myshape <- c(15, 16, 17, 18, 0:14)

pdf(output)

if(ncol(group_data) == 2) p <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = Group)) + geom_point(size = size) + scale_color_manual(values = mycol[1:length(levels(pca_data$Group))]) + guides(color = guide_legend(title = group_col_names[2]))
# if(ncol(group_data) == 3) p <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = Group, shape = Group2)) + geom_point(size = size)  + scale_color_manual(values = mycol[1:length(levels(pca_data$Group))]) + scale_shape_manual(values=myshape[1:length(levels(pca_data$Group2))])  + guides(color = guide_legend(title = group_col_names[2]), shape = guide_legend(title = group_col_names[3]))

p
dev.off()
