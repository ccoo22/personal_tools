#!/usr/bin/env Rscript
.libPaths("/home/genesky/software/r/3.5.1/lib64/R/library/")
library(docopt)
"Usage: pca_plot.r [options] 
Options:
   --pca_file <file>           pca file with 3 column:sample pc1 pc2. must with head
   --group_file <file>         group file with 2 column:sample group. must with head
   --output <pdf>              pdf output file" -> doc

opts                     <- docopt(doc, version = 'Program : plot pca v1.0 \n          甘斌 129\n')
pca_file                 <- opts$pca_file
group_file               <- opts$group_file
output                   <- opts$output


# # 测试用参数
# pca_file                 <- '/home/ganb/work/tmp/17B1116A_18B0511B/18B0511B/output/sample_pca.txt'
# group_file               <- '/home/ganb/work/tmp/17B1116A_18B0511B/18B0511B/group.txt'
# output                   <- './test.pdf'

library(ggplot2) # 正态分布检测要用

pca_data   = read.table(pca_file, head = T, check.names = F, stringsAsFactors=F)
group_data = read.table(group_file, head = T, check.names = F, stringsAsFactors=F)
colnames(pca_data)[1:3] = c('Sample', 'PC1', 'PC2')
colnames(group_data)[1:2] = c('Sample', 'Group')


# 取出分组名，没有的话，定义为NO_GROUP
group_name = unlist(lapply(pca_data$Sample, function(sample_name){ 
        group = 'NO_GROUP'
        if(sample_name %in% group_data$Sample) group = group_data[which(group_data$Sample == sample_name), 'Group']
        group  
    } ) )
pca_data$Group = factor(group_name)

# 颜色模版
mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)  
mycol <- colors()[rep(mycol, 50)] 

pdf(output)
p <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = Group)) + geom_point() + scale_color_manual(values = mycol[1:length(levels(pca_data$Group))])
p
dev.off()