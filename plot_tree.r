#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: pca_tree.r -t <file>  -o <file> [-g <file> --text_size <num>  ]
Options:
   -t, --tree_file <file>          tree 文件，newick格式
   -g, --group_file <file>         样本分组文件，两列，第一列样本名，第二列样本分组，必须含有表头，表头名字随意。 [default: NA]
   -o, --output <pdf>              pdf output file
   --text_size <num>               sample name size, default: auto calculate" -> doc

opts                     <- docopt(doc, version = 'Program : plot tree v1.0 \n          甘斌 129\n')
tree_file                <- opts$tree_file
group_file               <- opts$group_file
output                   <- opts$output
text_size                <- opts$text_size


# # 测试用参数
# tree_file                 <- '/home/ganb/work/tmp/17B1116A_18B0511B/17B1116A/output/sequence.fa.tre'
# group_file               <- '/home/ganb/work/tmp/17B1116A_18B0511B/17B1116A/group.txt'
# output                   <- './test.pdf'
 
options(warn=-1)
library(ggplot2, quietly = TRUE)
library(ggtree, quietly = TRUE)
library(ape, quietly = TRUE)


# 读入树
message("load tree")
tree = read.tree(tree_file)
samples = tree$tip.label

# 读入分组数据
if(group_file != 'NA')
{
  message("load group")
  group_data = read.table(group_file, head = T, check.names = F, stringsAsFactors=F)
  colnames(group_data)[1:2] = c('Sample', 'Group')
  rownames(group_data) = group_data$Sample

  # 取出分组名，没有的话，定义为NO_GROUP
  group_names = group_data[samples, 'Group']
  group_names[is.na(group_names)] = 'NO_GROUP'

  # 进化树上添加分组
  group_list = lapply(unique(group_names), function(group_name){  samples[group_names == group_name]   })
  names(group_list) = unique(group_names)
  tree <- groupOTU(tree, group_list, 'Group') 
}


# 颜色模版
mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)  
mycol <- colors()[rep(mycol, 50)] 

# 自适应大小
resizeTXT = round(300 / length(samples), 4) # 字体大小根据样本数而变
if(resizeTXT > 2) resizeTXT = 2 # 字体大小根据样本数而变
if(! is.null(text_size) ) resizeTXT = as.numeric(text_size) # 自定义字体大小

resizePDF = round(length(samples) / 20) # PDF尺寸大小根据样本数而变
if(resizePDF < 8) resizePDF = 8  # 最小尺寸是8


# 绘图
pdf(output, width = resizePDF, height = resizePDF)
par(mar = c(5,5, 5, 15))
message("plot tree")
pic = list()
if(group_file != 'NA')
{
  pic[[1]] = ggtree(tree, aes(color=Group), layout="circular", branch.length="none" ) + geom_tiplab2(size=resizeTXT, aes(angle = angle) ) + geom_text2(aes(label=as.numeric(label)*100,subset=!isTip), hjust=-.3,size = resizeTXT * 0.6) + theme(legend.position = "right") + ggtitle("(Phylogram) circular layout and no branch.length")    +  scale_color_manual(values = mycol[1:length(group_list)]) 
  pic[[2]] = ggtree(tree, aes(color=Group), layout="circular" )                       + geom_tiplab2(size=resizeTXT, aes(angle = angle))  + geom_text2(aes(label=as.numeric(label)*100,subset=!isTip), hjust=-.3,size = resizeTXT * 0.6) + theme(legend.position = "right") + ggtitle("(Phylogram) circular layout")                         +  scale_color_manual(values = mycol[1:length(group_list)])  
  pic[[3]] = ggtree(tree, aes(color=Group))                                        + geom_tiplab(size=resizeTXT)                       + geom_text2(aes(label=as.numeric(label)*100,subset=!isTip), hjust=-.3,size = resizeTXT * 0.6) + theme(legend.position = "right") + ggtitle("(Phylogram) rectangular layout")                      +  scale_color_manual(values = mycol[1:length(group_list)]) 
  pic[[4]] = ggtree(tree, aes(color=Group), branch.length="none")                     + geom_tiplab(size=resizeTXT)                       + geom_text2(aes(label=as.numeric(label)*100,subset=!isTip), hjust=-.3,size = resizeTXT * 0.6) + theme(legend.position = "right") + ggtitle("(Phylogram) rectangular layout and no branch.length") +  scale_color_manual(values = mycol[1:length(group_list)]) 
}else {
  pic[[1]] = ggtree(tree, layout="circular", branch.length="none" ) + geom_tiplab2(size=resizeTXT, aes(angle = angle) ) + geom_text2(aes(label=as.numeric(label)*100,subset=!isTip), hjust=-.3,size = resizeTXT * 0.6) + theme(legend.position = "right") + ggtitle("(Phylogram) circular layout and no branch.length")      
  pic[[2]] = ggtree(tree, layout="circular" )                        + geom_tiplab2(size=resizeTXT, aes(angle = angle))  + geom_text2(aes(label=as.numeric(label)*100,subset=!isTip), hjust=-.3,size = resizeTXT * 0.6) + theme(legend.position = "right") + ggtitle("(Phylogram) circular layout")                           
  pic[[3]] = ggtree(tree)                                            + geom_tiplab(size=resizeTXT)                       + geom_text2(aes(label=as.numeric(label)*100,subset=!isTip), hjust=-.3,size = resizeTXT * 0.6) + theme(legend.position = "right") + ggtitle("(Phylogram) rectangular layout")                        
  pic[[4]] = ggtree(tree, branch.length="none")                      + geom_tiplab(size=resizeTXT)                       + geom_text2(aes(label=as.numeric(label)*100,subset=!isTip), hjust=-.3,size = resizeTXT * 0.6) + theme(legend.position = "right") + ggtitle("(Phylogram) rectangular layout and no branch.length")   

}

print(pic[[1]]) 
print(pic[[2]]) 
print(pic[[3]]) 

print(pic[[4]]) 
# 创建图例、分组设置不同颜色、设置图形类型，设置分支长度                         文字大小 geio_tiplab2会产生倾斜角度                   显示bootstrap信息                                                                               分组提示位置                         标题                                                             设置分组颜色
dev.off()
