library(docopt)

# TEST RUN:
# /home/genesky/software/r/3.5.1/bin/Rscript TreeBar.r -t /home/yucy/TEST_DATASET/species.taxon.Abundance.txt -m AQ -g /home/yucy/TEST_DATASET/sample_groups.txt --output /home/yucy/TEST_DATASET/TEST_OUTPUT/

"Usage: TreeBar.r  -t <file> -g <file> -m <string> --output <dir> [-c <string> --lab <string> --lib <string>]
Options:
   -t , --taxonfile <file>        丰度文件， reads数量矩阵，第一列分类名称，第二列往后为每个样本的reads，第一行为样本名
   -m , --method <string>         数据类型：  AQ 或 RQ
   --lab <string>                 指定legend title [default: NA]
   -g , --groupfile <file>        样本分组文件，两列，第一列样本名，第二列样本分组，没有表头。
   -c , --color <string>          颜色列表，指定使用的色系，颜色之间用逗号分隔。 例如： #1F77B4FF,#FF7F0EFF,#2CA02CFF [default: NA]
   --lib <string>                 R包安装路径 lib [default: /home/genesky/software/r/4.0.3/lib64/R/library/]   
   --output <dir>                 数据结果目录" -> doc

opts                     <- docopt(doc, version = '群落分析TreeBar\n')
taxonfile                <- opts$taxonfile
method                   <- opts$method
lab                      <- opts$lab
groupfile                <- opts$groupfile
color                    <- opts$color
lib                      <- opts$lib
output                   <- opts$output

OUT_COMMAND <- function(){
    Rscript_str = paste(unlist(strsplit(commandArgs(trailingOnly = FALSE)[1], "/lib64/R/bin/exec/R"))[1], "/bin/Rscript", sep = "")
    library(funr); script_dir_str = sys.script()
    args = commandArgs(trailingOnly = T); parameter = ""; 
    for (arg in args) { 
      if(length(grep("#", arg)) > 0){
        parameter = paste(parameter, paste("'", arg, "'", sep = ""), sep = " ")
      }else{
        parameter = paste(parameter, arg, sep = " ")
      }
    }
    Command = paste("[Command]:", Sys.time(), Rscript_str, script_dir_str, parameter, sep = " ")
    if(exists("output")){
      log_file = paste(output, "command.info", sep = "/")
      if(file.exists(log_file)){
        write(Command, file = log_file, append = TRUE)
      }else{
        write(Command, file = log_file)
      }
    }else{
      write(Command, file = "command.info")
    }
}
OUT_COMMAND()


taxonfile = normalizePath(taxonfile)
groupfile = normalizePath(groupfile)
output    = normalizePath(output)
setwd(output)


.libPaths(lib)
library(vegan)
library(dplyr)
library(tidyr)
library(ggtree)
library(ggpubr)
library(cowplot)


# taxonfile                <- 'class.taxon.Abundance.txt'
# method                   <- 'RQ'
# lab                      <- 'NA'
# groupfile                <- 'sample_group.txt'
# color                    <- 'NA'
 
# ====================== prepare data =======================
if (method == "RQ"){
    if (lab == "NA") lab = "Relative abundance (%)"
}else if (method == "AQ"){
    if (lab == "NA") lab = "Absolute abundance"
}else{
    message("method 只能是AQ或RQ")
    q()
}

# 分组读入
SampleInfo = read.table(groupfile, header = F, sep = "\t", quote="", colClasses = "character")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo$Sample

# 丰度读入
taxondata = read.table(taxonfile, header = T, fill = T, check.names = F, sep = "\t", row.names = 1, comment.char = "", quote = "")

# 缺失样本检测
lost_samples = SampleInfo$Sample[!SampleInfo$Sample %in% colnames(taxondata)]
if (length(lost_samples) > 0){
    message("[Error] 部分分组文件中的样本在OTU文件中没有找到 : ", paste(lost_samples, collapse = ','))
    q()
}
taxondata = taxondata[, SampleInfo$Sample, drop = F]

# 处理掉行名特殊字符
rownames(taxondata) = gsub("\"","",rownames(taxondata))
colnames(taxondata) = gsub("\"","",colnames(taxondata))
rownames(taxondata) = as.character(strsplit(rownames(taxondata),split="{.*}",perl=T))
# 处理掉特殊的行
loc = grep("All|p__unclassified|c__unclassified|o__unclassified|f__unclassified|g__unclassified|s__unclassified", rownames(taxondata))
if(length(loc) > 0) data_clean = taxondata[-loc, ,drop = F]


# 最终 clean data 
data_clean = data_clean[rowSums(data_clean) != 0, ,drop = F]
if(nrow(data_clean) == 0){
    q()
}

# 绘图准备
# 相对量计算
data_relative_abundance = apply(data_clean, 2, function(x){ x / sum(x) })
if(nrow(data_clean) == 1) {  # 如果只有一行数据
    data_relative_abundance = data_clean
    data_relative_abundance[1, ] = 1
}

# RQ/AQ绘图采用的数据
data_plot_tmp = data_clean  # 绝对量
if(method == 'RQ') data_plot_tmp = data_relative_abundance  # 相对量

# 物种数量是否太多，做一些缩减处理
if(nrow(data_plot_tmp) > 30){
    index = which(apply(data_relative_abundance > 0.01, 1, any)==T)  # 取相对丰度 > 0.01的物种绘图
    data_high_abundance = data_plot_tmp[index,,drop = F]
    data_high_abundance = data_high_abundance[order(rowSums(data_high_abundance),decreasing = T),,drop =F]

    data_for_plot = rbind(data_high_abundance, apply(data_plot_tmp[-index,,drop = F], 2, sum))  # 低于 0.01的单独汇总为Others
    rownames(data_for_plot)[nrow(data_for_plot)] = "Others"
}else{
    data_for_plot = data_plot_tmp[order(rowSums(data_plot_tmp),decreasing = T),,drop =F]
} 
# No_Rank如果存在，把它放在末尾      
index_no_rank = grep("No_Rank", rownames(data_for_plot))
if(length(index_no_rank) > 0) {
    data_for_plot = rbind(as.matrix(data_for_plot[-index_no_rank,]), data_for_plot[index_no_rank,])
    rownames(data_for_plot)[nrow(data_for_plot)] = "No_Rank"
}
data_for_plot = as.data.frame(data_for_plot, check.names=F)  # 转成data.frame

# 开始准备绘图
####
# 进化树
###
fga.dist = vegdist(t(data_plot_tmp), method = "bray")
hc       = hclust(fga.dist, method = "average")
ddgram   = as.dendrogram(hc)
 
ggtree_plot <- ggtree::ggtree(ddgram) + ggtree::theme_tree2()  # + ggtree::geom_tiplab() 加样本标签，测试用
 
# 根据进化树中的样本顺序，调整绘图矩阵中的样本顺序
samples_with_order <- ggtree_plot$data %>% filter(!is.na(label)) %>% arrange(y) %>% .$label %>% as.character() # 提取ggtree绘图的样本名（已排序）
data_for_plot = data_for_plot[, samples_with_order , drop = F]
 
# 柱状图准备
# 颜色模版
if (color != 'NA'){
    mycol = c(unlist(strsplit(color, ',')))
    if (length(mycol) < nrow(data_for_plot)){
        message("指定颜色数量不足，需要数量：", nrow(data_for_plot))
        q()
    }
}else{
    mycol = c(119,132,147,454,89,404,123,529,463,552,28,54,84,100,558,43,652,31,610,477,256,588,99,81,503,104,562,76,96,495)
    mycol = colors()[rep(mycol,20)]
}
color = mycol[1:nrow(data_for_plot)] 


ggbarplot_plot = data_for_plot %>% 
                 mutate(anno_id = factor(rownames(.), levels = rownames(.)) ) %>% 
                 pivot_longer(cols = all_of(samples_with_order), names_to = 'sample', values_to = 'abundance') %>% # 检查一下，sample列是否是factor，且按照需要的顺序排列
                 ggbarplot(x = 'sample', y = 'abundance', fill = 'anno_id', color = "anno_id", legend = "right", position = position_stack(reverse = TRUE), legend.title = '', order = samples_with_order ) +   # 柱状图; , font.ytickslab = c(14) 修改样本名大小
                 rotate() +  # 旋转
                 scale_fill_manual(values = color[1:nrow(data_for_plot)]) + scale_color_manual(values = color[1:nrow(data_for_plot)]) + # 指定填充颜色
                 theme(axis.line.y = element_blank(), axis.ticks.y = element_blank() ) + xlab('') + ylab(lab) # 取消y轴刻度、title
 

# 进化树 + 柱状图 (样本名永远不会与其他图有交集)
width = 8 + 0.03 * length(samples_with_order)
height = 0.20 * length(samples_with_order)
if(height < 8) height = 8

name = unlist(strsplit(taxonfile,"/"))
pdf(paste(strsplit(name[length(name)], "\\.\\w+\\.txt"), "Community.TreeBar.pdf", sep = "."), width = width, height = height)
plot_grid(ggtree_plot, ggbarplot_plot, nrow = 1, rel_widths = c(0.5,4), align = 'h')
dev.off()

 
