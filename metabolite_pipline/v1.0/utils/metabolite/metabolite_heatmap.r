#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)
"Usage: plot_heatmap.r  -i <file> -o <pdf file>  --sample_group_file <file> [--id_file <file>- --scale <string>  --title <string> --rm_legend --width <int>  --height <int> --ann_colors <string> --display_number --cluster_row --cluster_col --show_row --show_col  --rlib <dir>]

Options:
    -i, --input <file>              输入文件，含有特征值和样本表达量的原始文件
    -o, --output <pdf file>         输出pdf文件路径,示例：./a.pdf
    --sample_group_file <file>      样本分组文件，第一列：样本名，第二列：分组名
    --id_file <file>                需要绘制的特征值ID文件[default: all]
    --scale <string>                归一化方式，none/row/column [default: none]
    --cluster_row                   进行行聚类
    --cluster_col                   进行列聚类
    --show_row                      显示行名
    --show_col                      显示列名
    --display_number                热图方框里显示数字
    --ann_colors <string>           分组颜色设定,数量务必等于分组数量,示例： \'group1,red;group2,blue\' [default: NA]
    --rm_legend                     去掉图例
    --width <int>                   pdf宽度 [default: 10]
    --height <int>                  pdf高度 [default: auto]
    --title <string>                标题    [default: heatmap plot]
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version = 'Program : metabolite heatmap analysis v1.0 \n          朱喜 320\n')
input              <- opts$input
output             <- opts$output
groupfile          <- opts$sample_group_file
id                 <- opts$id_file
scale              <- opts$scale
cluster_row        <- opts$cluster_row
cluster_col        <- opts$cluster_col
show_row           <- opts$show_row
show_col           <- opts$show_col
display_number     <- opts$display_number
width              <- as.integer(opts$width)
height             <- opts$height
title              <- opts$title
ann_colors         <- opts$ann_colors
rm_legend          <- opts$rm_legend
rlib               <- opts$rlib
 
# input              <- '/home/zhuxi/reserch/metabolite/MetaboAnalystR/script/normalized_data.txt'
# output             <- '/home/zhuxi/reserch/metabolite/MetaboAnalystR/script/tmp/heatmap.pdf'
# groupfile          <- '/home/zhuxi/reserch/metabolite/MetaboAnalystR/script/groupfile.txt'
# id                 <- '/home/zhuxi/reserch/metabolite/MetaboAnalystR/script/id.txt'
# scale              <- 'none'
# cluster_row        <- TRUE
# cluster_col        <- TRUE
# show_row           <- TRUE
# show_col           <- TRUE
# display_number     <- FALSE
# width              <- as.integer('8')
# height             <- as.integer('8')
# title              <- 'heatmap plot'
# ann_colors         <- 'case,orange;control,green'
# rm_legend          <- FALSE
# rlib               <- "/home/genesky/software/r/3.5.1/lib64/R/library"

# 加载R包
message("start metabolite_heatmap_plot.r")

.libPaths(rlib)
library(pheatmap)
library(RColorBrewer)

# 读入数据
data            <- read.table(input, sep = '\t', header = T, row.names = 1, check.names = F, quote = "") # 读取数据，行为样本
sample.groups   <- read.table(groupfile, header = F, sep = "\t", quote = "", colClasses = 'character')

samples          = as.character(sample.groups[, 1])
sample_group     = as.character(sample.groups[, 2])
feature_id       = rownames(data)

# 数据检测和准备

if (sum(samples %in% colnames(data)) != length(samples) )
{   
    lost_samples = samples[!samples %in% colnames(data)]
    message("[Error] 分组文件中，部分样本编号在丰度文件中不存在! 请仔细检查分组文件 : ", lost_samples)
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
if(ann_colors != 'NA')
{
    define_colors = unlist(strsplit(unlist(strsplit(ann_colors, ';')), ','))
    if(length(define_colors) != 2 * length(unique(sample_group)))
    {
        message("[Error] 自定义的分组颜色数量与实际分组数量不符")
        q()
    }
    define_colors    = t(matrix(define_colors, nrow = 2))
    mycol            = define_colors[, 2]
    names(mycol)     = define_colors[, 1] # 必须添加组名，进行对应
    ann_colors = list(Group = mycol)  # 制作pheatmap所要求的分组设定文件

}else {
    # 颜色模版
    mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)  
    mycol <- colors()[rep(mycol, 50)] 
    group_count = length(unique(sample_group))
    define_colors = mycol[1:group_count]
    names(define_colors) = unique(sample_group)  # 必须添加组名，进行对应
    ann_colors = list(Group = define_colors)
}

################################################# 开始绘图 #################################################

# 颜色模版,使用color_navy
color_default     = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
color_navy        = colorRampPalette(c("navy", "white", "firebrick3"))(100)
color_self_define = colorRampPalette(c("#0072c4", "white", "red"))(100)
data_plot = data[feature_id ,samples]  # 提取绘图数据

# 样本注释
annotation_group <- data.frame(Group = factor(sample_group, levels = unique(sample_group)))
rownames(annotation_group) <- samples
# 
row_name <- as.matrix(data[feature_id, 1])
row_name <- substr(row_name, 1, 24)

if(height == 'auto')
{
    height = 0.2 * nrow(data_plot) + 1  # 设定热图pdf高度
}else{
    height = as.numeric(height)
}

pdf(file = output, width = width, height = height)

pheatmap(data_plot,
    cluster_cols = cluster_col,
    cluster_rows = cluster_row,
    show_colnames = show_col,
    show_rownames = show_row,
    labels_row    = row_name,
    annotation_col = annotation_group,
    color = color_navy,
    fontsize = 10,
    cellheight = 12,
    scale = scale,
    display_numbers = display_number,
    annotation_colors = ann_colors,
    main = title,
    legend = !(rm_legend),
    silent = FALSE
) 
dev.off()
