#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: plot_heatmap.r  -i <file> -o <pdf file> --sample_group <file> [--reverse_sample_anno --row_rename_file <file> --break_min <numeric> --break_max <numeric> --scale <string> --gene_list <string> --title <string> --rm_legend --width <int>  --height <int> --ann_colors <string> --display_number --cluster_row --cluster_col --show_row --show_col --rm_annotation_legend --rlib <dir>]

Options:
    -i, --input <file>              输入文件，第一列为基因名，第二列及之后为每一个样本的表达量
    -o, --output <pdf file>         输出pdf文件路径,示例：./a.pdf
    --sample_group <file>           样本分组文件，第一列样本名，第二列及之后的列都是样本分组，包含表头。仅对该文件中包含的样本绘制热图。可以只有第一列样本信息，不给样本分组。如果样本有分组，且分组是空值，空值也被认为是一个分组符号。
    --gene_list <string>            绘制的基因列表，用“逗号”分隔，默认全部绘制 [default: NA]
    --scale <string>                归一化方式，none/row/column [default: none]
    --cluster_row                   进行行聚类
    --cluster_col                   进行列聚类
    --show_row                      显示行名
    --show_col                      显示列名
    --display_number                热图方框里显示数字
    --ann_colors <string>           sample_group第二列分组对应的颜色设定,数量务必等于该列分组数量,示例： red,blue [default: NA]
    --rm_legend                     去掉图例, scale图例
    --rm_annotation_legend          去掉样本分组图例
    --break_min <numeric>           颜色范围划分，最小值设定 [default: NA]
    --break_max <numeric>           颜色范围划分，最大值设定 [default: NA]
    --width <int>                   pdf宽度 [default: 7]
    --height <int>                  pdf高度 [default: 7]
    --title <string>                标题    [default: heatmap plot]
    --row_rename_file <file>        在显示热图时，对热图的行名进行重命名。两列数据，第一列基因名，要与input保持一致，且数量一致；第二列是新的名称，允许重复。包含表头。 [default: NA]
    --reverse_sample_anno           反转样本注释信息顺序。默认 sample_group 最后一列注释内容显示在最上方，但是有时候不是很好看，可以通过这个参数给它反过来。
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，热图\n')
input              <- opts$input
output             <- opts$output
sample_group       <- opts$sample_group
gene_list          <- opts$gene_list
cluster_row        <- opts$cluster_row
cluster_col        <- opts$cluster_col
show_row           <- opts$show_row
show_col           <- opts$show_col
scale              <- opts$scale
display_number     <- opts$display_number
width              <- as.integer(opts$width)
height             <- as.integer(opts$height)
title              <- opts$title
ann_colors         <- opts$ann_colors
rm_legend          <- opts$rm_legend
rm_annotation_legend          <- opts$rm_annotation_legend
break_min          <- opts$break_min
break_max          <- opts$break_max
row_rename_file    <- opts$row_rename_file
reverse_sample_anno    <- opts$reverse_sample_anno
rlib               <- opts$rlib
 

# 加载R包
message('加载pheatmap')
.libPaths(rlib)
library(pheatmap, quietly = TRUE)
library(RColorBrewer)

# 读入数据
message('读入数据')
data       <- read.table(input, sep='\t', header = TRUE, row.names = 1, check.names=FALSE, quote = "", comment.char = "") # 读取数据，行为样本
# 读入分组
data_group <- read.table(sample_group, sep='\t', header = TRUE, row.names = 1, check.names=FALSE, colClasses = 'character')  # 读入分组
sample_names = rownames(data_group)
choose_gene   = rownames(data)

# 读入行重命名标签
labels_row = NULL
if(row_rename_file != 'NA')
{
    data_row_rename <- read.table(row_rename_file, sep='\t', header = TRUE, row.names = 1, check.names=FALSE, stringsAsFactors=F) 
    labels_row = data_row_rename[, 1]
    names(labels_row) = rownames(data_row_rename)
}

if (sum(sample_names %in% colnames(data)) != length(sample_names) )
{   
    lost_samples = sample_names[!sample_names %in% colnames(data)]
    message("[Error] 部分输入的样本名在输入文件中没有找到 : ", lost_samples)
    q()
}

if(gene_list != 'NA')
{
    choose_gene = unlist(strsplit(gene_list, ','))
    if(sum(choose_gene %in% rownames(data)) != length(choose_gene))
    {
        lost_genes = choose_gene[!choose_gene %in% rownames(data)]
        message("[Error] 部分输入的基因名在输入文件中没有找到 : ", lost_genes)
        q()
    }
}
# 定义颜色/且有分组
annotation_colors = NA
if(ann_colors != 'NA' & ncol(data_group) > 0)  
{
    define_colors = c(unlist(strsplit(ann_colors, ',')))
    if(length(define_colors) != length(unique(data_group[, 1])))
    {
        message("[Error] 自定义的分组颜色数量不等于实际分组数量")
        q() 
    }
    names(define_colors) = unique(data_group[, 1])  # 必须添加组名，进行对应
    annotation_colors = list(Group = define_colors)  # 制作pheatmap所要求的分组设定文件
    names(annotation_colors)[1] = colnames(data_group)[1]  # 保证列名一致 
}else if (ncol(data_group) > 0) {  # 只有分组，用默认颜色
    # 颜色模版
    mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)  
    mycol <- colors()[rep(mycol, 50)] 
    group_count = length(unique(data_group[, 1]))
    define_colors = mycol[1:group_count]
    names(define_colors) = unique(data_group[, 1])  # 必须添加组名，进行对应
    annotation_colors = list(Group = define_colors)
    names(annotation_colors)[1] = colnames(data_group)[1]  # 保证列名一致 
}else{
    annotation_colors = NA  # 没有颜色
}
# 开始绘图
# 颜色模版,现在临时固定使用color_navy,后期再个性化调整
color_default     = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(200)
color_navy        = colorRampPalette(c("navy", "white", "firebrick3"))(200)
color_self_define = colorRampPalette(c("#0072c4", "white", "red"))(200)

message('开始绘图')
data_plot = data[choose_gene ,sample_names]  # 提取绘图数据


# 样本注释
annotation_group = NA
if(ncol(data_group) > 0)
{   
    annotation_group = data_group;
    annotation_group[,1] = factor(annotation_group[,1], levels = unique(annotation_group[,1]))
    
    # 样本注释顺序反过来
    if(reverse_sample_anno)
    {
        annotation_group = annotation_group[, rev(colnames(annotation_group)), drop=F]
    }
}

# breaks
breaks = NA
if(break_min != 'NA' & break_max != 'NA')
{   
    breaks = seq(as.numeric(break_min), as.numeric(break_max), length.out = 200)
}

# 调整label rows的顺序，与data_plot保持一致
if( !is.null(labels_row) ) labels_row = labels_row[rownames(data_plot)]

pdf(file=output, width=width, height=height)

pheatmap(data_plot,
      cluster_cols = cluster_col,
      cluster_rows = cluster_row,
      show_colnames = show_col,
      show_rownames = show_row,
      annotation_col = annotation_group,
      color = color_navy,
      scale = scale,
      display_numbers = display_number,
      annotation_colors = annotation_colors,
      main = title,
      legend = !(rm_legend),
      annotation_legend = !(rm_annotation_legend),
      breaks = breaks,
      labels_row = labels_row,
      silent = FALSE
) 
dev.off()


# 自定义数字显示方式 display_numbers = matrix(ifelse(test > 5, "*", "")
# 设定每个格子的宽度、高度 cellwidth = 15, cellheight = 12
# 去掉注释图例 annotation_legend = FALSE
# 对图例的颜色进行设定 annotation_colors = ann_colors，制作方式如下
# ann_colors = list(
#  Time = c("white", "firebrick"),
#  CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
#  GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
# )
# 把热图分块切割显示 gaps_row = c(10, 14) 在第10和14行处添加gap, 要求对行不进行聚类
# cutree_col = 2参数将列按聚类树的结果分成两部分, 要求对列进行聚类
# 自定义行的标签名 labels_row = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "Il10", "Il15", "Il1b")
# order_row = aa$tree_row$order  #记录热图的行排序
# order_col = aa$tree_col$order    #记录热图的列排序
# legend_breaks参数设定图例显示范围，legend_labels参数添加图例标签 legend_breaks = c(1:5), legend_labels = c("1.0","2.0","3.0","4.0","5.0")
