#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: plot_heatmap.r  -i <file> -o <pdf file> --sample_list <string> --group_list <string> [--scale <string> --gene_list <string> --title <string> --rm_legend --width <int>  --height <int> --ann_colors <string> --display_number --cluster_row --cluster_col --show_row --show_col  --rlib <dir>]

Options:
    -i, --input <file>              输入文件，第一列为基因名，第二列及之后为每一个样本的表达量
    -o, --output <pdf file>         输出pdf文件路径,示例：./a.pdf
    --sample_list <string>          样本列表，用“逗号”分隔,例如：C1,C2,C3,C4,S1,S2,S3,S5
    --group_list <string>           样本分组列表，与样本列表中的样本名一一对应，用“逗号”分隔，例如：C,C,C,C,S,S,S,S
    --gene_list <string>            绘制的基因列表，用“逗号”分隔，默认全部绘制 [default: NA]
    --scale <string>                归一化方式，none/row/column [default: none]
    --cluster_row                   进行行聚类
    --cluster_col                   进行列聚类
    --show_row                      显示行名
    --show_col                      显示列名
    --display_number                热图方框里显示数字
    --ann_colors <string>           分组颜色设定,数量务必等于分组数量,示例： red,blue [default: NA]
    --rm_legend                     去掉图例
    --width <int>                   pdf宽度 [default: 7]
    --height <int>                  pdf高度 [default: 7]
    --title <string>                标题    [default: heatmap plot]
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，热图\n')
input              <- opts$input
output             <- opts$output
sample_list        <- opts$sample_list
group_list         <- opts$group_list
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
rlib               <- opts$rlib
 

# 加载R包
message('加载pheatmap')
.libPaths(rlib)
library(pheatmap, quietly = TRUE)
library(RColorBrewer)

# 读入数据
message('读入数据')
data <- read.table(input, sep='\t', header = TRUE, row.names = 1, check.names=FALSE) # 读取数据，行为样本
sample_names = unlist(strsplit(sample_list, ','))
sample_groups = unlist(strsplit(group_list, ','))
choose_gene   = rownames(data)



if (length(sample_names) != length(sample_groups))
{
    message("[Error] 输入的样本数量与分组数量不一致")
    q()
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
if(ann_colors != 'NA')
{
    define_colors = c(unlist(strsplit(ann_colors, ',')))
    if(length(define_colors) != length(unique(sample_groups)))
    {
        message("[Error] 自定义的分组颜色数量不等于实际分组数量")
        q() 
    }
    names(define_colors) = unique(sample_groups)  # 必须添加组名，进行对应
    ann_colors = list(Group = define_colors)  # 制作pheatmap所要求的分组设定文件
}else {
    # 颜色模版
    mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)  
    mycol <- colors()[rep(mycol, 50)] 
    group_count = length(unique(sample_groups))
    define_colors = mycol[1:group_count]
    names(define_colors) = unique(sample_groups)  # 必须添加组名，进行对应
    ann_colors = list(Group = define_colors)
}
# 开始绘图
# 颜色模版,现在临时固定使用color_navy,后期再个性化调整
color_default     = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
color_navy        = colorRampPalette(c("navy", "white", "firebrick3"))(100)
color_self_define = colorRampPalette(c("#0072c4", "white", "red"))(100)

message('开始绘图')
data_plot = data[choose_gene ,sample_names]  # 提取绘图数据


# 样本注释
annotation_group <- data.frame(Group = factor(sample_groups, levels = unique(sample_groups)))
rownames(annotation_group) <- sample_names
# 
 
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
      annotation_colors = ann_colors,
      main = title,
      legend = !(rm_legend),
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