#!/home/genesky/software/r/4.1.2/bin/Rscript

library(docopt)

"Usage: plot_heatmap_complexheatmap.r  -i <file> -o <pdf file> --sample_group <file> [--ann_colors <string>  --rlib <dir>  --pdf_width <int> --pdf_height <int> --display_number --display_number_size <numeric>  --row_rename_file <file> --col_rename_file <file>  --row_font_size <numeric> --col_font_size <numeric>  --cell_height <numeric> --cell_width <numeric>  --cell_border_width <numeric>  --cell_border_color <string> --title <string>  --rm_legend --gene_list <string> --col_temp <string> --col_temp_defined <string> --legend_breaks <string> --legend_name <string> --scale <string> --cluster_row  --cluster_col --show_row --show_col]

Options:
    -i, --input <file>              输入文件，第一列为基因名，第二列及之后为每一个样本的表达量
    --sample_group <file>           样本分组文件，第一列样本名，第二列及之后的列都是样本分组，包含表头。仅对该文件中包含的样本绘制热图。可以只有第一列样本信息，不给样本分组。如果样本有分组，且分组是空值，空值也被认为是一个分组符号。
                                    注：如果有很多注释列，complexheatmap 会自动把legend换行显示
                                    注：legend的颜色都是随机搭配生成的，每一次运行，配色结果都不相同，除非你主动声明颜色搭配
    --ann_colors <string>           设定sample_group第二列分组对应的颜色，数量务必等于该列分组数量,示例： red,blue,#648c11 
                                    注：该方法可以阻止软件自动配色
                                    注：颜色添加的顺序与sample_group第二列中的分类名称 unique 顺序对应
    -o, --output <pdf file>         输出pdf文件路径,示例：./a.pdf
    --gene_list <string>            绘制的基因列表，用“逗号”分隔，默认全部绘制 [default: NA]
    --col_temp <string>             热图颜色模版, 目前只提供了 2个配色模版，支持： navy_white_red / navy_white_firbrick3 [default: navy_white_firbrick3]
    --col_temp_defined <string>     自定义热图颜色模版, 只接收16进制颜色类型，多个颜色之间用逗号分隔， 例如 #0047ab,#e9967a,#648c11 [default: NA]
                                    当定义了该参数，--col_type 参数会被忽略
    --legend_breaks <string>        控制热图数值与颜色的对应关系（限制数值显示范围），至少填写三个数值，且从小到大，逗号分隔。 例如： '-1,0,1' [default: NA]
    --legend_name <string>          热图legend的名称 [default: none]
    --scale <string>                归一化方式，none/row/column [default: none]


    --cluster_row                   进行行聚类
    --cluster_col                   进行列聚类
    --show_row                      显示行名
    --show_col                      显示列名
    --display_number                热图方框里显示数字
    --display_number_size <numeric> 修改热图方框里字体大小, 例如 10

    --rm_legend                     去掉图例, scale图例
    --pdf_width <int>               pdf宽度 [default: 7]
    --pdf_height <int>              pdf高度 [default: 7]
    --title <string>                标题    [default: heatmap plot]
    --row_rename_file <file>        在显示热图时，对热图的行名进行重命名。两列数据，第一列基因名，要与input保持一致，且数量一致；第二列是新的名称，允许重复。包含表头。
    --col_rename_file <file>        在显示热图时，对热图的列名进行重命名。两列数据，第一列样本名，要与input保持一致，且数量一致；第二列是新的名称，允许重复。包含表头。
    --row_font_size <numeric>       行名字体大小， 例如 10
    --col_font_size <numeric>       列名字体大小， 例如 10
    --cell_border_color <string>    热图边框的颜色， 16进制字符， 例如 #0047ab [default: NA]
    --cell_border_width <numeric>   热图边框的宽度，例如 2 [default: NA]
    --cell_width <numeric>          cell宽度, 例如 0.5 ，控制每个cell的宽度
    --cell_height <numeric>         cell高度, 例如 0.5
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/4.1.2/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，complexheatmap 热图\n')
input              <- opts$input
output             <- opts$output
sample_group       <- opts$sample_group
ann_colors         <- opts$ann_colors
gene_list          <- opts$gene_list
scale              <- opts$scale

col_temp           <- opts$col_temp
col_temp_defined   <- opts$col_temp_defined
legend_breaks      <- opts$legend_breaks
legend_name        <- opts$legend_name

title              <- opts$title

pdf_width          <- as.numeric(opts$pdf_width)
pdf_height         <- as.numeric(opts$pdf_height)

rm_legend          <- opts$rm_legend

cluster_row        <- opts$cluster_row
cluster_col        <- opts$cluster_col
show_row           <- opts$show_row
show_col           <- opts$show_col

display_number     <- opts$display_number
display_number_size      <- opts$display_number_size

row_font_size      <- opts$row_font_size
col_font_size      <- opts$col_font_size
row_rename_file    <- opts$row_rename_file
col_rename_file    <- opts$col_rename_file

cell_width         <- opts$cell_width
cell_height        <- opts$cell_height
cell_border_width  <-  opts$cell_border_width
cell_border_color <- opts$cell_border_color



rlib               <- opts$rlib



# 加载R包
message('加载ComplexHeatmap')
.libPaths(rlib)
library(ComplexHeatmap, quietly = TRUE)
library(circlize)
library(RColorBrewer)
# input =  "/home/zhangshuang/work/other_project/21B0420A/cluster/tcga_brca.fpkm.diff.txt"
# sample_group = "/home/zhangshuang/work/other_project/21B0420A/heatmap/group.txt"

# 读入数据
message('读入数据')
data_raw       <- read.table(input, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "") # 第一列为基因名，第二列及之后为每一个样本的表达量

## 读入分组
data_group <- read.table(sample_group, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "", colClasses = 'character')  # 读入分组

choose_samples = rownames(data_group)
choose_genes   = rownames(data_raw)


## 分组文件是否丢失数据
lost_samples = setdiff(choose_samples, colnames(data_raw))
if (length(lost_samples) > 0 )
{   
    message("[Error] sample_group 中的样本 在 input 中有缺失 : ", paste(lost_samples, collapse=','))
    q()
}

## 选定要分析的基因
if(gene_list != 'NA')
{
    choose_genes = unlist(strsplit(gene_list, ','))
    lost_genes = setdiff(choose_genes, rownames(data_raw))
    if(length(lost_genes) > 0)
    {
        message("[Error] gene_list 中的基因名在input中没有找到 : ", paste(lost_genes, collapse=','))
        q()
    }
}

## 原始热图数据准备
data_choose = data_raw[choose_genes, choose_samples, drop=F]
data_choose = as.matrix(data_choose)
## scale 归一化 处理
message('确认是否归一化数据')
data_plot = data_choose
if(scale == 'row')
{   
    message('    row归一化数据')
    data_plot = t(scale(t(data_choose)))
}
if(scale == 'col')
{
    message('    col归一化数据')
    data_plot = scale(data_choose)
}

## 参考 pheatmap 代码，确定当前绘图的breaks
message('确定热图颜色、数值 breaks')
if (!identical(legend_breaks, 'NA'))
{   
    # 自定义breaks
    message('    使用自定义 breaks')
    legend_breaks = as.numeric(unlist(strsplit(legend_breaks, ',')))
}
## 自动 breaks
if (identical(scale, "row") || identical(scale, "column")) {
    if (identical(legend_breaks, 'NA')) {
        message('    使用自动化 breaks')
        lim = quantile(abs(data_plot), 0.975)
        le = pretty(c(-lim, lim), n = 3)
        if (length(le) == 7 && le[1] == -3) {
            le = c(-3, -1.5, 0, 1.5, 3)
        }
        else if (!0 %in% le) {
            le = c(le[1], le[1]/2, 0, le[length(le)]/2, le[length(le)])
        }
        legend_breaks = le
    }
}

## 热图颜色模版, 注意 -2,0,2 同时起到了限制热图数值显示范围的作用
message('颜色模版确认')
col_fun = NA

if(identical(legend_breaks, 'NA'))
{
    if(col_temp == 'navy_white_red')
    {
        col_fun = colorRampPalette(c("navy", "white", "red"))(200)
    }
    if(col_temp == 'navy_white_firbrick3')
    {
        col_fun = colorRampPalette(c("navy", "white", "firebrick3"))(200)
    } 
    if(col_temp_defined != 'NA')
    {   
        # 自定义配色方案
        col_temp_defined = unlist(strsplit(col_temp_defined,','))
        col_fun = colorRampPalette(col_temp_defined)(200)
    }
}else{
    if(col_temp == 'navy_white_red')
    {   
        col_fun = colorRamp2(legend_breaks, colorRampPalette(c("navy", "white", "red"))(length(legend_breaks)))
    }
    if(col_temp == 'navy_white_firbrick3')
    {
        col_fun = colorRamp2(legend_breaks, colorRampPalette(c("navy", "white", "firebrick3"))(length(legend_breaks)))
    } 
    if(col_temp_defined != 'NA')
    {   
        # 自定义配色方案
        col_temp_defined = unlist(strsplit(col_temp_defined,','))
        col_fun = colorRamp2(legend_breaks, colorRampPalette(col_temp_defined)(length(legend_breaks)))
    }
}

## cell 边框颜色
if(identical(cell_border_color, 'NA'))
{
    cell_border_color = ifelse(nrow(data_plot) <100 & ncol(data_plot) < 100, "grey60", NA)
}
## cell  边框大小
cell_border_width = ifelse(identical(cell_border_width, 'NA'), NA, as.numeric(cell_border_width))

## 列注释（根据 sample_group 文件制作）， 绘图矩阵已经按照sample_group 矩阵排过序了，这里可以直接用
top_annotation = NULL
if(ncol(data_group) > 0)
{   
    if(!is.null(ann_colors))
    {   
        # 自定义配色
        ann_colors = unlist(strsplit(ann_colors, ','))
        # 数量是否与第二列注释相同
        class_second_column = unique(data_group[,1])
        if(length(ann_colors) != length(class_second_column))
        {
            message("[Error] ann_colors 参数中声明的颜色数量与 sample_group 文件的第二列分类数量不符，请重新调整参数: ann_colors = ", ann_colors)
            q()
        }
        color = list(structure(ann_colors, names=class_second_column))
        names(color) = colnames(data_group)[1]
        top_annotation = HeatmapAnnotation(df = data_group, col=color)

    }else{
        # 自动随机配色
        top_annotation = HeatmapAnnotation(df = data_group)
    }
    
}

## 设置热图区域的宽度、高度
width = NULL
if(!is.null(cell_width))
{
    width = unit(as.numeric(cell_width) * ncol(data_plot), "cm")
}

height = NULL
if(!is.null(cell_height))
{
    height = unit(as.numeric(cell_height) * nrow(data_plot), "cm")
}

## 行、列名称字体大小
if(!is.null(row_font_size))
{
    row_font_size = as.numeric(row_font_size)
}
if(!is.null(col_font_size))
{
    col_font_size = as.numeric(col_font_size)
}


## 重设行、列名称
row_labels = rownames(data_plot)
if(!is.null(row_rename_file))
{
    row_rename_data = read.table(row_rename_file, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "", colClasses = 'character') # 两列信息
    # 检查
    lost_rows = setdiff(choose_genes, rownames(row_rename_data))
    if(length(lost_rows) > 0)
    {
        message("[Error] row_rename_file 中必须包含input文件中的所有要分析的基因名称 : ", paste(lost_rows, collapse=','))
        q()
    }
    row_labels = structure(row_rename_data[choose_genes, 1], names = choose_genes)
}

column_labels = colnames(data_plot)
if(!is.null(col_rename_file))
{
    col_rename_data = read.table(col_rename_file, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "", colClasses = 'character') # 两列信息
    # 检查
    lost_cols = setdiff(choose_samples, rownames(col_rename_data))
    if(length(lost_cols) > 0)
    {
        message("[Error] col_rename_file 中必须包含input文件中的所有要分析的样本名称 : ", paste(lost_cols, collapse=','))
        q()
    }
    column_labels = structure(col_rename_data[choose_samples, 1], names = choose_samples)
}


## 热图里显示数字、以及数字大小
if(!is.null(display_number_size))
{
    display_number_size = as.numeric(display_number_size)
}else{
    display_number_size = NA
}

cell_fun = NULL
if(display_number)
{
    cell_fun <- function(j, i, x, y, width, height, fill) 
        { 
            grid.text(sprintf("%.1f", data_plot[i, j]), x, y, gp=gpar(fontsize = display_number_size)) 
        }
}




message('开始绘图')
pdf(file=output, width=pdf_width, height=pdf_height)
Heatmap(matrix = data_plot,
    col = col_fun, # 热图颜色模版
    name = ifelse(legend_name == 'none', ' ', legend_name),  # key legend名称
    cluster_rows = cluster_row, # 行聚类
    cluster_columns = cluster_col, # 列聚类
    show_row_names = show_row, # 显示行名称
    show_column_names = show_col, # 显示列名称
    clustering_distance_rows = "euclidean",  # 行距离计算方法
    clustering_method_rows = "complete",  # 行聚类方法
    clustering_distance_columns = "euclidean",  # 列距离计算方法
    clustering_method_columns = "complete",  # 列聚类方法
    rect_gp = gpar(col = cell_border_color, lwd = cell_border_width), # cell 边框颜色、边框大小
    top_annotation = top_annotation, # 添加列注释
    show_heatmap_legend = !rm_legend, # 去掉热图legend
    column_title = title, # 列标题
    width = width,  # 设置热图区域的宽度： 不包括进化树、legend
    height = height,  # 设置热图区域的高度： 不包括进化树、legend
    row_names_gp = gpar(fontsize = row_font_size),  # 行名颜色、大小调整
    column_names_gp = gpar(fontsize = col_font_size), # 列名颜色、大小调整
    row_labels = row_labels,  # 重设 行名
    column_labels = column_labels, # 重设列名
    cell_fun = cell_fun,  # 热图里添加数字
)
dev.off()

