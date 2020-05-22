#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: plot_volcano.r  -i <file> -o <pdf file> [--position <string> --plot_name <string> --color_name <string> --x_lab <string> --y_lab <string> --title <string> --rlib <dir> --horizontal --sort_value <string> --width <int> --height <int> --label_size <int>]

Options:
    -i, --input <file>        输入文件，第一列为数据名称,第二列及后面的列为要绘制的数据，默认使用第一二列数据绘制，包含表头
    -o, --output <pdf file>   输出pdf文件路径,示例：./a.pdf
    --plot_name <string>      指定要绘制输入文件中的列名,默认绘制第二列 [default: None]
    --color_name <string>     指定分组颜色对应的列名，默认统一红色 [default: None]
    --x_lab <string>          x轴名称 [default: feature]
    --y_lab <string>          y轴名称 [default: None]
    --title <string>          标题    [default:  barplot]
    --width <int>             pdf宽度    [default:  7]
    --height <int>            pdf高度    [default:  7]
    --label_size <int>        文字大小   [default:  12]
    --position <string>       当指定color_name是，x轴会遇到柱子叠加显示的问题， 支持 stack/split 两种模式 [default: split]
    --horizontal              水平绘图
    --sort_value <string>     绘图数据排序, 默认不排序，None/Ascend/Descend [default: None]
    --rlib <dir>              R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，小提琴图\n')
input             <- opts$input
output            <- opts$output
plot_name         <- opts$plot_name
color_name        <- opts$color_name
x_lab             <- opts$x_lab
y_lab             <- opts$y_lab
title             <- opts$title
horizontal        <- opts$horizontal
sort_value        <- opts$sort_value
position          <- opts$position
width             <- as.integer(opts$width)
height            <- as.integer(opts$height)
label_size        <- as.integer(opts$label_size)
rlib              <- opts$rlib
 
direction = 'vertical'
if (horizontal) direction = 'horizontal'

# 加载R包
message('加载ggpubr')
.libPaths(rlib)
library(ggpubr, quietly = TRUE)

# 读入数据
message('读入数据')
data <- read.table(input, sep='\t', header = TRUE, check.names=FALSE) # 读取数据，行为样本
colnames(data)[1] <- c('plot_feature')

# 检查绘图信息
fill_color = 'red'
if(plot_name == 'None') plot_name   = colnames(data)[2]
if(color_name != 'None') fill_color = color_name
if(y_lab     == 'None') y_lab       = plot_name

# 提取绘图数据、并修改列名（特殊列名会导致绘图出错，所以需要重命名）
data_plot = data[, c('plot_feature', plot_name)]
if(color_name != 'None')  data_plot[[color_name]] = data[, color_name]  # 把分组颜色标记信息加入
colnames(data_plot)[2] = 'value'

# 固定绘图顺序
if(sort_value == 'Ascend')   data_plot <- data_plot[sort(data_plot$value, decreasing = FALSE, index.return = T)$ix, ]   # 升序
if(sort_value == 'Descend')  data_plot <- data_plot[sort(data_plot$value, decreasing = TRUE, index.return = T)$ix, ]    # 降序

data_plot$plot_feature = factor(data_plot$plot_feature, levels = unique(data_plot$plot_feature))
 
# 开始绘图
message('开始绘图：', plot_name)
# 设定图像宽高，以保证字体不会重叠

# feature_width = 0.2 * nrow(data_plot)  # 预估占用PDF宽度, 暂定每一个feature占用0.2
# if (direction == 'vertical') width = feature_width
# if (direction == 'horizontal') height = feature_width

pdf(output, width = width, height = height)

if(position == 'split')
{
    p <- ggbarplot(data_plot, x = "plot_feature", y = 'value',
        fill = fill_color, color = "white", orientation = direction, xlab = x_lab, ylab = y_lab, position = position_dodge(0.9)
    ) + theme_pubr(base_size = label_size)
}else {
    p <- ggbarplot(data_plot, x = "plot_feature", y = 'value',
        fill = fill_color, color = "white", orientation = direction, xlab = x_lab, ylab = y_lab
    ) + theme_pubr(base_size = label_size)
}


print(p)
dev.off()
