#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: plot_volcano.r  -i <file> -o <pdf file> [--x_lab <string> --y_lab <string> --title <string> --rlib <dir> --horizontal --width <int> --height <int>]

Options:
    -i, --input <file>        输入文件，第一列为数据名称,第二列为数字，包含表头
    -o, --output <pdf file>   输出pdf文件路径,示例：./a.pdf
    --x_lab <string>          x轴名称 [default: feature]
    --y_lab <string>          y轴名称 [default: -log10(p value)]
    --title <string>          标题    [default:  barplot]
    --width <int>             pdf宽度    [default:  7]
    --height <int>            pdf高度    [default:  7]
    --horizontal              水平绘图
    --rlib <dir>              R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，小提琴图\n')
input             <- opts$input
output            <- opts$output
x_lab             <- opts$x_lab
y_lab             <- opts$y_lab
title             <- opts$title
horizontal        <- opts$horizontal
width             <- as.integer(opts$width)
height            <- as.integer(opts$height)
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
colnames(data) <- c('feature', 'value')
data$feature = factor(data$feature, levels = data$feature)

# 开始绘图
message('开始绘图')
pdf(output, width = width, height = height)
p <- ggbarplot(data, "feature", "value",
 fill = "red", color = "red", orientation = direction, xlab = x_lab, ylab = y_lab
 )

print(p)
dev.off()
