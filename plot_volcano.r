#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: plot_volcano.r  -i <file> -o <pdf file> [--x_lab <string> --y_lab <string> --title <string> --rlib <dir> --xmin <int> --xmax <int> --ymin <int> --ymax <int>]

Options:
    -i, --input <file>        输入文件，第一列为基因名，第二列为log2FoldChange,第三列为Pvalue,第四列为分组(推荐关键字：Up,Down, Not DEG),包含表头
    -o, --output <pdf file>   输出pdf文件路径,示例：./a.pdf
    --x_lab <string>          x轴名称 [default: log2(fold change)]
    --y_lab <string>          y轴名称 [default: -log10(p value)]
    --title <string>          标题    [default: volcano plot]
    --xmin <int>              x轴显示最小值    [default: -10]
    --xmax <int>              x轴显示最大值    [default: 10]
    --ymin <int>              y轴显示最小值    [default: 0]
    --ymax <int>              y轴显示最大值    [default: NA ]
    --rlib <dir>              R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，小提琴图\n')
input             <- opts$input
output            <- opts$output
xmin              <- as.integer(opts$xmin)
xmax              <- as.integer(opts$xmax)
ymin              <- as.integer(opts$ymin)
ymax              <- as.integer(opts$ymax)
x_lab             <- opts$x_lab
y_lab             <- opts$y_lab
title             <- opts$title
rlib              <- opts$rlib
 

# 加载R包
message('加载ggpubr')
.libPaths(rlib)
library(ggpubr, quietly = TRUE)

# 读入数据
message('读入数据')
data <- read.table(input, sep='\t', header = TRUE, check.names=FALSE) # 读取数据，行为样本
colnames(data) <- c('gene', 'fc', 'pvalue', 'group')
data$pvalue_log = -log10(data$pvalue)

# 开始绘图
message('开始绘图')
pdf(output)
# tiff(filename = output)
p <- ggplot(data, aes(x = fc, y = pvalue_log, color = group)) + 
    geom_point() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(x = x_lab, y = y_lab, title = title) +
    scale_x_continuous(limits=c(xmin, xmax)) +
    scale_y_continuous(limits=c(ymin, ymax)) +
    scale_colour_manual(values=c('green', 'black', 'red'))
print(p)
dev.off()
