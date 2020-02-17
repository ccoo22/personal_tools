#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: dotplot.r  -i <file> --pdf <file> [--Rlib <dir> --width <int> --height <int> --title <string> ]
Options:
   -i, --input <file>    输入文件，两列或三列，包含表头，第一列x，第二列y，第三列分组，如果没有分组信息，第三列可不提供
   --pdf <file>          输出pdf文件路径
   --width <int>         pdf宽度 [default: 7]
   --height <int>        pdf高度 [default: 7]
   --title <string>      图片标题 [default: scatter]
   --Rlib <dir>          R包路径  [default: /home/genesky/software/r/3.5.1/lib64/R/library] " -> doc

opts                     <- docopt(doc, version = '点图绘制 \n          甘斌 129\n')
input                    <- opts$input
output_pdf               <- opts$pdf
width                    <- opts$width
height                   <- opts$height
title                    <- opts$title
Rlib                     <- opts$Rlib
 
.libPaths(Rlib)

# 加载R包
library(ggpubr, quietly = TRUE)

# 读入数据
data <- read.table(input, sep='\t', header = TRUE, check.names=FALSE) # 读取数据，行为样本
col_names <- colnames(data)

set_color = "black"
if(length(col_names) == 3) set_color = col_names[3]
# pdf
pdf(output_pdf, width = as.numeric(width), height = as.numeric(height)) # 一组2.5宽度
    
p <- ggscatter(data, x=col_names[1], y=col_names[2], color = set_color, 
          palette = "npg", #杂志nature的配色
          # palette = "aaas", #杂志Science的配色
          # palette = "jco", #按jco杂志配色方案
          title = title
          )  + ylab(col_names[2])  + theme(plot.title = element_text(hjust = 0.5))  # 去掉legend,标题居中  
print(p)

dev.off()

