#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: boxplot.r  -i <file> --pdf <file> [--Rlib <dir> --width <int>]
Options:
   -i, --input <file>    输入文件，两列，包含表头，第一列分组，第二列数值型数据
   --pdf <file>          输出pdf文件路径
   --width <int>         pdf绘图时，每一个分组的宽度 [default: 2.5]
   --Rlib <dir>          R包路径  [default: /home/genesky/software/r/3.5.1/lib64/R/library] " -> doc

opts                     <- docopt(doc, version = '点图绘制 \n          甘斌 129\n')
input                    <- opts$input
output_pdf               <- opts$pdf
width                    <- opts$width
Rlib                     <- opts$Rlib
 
.libPaths(Rlib)

# 加载R包
library(ggpubr, quietly = TRUE)

# 读入数据
data <- read.table(input, sep='\t', header = TRUE, check.names=FALSE) # 读取数据，行为样本
col_names <- colnames(data)

# 分组数量
group_count = length(unique(data[, 1]))


# 颜色模版
mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)  
mycol <- colors()[rep(mycol, 50)] 

# pdf
pdf(output_pdf, width = group_count * as.numeric(width), height = 7) # 一组2.5宽度
        
p <- ggboxplot(data, x=col_names[1], y=col_names[2], color = col_names[1], 
          palette = "npg", #杂志nature的配色
          # palette = "aaas", #杂志Science的配色
          # palette = "jco", #按jco杂志配色方案
          add = c("jitter"),  # 增加点
          add.params = list(color = col_names[1]),
          outlier.size = -1, # 隐藏离群点，否则离群点会在图上出现两次（boxplot绘制离群点，jitter再绘制所有点，从而使离群点绘制两次）
          title = 'boxplot'
          )  + ylab(col_names[2]) + guides(fill = FALSE, color = FALSE) + theme(plot.title = element_text(hjust = 0.5))  # 去掉legend,标题居中  
print(p)

dev.off()

