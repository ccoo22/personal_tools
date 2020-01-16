#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: dotplot.r  -i <file> --pdf <file> [--Rlib <dir> --width <int>]
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
bin_width = (max(data[, 2]) - min(data[, 2])) / 40 # 设定点的大小
if(bin_width == 0) bin_width = 0.00001
        
p <- ggdotplot(data, x=col_names[1], y=col_names[2], fill = col_names[1], 
          binwidth = bin_width,
          palette = mycol, #杂志nature的配色
          # palette = "npg", #杂志nature的配色
          # palette = "aaas", #杂志Science的配色
          # palette = "jco", #按jco杂志配色方案
          add = c("boxplot"),  # 增加boxplot
          add.params = list(color = col_names[1]),
          title = 'dotplot'
          )  + ylab(col_names[2]) + guides(fill = FALSE, color = FALSE) + theme(plot.title = element_text(hjust = 0.5))  # 去掉legend,标题居中  
print(p)

dev.off()

