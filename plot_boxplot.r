#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: boxplot.r  -i <file> --pdf <file> [--Rlib <dir> --width <int> --y_name <string> ]
Options:
   -i, --input <file>    输入文件，至少三列，包含表头，第一列样本名，第二列分组，第三列及之后的列为数值型数据，每一列对应一张图
   --pdf <file>          输出pdf文件路径
   --y_name <string>     y轴名称
   --width <int>         pdf绘图时，每一个分组的宽度 [default: 2.5]
   --Rlib <dir>          R包路径  [default: /home/genesky/software/r/3.5.1/lib64/R/library] " -> doc

opts                     <- docopt(doc, version = '点图绘制 \n          甘斌 129\n')
input                    <- opts$input
output_pdf               <- opts$pdf
y_name                   <- opts$y_name
width                    <- opts$width
Rlib                     <- opts$Rlib
 
.libPaths(Rlib)

# 加载R包
library(ggpubr, quietly = TRUE)

# 读入数据
data <- read.table(input, sep='\t', header = TRUE, check.names=FALSE) # 读取数据，行为样本
colnames(data)[c(1,2)] = c('SampleName', 'SampleGroup')

# 分组数量
group_count = length(unique(data[, 'SampleGroup']))

# 基因数量
genes = colnames(data)[3:ncol(data)]


# 颜色模版
mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)  
mycol <- colors()[rep(mycol, 50)] 

# pdf
pdf(output_pdf, width = group_count * as.numeric(width), height = 7) # 一组2.5宽度

for(gene in genes)
{
   data_plot = data[, c('SampleName', 'SampleGroup', gene)]
   colnames(data_plot)[3] = 'GENE'  # 变量不能以数字开头，故替换名称
   data_plot = data_plot[complete.cases(data_plot), ]  # 去掉缺失


   message("plot : ", gene)
   p <- ggboxplot(data_plot, x='SampleGroup', y='GENE', color = 'SampleGroup', 
          palette = "npg", #杂志nature的配色
          # palette = "aaas", #杂志Science的配色
          # palette = "jco", #按jco杂志配色方案
          add = c("jitter"),  # 增加点
          add.params = list(color = 'SampleGroup'),
          outlier.size = -1, # 隐藏离群点，否则离群点会在图上出现两次（boxplot绘制离群点，jitter再绘制所有点，从而使离群点绘制两次）
          title = gene
          )  + ylab(y_name) + guides(fill = FALSE, color = FALSE) + theme(plot.title = element_text(hjust = 0.5))  # 去掉legend,标题居中  
   print(p)
}
dev.off()

