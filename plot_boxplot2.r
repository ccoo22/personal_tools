#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: boxplot.r  -i <file> --pdf <file> [--add_diff --Rlib <dir> --width <int> --y_name <string> ]
Options:
   -i, --input <file>    输入文件，至少三列，包含表头，第一列样本名，第二列为数值型数据，第三列及之后的列为分组数据，每一列对应一张图
   --pdf <file>          输出pdf文件路径
   --add_diff            增加差异分析
   --Rlib <dir>          R包路径  [default: /home/genesky/software/r/3.5.1/lib64/R/library] " -> doc

opts                     <- docopt(doc, version = '点图绘制 \n          甘斌 129\n')
input                    <- opts$input
output_pdf               <- opts$pdf
add_diff                 <- opts$add_diff
Rlib                     <- opts$Rlib
 
.libPaths(Rlib)

# 加载R包
library(ggpubr, quietly = TRUE)

# 读入数据
data <- read.table(input, sep='\t', header = TRUE, check.names=FALSE) # 读取数据，行为样本
colnames(data)[1] = 'SampleName'
value_name = colnames(data)[2]
class_names = colnames(data)[3:ncol(data)]

# pdf
pdf(output_pdf, width = 7, height = 7) 

for(class_name in class_names)
{
    data_plot = data[, c('SampleName', value_name, class_name)]
    data_plot = data_plot[complete.cases(data_plot), ]  # 去掉缺失


    message("plot : ", class_name)
    p <- ggboxplot(data_plot, x=class_name, y=value_name, fill = class_name, 
          palette = "npg", #杂志nature的配色
          # palette = "aaas", #杂志Science的配色
          # palette = "jco", #按jco杂志配色方案
          add = c("jitter"),  # 增加点
          add.params = list(color = 'black'),
          outlier.size = -1, # 隐藏离群点，否则离群点会在图上出现两次（boxplot绘制离群点，jitter再绘制所有点，从而使离群点绘制两次）
          title = class_name
          )  + ylab(value_name) + guides(fill = FALSE, color = FALSE) + theme(plot.title = element_text(hjust = 0.5))  # 去掉legend,标题居中  
    if(add_diff){
        p = p + stat_compare_means()
    }
   
    print(p)
}
dev.off()

