#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: barplot.r  -i <file> --pdf <file> --sample_group <file> [--y_lab <string>  --add_mark <file>  --width <int> --height <int> --legend_pos <string> --title <string> --position_dodge_value <numeric> --rotate_x_text <numeric> ]
Options:
    -i, --input <file>                   一行对应一个绘图特征，一列对应一个样本。第一列是特征名称，x轴标签，第二列及之后是每个样本的数据。
    --pdf <file>                         the pdf output file 
    --sample_group <file>                样本分组文件，第一列样本名，第二列样本分组，包含表头。仅对该文件中包含的样本绘图
    --add_mark <file>                    在特征的柱状图上方显示pvalue等标记信息，两列数据，第一列特征名，要与input保持一致，且数量一致；第二列是要添加的标签信息。包含表头。当前的版本仅支持2组的情况。 [default: NA]
    --width <int>                        pdf宽度 [default: 7]
    --height <int>                       pdf高度 [default: 7]
    --title <string>                     图片标题 [default: barplot]
    --rotate_x_text <numeric>            旋转x轴字体角度，例如 45 [default: 0]
    --position_dodge_value <numeric>     调整分组之间柱子的距离, 0.7靠在一起，>0.7分离 [default: 0.7]
    --legend_pos <string>                legend显示的位置，top, bottom, left, right, none [default: right]
    --y_lab <string>                     设定y轴标签 [default: Value]
   " -> doc

opts                    <- docopt(doc, version = 'Program : 对一个矩阵，绘制barplot，所有特征放在一张图上。barplot  甘斌 129\n')
input                   <- opts$input
output_pdf              <- opts$pdf
sample_group            <- opts$sample_group
width                   <- as.integer(opts$width)
height                  <- as.integer(opts$height)
position_dodge_value    <- as.numeric(opts$position_dodge_value)
rotate_x_text           <- as.numeric(opts$rotate_x_text)
legend_pos              <- opts$legend_pos
title                   <- opts$title
add_mark                <- opts$add_mark
y_lab                   <- opts$y_lab

library(ggpubr, quietly = TRUE)
library(tidyr, quietly = TRUE)

# 读取数据，行为样本
data <- read.table(input, sep='\t', header = TRUE, check.names=FALSE) 
colnames(data)[1] = 'feature_name'
# 获取样本、分组
data_group <- read.table(sample_group, sep='\t', header = TRUE, row.names = 1, check.names=FALSE, colClasses = 'character')  # 读入分组
sample_names = rownames(data_group)
sample_groups = data_group[,1]


################### 绘图
 
# 提取、转换数据
data_tmp          <- data[, c('feature_name', sample_names), drop = FALSE] # 取出样本与位置信息
data_plot         <- gather(data_tmp, Sample, Value, -feature_name) # 数据转换格式
# 增加样本分组列
data_plot$Group <- NA 
for(sample_col in 1:length(sample_names))
{
    data_plot$Group[ data_plot$Sample == sample_names[sample_col] ] <- sample_groups[sample_col]
}
data_plot$Group = factor(data_plot$Group, levels=unique(data_plot$Group))

# 去掉缺失值
data_plot <- data_plot[complete.cases(data_plot), ]

pdf(file=output_pdf, width=width, height=height)

p <- ggbarplot(data_plot, x="feature_name", y="Value", fill = "Group", 
            # palette = "npg", #杂志nature的配色
            palette = "aaas", #杂志Science的配色
            # palette = "jco", #按jco杂志配色方案
            add = "mean_se", # 取均值，并添加标准差曲线
            error.plot = "upper_errorbar", # 只显示上曲线
            position = position_dodge(position_dodge_value),
            title = title
        ) + ylab(y_lab) + theme(plot.title = element_text(hjust = 0.5)) # 标题居中    

p <- ggpar(p, legend = legend_pos) # legend位置
p <- p + rotate_x_text(rotate_x_text)


# 添加pvalue
if(add_mark != 'NA')
{   
    data_add_mark <- read.table(add_mark, sep='\t', header = TRUE, check.names=FALSE, stringsAsFactors=F) 

    feature = data_add_mark[, 1]
    mark    = data_add_mark[, 2]
    group1  = unique(sample_groups)[1]
    group2  = unique(sample_groups)[2]
    
    # 确定y轴绘图位置
    y.position = c()
    for(gene in feature)
    {
        gene_data = data_plot[data_plot$feature_name == gene, ]
        value1_75     <- quantile (gene_data[gene_data$Group == group1, 'Value' ], .75)[[1]]
        value2_75     <- quantile (gene_data[gene_data$Group == group2, 'Value' ], .75)[[1]]
        pos = max(value1_75, value2_75)
        y.position = c(y.position, pos)
    }
    tip_length = max(y.position) * 0.02
 
    y.position = y.position + max(y.position) * 0.04
    
    add_info = data.frame(feature, mark, group1, group2, y.position)

    p <- p + stat_pvalue_manual(add_info, x='feature', xmin = 'group1', xmax = 'group2', y.position = 'y.position', label = "mark", tip.length = rep(tip_length, nrow(add_info))) 

}

print(p)
dev.off()

