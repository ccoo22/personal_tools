#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: plot_kegg_dotplot.r  -i <file> -o <pdf file>  [--rlib <string>]

Options:
    -i, --input <file>              输入文件，务必使用kegg原始的输出结果。每一行是一个通路，必须含有表头：Description、GeneRatio、p.adjust、Count
    -o, --output <pdf file>         输出pdf文件路径,示例：./a.pdf
    --rlib <string>                 r包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，自己绘制kegg点状图\n')
input              <- opts$input
output             <- opts$output
rlib               <- opts$rlib
 
.libPaths(rlib)
library(ggplot2)

message("读入数据")
data = read.table(input, head = T, row.names = 1, sep = "\t", stringsAsFactors=F)
if(sum(! c('Description', 'GeneRatio', 'p.adjust', 'Count') %in% colnames(data)) > 0)
{
    message("需要的Description、GeneRatio、p.adjust、Count列信息丢失，请仔细检查，或修改脚本")
    q()
}

message("数据预处理")
data$GeneRatio = unlist(lapply(as.character(data[, 'GeneRatio']), function(x){ res = as.integer( unlist(strsplit(x, split="/")) ); res[1]/res[2] })) # 计算比例
data = data[order(data$GeneRatio), ]  # 从小到大排序
data$Description = factor(data$Description, data$Description)  # 转换成factor格式，从而确保绘图顺序

pdf(output, width = 10)
p <- ggplot(data, aes(x=Description, y=GeneRatio)) + 
     geom_point(aes(size=Count, color=p.adjust)) + 
     scale_color_continuous(low="red", high="blue", name = 'p.adjust', guide=guide_colorbar(reverse=TRUE)) + 
     xlab(NULL) + scale_size(range=c(3, 8)) + 
     guides(color = guide_colourbar(order=1, reverse=TRUE), size = guide_legend(order=2)) +  # 固定legend显示顺序
     coord_flip()
print(p)

dev.off()

 