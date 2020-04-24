#!/home/genesky/software/r/3.5.1/bin/Rscript

.libPaths("/home/genesky/software/r/3.5.1/lib64/R/library/")
library(docopt)
"Usage: pca_plot.r --pca_file <file> --group_file <file> --output <html>  [--size <int> --alpha <float>] 
Options:
   --pca_file <file>           PCA坐标文件，第一列样本名，第二、三列分别是pc1/pc2坐标。含有表头。
   --group_file <file>         分组文件，含有表头。最多三列，可以只放入两列。第一列样本名，第二列分组1，图中用颜色区分。第三列分组2，图中用形状区分。样本顺序没有要求，但是必须包含pca_file中的所有样本
   --size <int>                点的大小 [default: 1]
   --alpha <float>             填充颜色透明度 0-1取值，0完全透明，1不透明。 [default: 1]
   --output <html>             html输出文件， 示例： ./pca.html" -> doc

opts                     <- docopt(doc, version = 'Program : plot pca v1.0 \n          甘斌 129\n')
pca_file                 <- opts$pca_file
group_file               <- opts$group_file
output                   <- opts$output
size                     <- as.numeric(opts$size)
alpha                    <- as.numeric(opts$alpha)
 
# # 测试用参数
# pca_file                 <- '/home/ganb/work/tmp/17B1116A_18B0511B/18B0511B/output/sample_pca.txt'
# group_file               <- '/home/ganb/work/tmp/17B1116A_18B0511B/18B0511B/group.txt'
# output                   <- './test.pdf'

library(plotly)  

pca_data   = read.table(pca_file, head = T, check.names = F, stringsAsFactors=F, row.names = 1)
group_data = read.table(group_file, head = T, check.names = F, stringsAsFactors=F, row.names = 1)

group_col_names = colnames(group_data)
colnames(pca_data)[1:2] = c('PC1', 'PC2')
colnames(group_data)[1] = c('Group')

if(ncol(group_data) == 2) colnames(group_data)[2] = 'Group2'


# 检查分组信息是否缺失
if(sum(rownames(pca_data) %in% rownames(group_data)) != nrow(pca_data))
{
    message("PCA文件中部分样本的分组信息缺失，请仔细确认")
    q()
}


# 准备绘图数据
group_data = group_data[rownames(pca_data), , drop = FALSE]  
data_plot = data.frame(pca_data, group_data, Sample=rownames(pca_data), Group2=paste0('a', group_data$Group), stringsAsFactors=F)

# 形状模版
myshape <- c('circle', 'square', 'diamond', 'cross','triangle', 'hexagon')  
if(ncol(group_data) == 2)
{
    data_plot$Shape = data_plot$Group2
    group2_unique = unique(data_plot$Shape)
    for(col in 1:length(group2_unique))
    {
        data_plot$Shape[data_plot$Shape == group2_unique[col]] = I(myshape[col])
    }
}

# 颜色板
color_map = c('darkmagenta', 'red', 'green', 'yellow', 'cyan', 'blue')

# 绘图
if(ncol(group_data) == 1) p <- plot_ly(data_plot, type = 'scatter', x = ~PC1, y = ~PC2, color =~Group,                 text =~Sample, size = size, alpha = alpha, colors = color_map)
if(ncol(group_data) == 2) p <- plot_ly(data_plot, type = 'scatter', x = ~PC1, y = ~PC2, color =~Group, symbol =~Shape, text =~Sample, size = size, alpha = alpha, colors = color_map)
htmlwidgets::saveWidget(p, output, selfcontained = T)
