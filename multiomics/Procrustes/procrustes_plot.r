#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)
"Usage: procrustes_plot.r --matrix1 <file> --matrix2 <file> --mid1 <string> --mid2 <string> -o <dir> -g <file> [--color <string> --trans <string> --lib <string>]
Options:
	--matrix1 <file>            矩阵1。每一行是一个样本，每一列是一个特征。
	--matrix2 <file>            矩阵2。每一行是一个样本，每一列是一个特征。
	--mid1 <string>             矩阵1名称。
	--mid2 <string>             矩阵2名称。
	-o, --output <dir>          输出文件路径。
	-g, --groupfile <file>      分组信息，有表头，第一列样本，第二列分组信息。
	--color <string>            指定分组颜色。 [default: #cf5641,#90900a,#87CEFA,#F0E68C]
	--trans <string>            矩阵数据是否需要转置，默认不转置，即每一行是一个样本，每一列是一个特征。可选 TRUE/FALSE [default: FALSE]
	--lib <string>              R包路径lib [default: /home/genesky/database/r/3.5.1/]"-> doc

opts       <- docopt(doc, version = 'procrustes plot 比较两组数据一致性\n')
matrix1    <- opts$matrix1
matrix2    <- opts$matrix2
mid1       <- opts$mid1
mid2       <- opts$mid2
output     <- opts$output
groupfile  <- opts$groupfile
color      <- opts$color
trans      <- opts$trans
lib        <- opts$lib

lib = c(unlist(strsplit(lib, ',')))
.libPaths(lib)

# /home/genesky/software/r/3.5.1/bin/Rscript /home/genesky/pipeline/tools_script/Procrustes/procrustes_plot.r --matrix1 /home/genesky/pipeline/tools_script/Procrustes/example/xq_all.xls --matrix2 /home/genesky/pipeline/tools_script/Procrustes/example/fb_all.xls --mid1 serum --mid2 faecal -g /home/genesky/pipeline/tools_script/Procrustes/example/group.txt -o /home/genesky/pipeline/tools_script/Procrustes/example/

options(scipen=200)
library(vegan)
library(ggplot2)


group <- read.table(groupfile, sep = '\t', header = T, stringsAsFactors = FALSE, check.names = FALSE)
colnames(group) <- c("sample","groups")

xq <- read.table(matrix1, row.names = 1, header = T, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
if (trans == "TRUE") xq = as.data.frame(t(xq)) #转置
# 判断指定的样本是否存在于矩阵数据1
if (length(group$sample[!(group$sample %in% rownames(xq))]) != 0)
{
	cat('[Error] 分组文件中，部分样本编号在文件中不存在! 请仔细检查分组文件 !\n')
	q()
}
xq = xq[group$sample, ]
xq_pca <- rda(xq, scale = TRUE)

fb <- read.table(matrix2, row.names = 1, header = T, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
if (trans == "TRUE") fb = as.data.frame(t(fb)) #转置
# 判断指定的样本是否存在于矩阵数据2
if (length(group$sample[!(group$sample %in% rownames(fb))]) != 0)
{
	cat('[Error] 分组文件中，部分样本编号在文件中不存在! 请仔细检查分组文件 !\n')
	q()
}
fb = fb[group$sample, ]
fb_pca <- rda(fb, scale = TRUE)


proc <- procrustes(X = xq_pca, Y = fb_pca, symmetric = TRUE)
# summary(proc)
set.seed(123)
prot <- protest(X = xq_pca, Y = fb_pca, permutations = how(nperm = 999))

#提取 Procrustes 分析的坐标
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

Y$samples <- rownames(Y)

mycol = unlist(strsplit(color, split=","))
col = c()
for (i in 1:length(unique(group$groups))){
    col[group$groups == unique(group$groups)[i]] = mycol[i]
}

#ggplot2 作图
p <- ggplot(Y) +
geom_point(aes(X1, X2, color = "#26479b"), size = 1.5, shape = 16) +
geom_point(aes(PC1, PC2, color = "#ed9859"), size = 1.5, shape = 15) +
scale_color_manual(name = "group", values = c('#26479b' = '#26479b', '#ed9859' = '#ed9859'), breaks = c("#26479b", "#ed9859"), labels = c(mid2, mid1)) +
scale_shape_manual(name = "group", values = c('#26479b' = 16, '#ed9859' = 15)) +
geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.1, 'cm')), color = col, size = 0.3) +
theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "bottom", panel.background = element_rect(color = 'black', fill = 'transparent'), 
    legend.key = element_rect(fill = 'transparent')) +
labs(x = 'Dimension 1', y = 'Dimension 2', color = '') + 
geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
annotate('text', label = sprintf(paste0("M^2 == ", prot$ss)), 
    x = -0.21, y = 0.42, size = 3, parse = TRUE) +
annotate('text', label = sprintf(paste0("P < ", prot$signif)), 
    x = -0.21, y = 0.40, size = 3, parse = TRUE) + 
guides(color = guide_legend(override.aes = list(shape = c(16,15))))

pdf(file=paste(output, "procrustes.pdf", sep="/"))
p
dev.off()