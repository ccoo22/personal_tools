#/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)
options(scipen = 200)
"Usage: manhattan.r -i <file> -o <file> [--size <width,height> -s <suggestiveline> -g <genomewideline>]

Options:
  -i, --input <file>          the input file name, always filter.correlation.result.txt
  -o, --output <file>         the output directory name
  --size <width,height>       the size of plot [default: 900,500], 逗号分隔
  -s <suggestiveline>         the suggestiveline [default: 0.05]
  -g <genomewideline>         the genomewideline [default: 0.01] "-> doc

opts           <- docopt(doc)
input          <- opts$input
output         <- opts$output
size           <- opts$size
sline          <- opts$s
gline          <- opts$g

library(ggplot2)
library(qqman)

data <- read.table(input, sep='\t', header = T, check.names=FALSE, stringsAsFactors=FALSE)
size <- unlist(strsplit(size, ','))
#png(output, width = as.numeric(size[1]), height = as.numeric(size[2]))

#散点图
scatter_diagram_png <- paste(output, "/scatter_diagram.png", sep="")
png(scatter_diagram_png)
ggplot(data, aes(x = data$METH_RNA_dist, y = -log10(data$Pvalue))) + 
		geom_point() + 
		labs(x = "DISTANCE", y = bquote(paste(-log[10],"(p value)",sep="")))
dev.off()

# 曼哈顿图
mycol = c(119,132,147,454,89,404,123,62,463,461,128,139,552,28,54,84,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,529)
mycol = colors()[rep(mycol,20)]

data$plot_id <- paste(data$RNA_ID, data$METH_ID, sep="-")
data_plot <- data[, c('plot_id', 'METH_chr', 'METH_Position', 'Pvalue')]
colnames(data_plot) <- c('SNP', 'CHR', 'BP', 'P')
data_plot$CHR[data_plot$CHR == 'X'] = '23'
data_plot$CHR[data_plot$CHR == 'Y'] = '24'
data_plot$CHR <- as.numeric(data_plot$CHR)
manhattan_png <- paste(output, "/manhattan.png", sep="")
png(manhattan_png, width = as.numeric(size[1]), height = as.numeric(size[2]))
manhattan(data_plot, chrlabs=sort(unique(data_plot$CHR)), main="Manhattan", col=mycol) 
dev.off()

# QQ图
qq_png <- paste(output, "/QQ.png", sep="")
png(qq_png, width = 960, height = 960)
par(mar = c(7,7,9,5), xaxs = "i", yaxs = "i")
qq(data_plot$P, main = "QQ_plot", xlim = c(0, 6), ylim = c(0, 9), cex = 1.5, cex.lab = 1.8, cex.axis = 1.8, cex.main = 2.5, col = "blue4" , las = 1) 
dev.off()
