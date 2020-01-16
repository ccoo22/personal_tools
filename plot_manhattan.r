#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: manhattan.r -i <file> -o <file> [--size <width,height> -s <suggestiveline> -g <genomewideline>]

Options:
  -i, --input <file>          the input file name, 输入文件包括表头，1到4列依次为染色体号、SNP、位置、P值
  -o, --output <file>         the output file name
  --size <width,height>       the size of plot [default: 900,500], 逗号分隔
  -s <suggestiveline>         the suggestiveline [default: 0.05]
  -g <genomewideline>         the genomewideline [default: 0.01] "-> doc

opts           <- docopt(doc)
input          <- opts$input
output         <- opts$output
size           <- opts$size
sline          <- opts$s
gline          <- opts$g

library(qqman)

data <- read.table(input, header = T)
colnames(data) <- c("CHR","SNP","BP","P")
data <- na.omit(data)
mycol = c(119,132,147,454,89,404,123,62,463,461,128,139,552,28,54,84,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,529)
mycol = colors()[rep(mycol,20)]
size <- unlist(strsplit(size, ','))
png(output, width = as.numeric(size[1]), height = as.numeric(size[2]))
manhattan(data, col=mycol, suggestiveline = -log10(as.numeric(sline)), genomewideline = -log10(as.numeric(gline)))

dev.off()


