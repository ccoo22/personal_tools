#!/home/genesky/software/r/4.0.2/bin/Rscript

library(docopt)

"Usage: plot_circlize_label.r -i <file> --species <string> -o <pdf>  [--plot_type <string> --pdf_width <numeric> --pdf_heigth <numeric> --cex <numeric>  --rlib <dir> ]

Options:
    -i, --input <file>              输入文件，必须含有表头 chr start end text。 染色体要添加chr前缀
    --species <string>              物种名称，暂时支持： hg19/hg38/mm10 
    -o, --output <pdf>              输出pdf文件
    --plot_type <string>            ideogram显示的内容组合 ideogram,axis,labels [default: ideogram,axis,labels]
    --cex <numeric>                文字大小 [default: 0.8]
    --pdf_width <numeric>           PDF宽度 [default: 7]
    --pdf_height <numeric>          PDF高度 [default: 7]
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/4.0.2/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，圈图添加文字\n')
input               <- opts$input
species             <- opts$species
output              <- opts$output
plot_type           <- opts$plot_type
cex                 <- as.numeric(opts$cex)
pdf_width           <- as.numeric(opts$pdf_width)
pdf_height          <- as.numeric(opts$pdf_height)

rlib                <- opts$rlib
.libPaths(rlib)

# 导入 -> package
library("circlize")
set.seed(91)

###################################################################### 主程序
message("read input")
data_input      = read.table(input, head = T , check.names = F, sep = "\t", stringsAsFactors=F)

pdf(output, width=pdf_width, height=pdf_height)
# 90度开始
circos.par(start.degree=90, cell.padding = c(0, 0, 0, 0))

# ideogram 绘制
message("initializeWithIdeogram")
plot_type = unlist(strsplit(plot_type, ','))
circos.initializeWithIdeogram(species = species, plotType = plot_type)

# 添加Labels
message("genomicLabels")
circos.genomicLabels(data_input, labels.column = 4, side = "inside", cex = cex )

message("clear")
circos.clear()

dev.off()

