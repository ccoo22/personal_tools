#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: maf_plot.r  -i <file> -o <pdf file> [ --show_gene_count <int> --rlib <dir> ]

Options:
    -i, --input <file>        输入maf文件
    -o, --output <pdf file>   输出pdf文件路径,示例：./a.pdf
    --show_gene_count <int>   展示前n个基因 [default: 50]
    --rlib <dir>              R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，肿瘤MAF绘图\n')
input             <- opts$input
output            <- opts$output
show_gene_count  <- opts$show_gene_count
no_summary_result <- opts$no_summary_result
rlib              <- opts$rlib

show_gene_count <- as.integer(show_gene_count)

if(!is.null(rlib)) .libPaths(rlib)

# input = '/home/ganb/work/Tumor/19B0708B/report/test.maf.txt'
message("加载maftools包，1-2分钟")
library(maftools) 
set.seed(91)

# 读入maf
maf = read.maf(maf = input)
# 文件中样本数量
sampleNum = length(levels(maf@data$Tumor_Sample_Barcode))

if(sampleNum == 0)
{
    message("MAF文件中，Tumor_Sample_Barcode列为空，不能绘图")
    q()
}



pdf(output, width = 12, height = 10 )

message("基因展示图")
oncoplot(maf = maf, top = show_gene_count)

message("基因展示图，显示样本名")
plot.new()
oncoplot(maf = maf, top = show_gene_count, showTumorSampleBarcodes = TRUE)

message("汇总统计图")
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()


# oncoplot(maf, top = 20, genes = NULL, mutsig = NULL,
#        mutsigQval = 0.1, drawRowBar = TRUE, drawColBar = TRUE,
#        clinicalFeatures = NULL, annotationDat = NULL,
#        annotationColor = NULL, genesToIgnore = NULL,
#        showTumorSampleBarcodes = FALSE, removeNonMutated = TRUE,
#        colors = NULL, sortByMutation = FALSE, sortByAnnotation = FALSE,
#        annotationOrder = NULL, keepGeneOrder = FALSE,
#        GeneOrderSort = TRUE, sampleOrder = NULL, writeMatrix = FALSE,
#        fontSize = 10, SampleNamefontSize = 10, titleFontSize = 15,
#        legendFontSize = 12, annotationFontSize = 12,
#        annotationTitleFontSize = 12, bgCol = "#CCCCCC", borderCol = NA,
#        colbar_pathway = FALSE)

# plotmafSummary(maf, rmOutlier = TRUE, dashboard = TRUE,
#   titvRaw = TRUE, addStat = NULL, showBarcodes = FALSE, fs = 1,
#   textSize = 0.8, color = NULL, titleSize = c(1, 0.8),
#   titvColor = NULL, top = 10)


