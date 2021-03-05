#!/home/genesky/software/r/4.0.2/bin/Rscript

library(docopt)

"Usage: plot_maf_oncoplot.r  -i <file> -o <pdf file> [--clinical <file> --sort_by_annotation --gene_mar <int> --barcode_mar <int> --show_tumor_sample --genes <string>  --draw_titv --gene_count <int> --pdf_width <int> --pdf_height <int>   --rlib <dir> ]

Options:
    -i, --input <file>        输入maf文件
    -o, --output <pdf file>   输出pdf文件路径,示例：./a.pdf
    -c, --clinical <file>     带有表头的样本临床信息文件。至少包含两列数据，第一列 Tumor_Sample_Barcode , 第二列及以后的列是样本分类。
                              第一列的名称必须是 Tumor_Sample_Barcode，
                              临床文件中的样本数量可以比maf中的多
                              表头名称仅允许字母、数字、下划线
    --sort_by_annotation      根据第一个临床表型，对样本进行排序
    --genes <string>          绘制指定基因, 基因名用逗号分隔
    --show_tumor_sample       显示肿瘤样本名称
    --draw_titv               在oncoplot图中添加titv图
    --gene_count <int>        展示前n个基因 [default: 50]
    --pdf_width <int>         pdf宽度 [default: 12]
    --pdf_height <int>        pdf高度 [default: 10]
    --gene_mar <int>          基因名区域的宽度(部分基因长度太大，显示不全，可以调高这个值) [default: 5]
    --barcode_mar <int>       样本名区域的宽度(部分样本名长度太大，显示不全，可以调高这个值) [default: 5]
    --rlib <dir>              R包路径 [default: /home/genesky/software/r/4.0.2/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，肿瘤MAF绘图\n')
input             <- opts$input
output            <- opts$output
clinical          <- opts$clinical
sort_by_annotation <- opts$sort_by_annotation
draw_titv         <- opts$draw_titv
show_tumor_sample <- opts$show_tumor_sample
genes             <- opts$genes
gene_count        <- as.numeric(opts$gene_count)
pdf_width         <- as.numeric(opts$pdf_width)
pdf_height        <- as.numeric(opts$pdf_height)
gene_mar          <- as.numeric(opts$gene_mar)
barcode_mar       <- as.numeric(opts$barcode_mar)
rlib              <- opts$rlib

if(!is.null(rlib)) .libPaths(rlib)
if(!is.null(genes)) genes = unlist(strsplit(genes, ','))

# input = '/home/ganb/work/Tumor/19B0708B/report/test.maf.txt'
message("加载maftools包，1-2分钟")
library(maftools) 
set.seed(91)

# 读入maf
clinical_name = NULL
if(is.null(clinical))
{
    maf = read.maf(maf = input)
}else{
    maf = read.maf(maf = input, clinicalData = clinical)
    clinical_name = colnames(maf@clinical.data)[2:ncol(maf@clinical.data)]
}

# 文件中样本数量
sampleNum = length(levels(maf@data$Tumor_Sample_Barcode))

if(sampleNum == 0)
{
    message("MAF文件中，Tumor_Sample_Barcode列为空，不能绘图")
    q()
}

# # 对突变类型的颜色进行指定
# vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
# names(vc_cols) = c(
#   'Frame_Shift_Del',
#   'Missense_Mutation',
#   'Nonsense_Mutation',
#   'Multi_Hit',
#   'Frame_Shift_Ins',
#   'In_Frame_Ins',
#   'Splice_Site',
#   'In_Frame_Del'
# )
# colors = vc_cols

pdf(output, width = pdf_width, height = pdf_height )
message("基因展示图")
oncoplot(maf = maf, top = gene_count, draw_titv = draw_titv, genes = genes, showTumorSampleBarcodes =  show_tumor_sample, gene_mar = gene_mar, barcode_mar = barcode_mar, clinicalFeatures = clinical_name, sortByAnnotation = sort_by_annotation)
dev.off()

# message("基因展示图")
# oncoplot(maf = maf, top = show_gene_count)

# message("基因展示图，显示样本名")
# plot.new()
# oncoplot(maf = maf, top = show_gene_count, showTumorSampleBarcodes = TRUE)

# message("汇总统计图")
# plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
# dev.off()


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


