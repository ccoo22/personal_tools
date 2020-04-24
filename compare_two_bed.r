#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: compare_two_bed.r  --bed1 <file> --bed2 <file> [--permutation <int> --genome <string> --output_file <file> --output_pdf <file> --pdf_width <int>]

Options:
    --bed1 <file>          第一个bed文件
    --bed2 <file>          第二个bed文件
    --output_file <file>   检验结果输出文件
    --output_pdf <file>    置换检验正态分布图输出
    --permutation <int>    置换检验频次 [default: 1000]
    --genome <string>      基因组符号,例如 hg19,mm9, dm2 [default: hg19]
    --pdf_width <int>      pdf宽度 [default: 14]
    --rlib <dir>           R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，比较两个bed文件，计算是否是随机的两个区域。p值越显著，关联性越大\n')
bed1            <- opts$bed1
bed2            <- opts$bed2
genome          <- opts$genome
permutation     <- as.integer(opts$permutation)
output_file     <- opts$output_file
output_pdf      <- opts$output_pdf
pdf_width       <- as.numeric(opts$pdf_width)
rlib            <- opts$rlib
# https://bioconductor.org/packages/release/bioc/vignettes/regioneR/inst/doc/regioneR.html
# 加载R包
message('加载regioneR')
.libPaths(rlib)
library(regioneR, quietly = TRUE)

# 读入数据
bed1_grange <- toGRanges(bed1, genome = genome)
bed2_grange <- toGRanges(bed2, genome = genome)


message("permutation test with times : ", permutation)
message("根据区域数量、迭代次数，运行时间可能会比较长")
pt <- overlapPermTest(A=bed1_grange, B=bed2_grange, ntimes=permutation)
result = summary(pt)

message("overlap check")

result$OverlapRegion  <- numOverlaps(bed1_grange, bed2_grange)  # The numOverlaps function receives 2 RS A and B and returns the number of regions in A that overlap a region in B.
result$RegionCountIn1 <- length(bed1_grange)
result$RegionCountIn2 <- length(bed2_grange)
result$MeanDistance   <- meanDistance(bed1_grange, bed2_grange) # given two RS A and B, computes the mean of the distance from every region in A to the closest region in B. It is useful to answer question of the type “Are my highly expressed genes closer to a certain TFBS than expected by chance?”


# 文件输出
if(!is.null(output_file)) write.table(result, output_file, quote = F, sep = "\t", row.names = FALSE, col.names = TRUE)
if(!is.null(output_pdf))
{
    pdf(output_pdf, width = pdf_width)
    plot(pt)
    dev.off()
}