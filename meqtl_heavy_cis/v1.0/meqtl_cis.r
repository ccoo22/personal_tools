#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: meqtl_cis.r --output <file> --exp_file <file> --exp_loc_file <file> --snv_file <file> --snv_loc_file <file>  [--cisdist <string> --cis_p <string> --cov_file <file> --rlib <dir> --thread <int> --exp_anno <file> --snp_anno <file>]
Options:
    --output <file>              CIS 结果输出文件
    --exp_file <file>            表达量数据矩阵，第一列是id，后面n列都是样本，数值型数据，如果缺失，空着。含有表头。务必保证样本名和顺序与snp_file中的一致。
    --exp_loc_file <file>        表达量矩阵id的坐标信息，四列数据，含有表头。分别是： id chr start end
    --snv_file <file>            突变量化矩阵，第一列是id, 后面n列都是样本，数值型数据，如果缺失，空着。含有表头。通常是 0、1、2
    --snv_loc_file <file>        突变量化矩阵id的坐标信息，三列数据，含有表头，分别是：id chr pos
    --cisdist <string>           CIS分析的距离 [default: 1e6]
    --cis_p <string>             CIS分析的P值限制，大于该P值，结果不输出 [default: 0.05]
    --cov_file <file>            协变量矩阵，可以不输入 
    --exp_anno <file>            基因注释文件，行数要与表达量矩阵相等，且ID一致，第一列为gene ID, 后面所有列都是注释，把该文件的所有列信息加入结果文件中。
    --snp_anno <file>            snp注释文件，行数要与突变矩阵相等，且ID一致，第一列为snp ID，把该文件的所有列信息加入结果文件中。
    --thread <int>               并行线程数量，加快后期数据统计 [default: 10]
    --rlib <dir>                 R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts            <- docopt(doc, version = 'Program : EQTL cis analysis based on mRNA and SNV data \n        甘斌 129\n')
exp_file        <- opts$exp_file
exp_loc_file    <- opts$exp_loc_file
snv_file        <- opts$snv_file
snv_loc_file    <- opts$snv_loc_file
output          <- opts$output
cov_file        <- opts$cov_file
cis_pval        <- opts$cis_p
cisDist_val     <- as.numeric(opts$cisdist)
thread          <- as.integer(opts$thread)
exp_anno        <- opts$exp_anno
snp_anno        <- opts$snp_anno
rlib            <- opts$rlib

if(!is.null(rlib)) .libPaths(rlib)

# 测试参数
# exp_file        <- "exp.raw.txt"
# exp_loc_file    <- "exp.pos.txt"
# snv_file        <- "snp.code.txt"
# snv_loc_file    <- "snp.pos.txt"
# output          <- "b.txt"
# cov_file        <- NULL
# cis_pval        <- 0.05
# cisDist_val     <- 1e6
 

#########
library(MatrixEQTL)
library(parallel)
library(CMplot)
library(ggpubr)

useModel = modelLINEAR

# Genotype file name
SNP_file_name <- snv_file
snps_location_file_name <- snv_loc_file

# Gene expression file name
expression_file_name <- exp_file
gene_location_file_name <- exp_loc_file

# Covariates file name
covariates_file_name <- cov_file

# Output file name
output_file_name_cis <- output

# Only associations significant at this level will be saved
# pvOutputThreshold_cis = 0.05
# pvOutputThreshold_tra = 0

pvOutputThreshold_cis = as.numeric(cis_pval)

# Set to numeric() for identity
errorCovariance = numeric()

# cisDist = 1e6
cisDist = as.numeric(cisDist_val)

## Load genotype data
message("Load genotype data")
snps = SlicedData$new()
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values;
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)

## Load gene expression data
message("Load gene expression data")
gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

## Load covariates
cvrt = SlicedData$new()
if(!is.null(cov_file)) {
    cvrt$fileDelimiter = "\t"      # the TAB character
    cvrt$fileOmitCharacters = "NA" # denote missing values;
    cvrt$fileSkipRows = 1          # one row of column labels
    cvrt$fileSkipColumns = 1       # one column of row labels
    if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name)
    }
}

## Run the analysis 临时不输出
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)
row.names(snpspos) = snpspos[, 1]
row.names(genepos) = genepos[, 1]

me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name      = NULL,
    pvOutputThreshold     = 0,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis  = NULL,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

################
# 补充数据
################
result = me$cis$eqtls
pair_count = nrow(result)
message("对 ", pair_count, " 个pair结果添加 chr/pos/distance等信息。会比较慢，多线程并行数量: ", thread)

# 转换字符型
result$snps = as.character(result$snps)
result$gene = as.character(result$gene)

###########
# (1) 表达量、样本数量统计
###########
gene_exp    = read.table(expression_file_name, header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t", check.names = F)
snp_val     = read.table(SNP_file_name, header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t", check.names = F)

# 并行线程，统计每个配对情况下的数据统计
cl <- makeCluster(thread)
clusterExport(cl, c("gene_exp", "snp_val", "result"))  # 通知并行线程的共享变量
tongji <-  parLapply(cl, 1:nrow(result), function(row){
        gene = as.numeric(gene_exp[result[row, 'gene'], ])
        snp = as.numeric(snp_val[result[row, 'snps'], ])
        data <- data.frame(gene, snp)
        data <- data[complete.cases(data), ] # 去掉缺失值

        zero <- sum(data$snp == 0)
        one <- sum(data$snp == 1)
        two <- sum(data$snp == 2)
        data_mean = mean(data$gene)
        data_median = median(data$gene)
        c(nrow(data), paste(zero, one, two, sep="|"), data_mean, data_median)
      })
stopCluster(cl)
tongji <- data.frame(matrix(unlist(tongji), nrow=nrow(result), byrow=T))

result$sampleNum     = tongji[,1]
result$snp_SampleNum = tongji[,2]

###########
# (2) 添加SNP注释
###########
result$snp_chr    = snpspos[result$snps, 2]
result$snp_pos    = snpspos[result$snps, 3]
if(!is.null(snp_anno))
{
    snp_anno_data    = read.table(snp_anno, header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t", check.names = F)
    snp_anno_colnames = colnames(snp_anno_data)
    start = ncol(result) +  1
    result <- cbind(result, snp_anno_data[result$snps, ])
    colnames(result)[start:ncol(result)] = snp_anno_colnames  # 添加列名
}

###########
# (3) 添加基因注释
###########
result$gene_chr       = genepos[result$gene, 2]
result$gene_left_pos  = genepos[result$gene, 3]
result$gene_right_pos = genepos[result$gene, 4]
result$gene_mean      = tongji[,3]  
result$gene_median    = tongji[,4]
if(!is.null(exp_anno))
{
    exp_anno_data    = read.table(exp_anno, header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t", check.names = F)
    exp_anno_colnames = colnames(exp_anno_data)
    start = ncol(result) +  1
    result <- cbind(result, exp_anno_data[result$gene, ])
    colnames(result)[start:ncol(result)] = exp_anno_colnames  # 添加列名
}
result$snp_gene_distance = result$snp_pos - result$gene_left_pos

write.table(result, output, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

###########
# (4) 绘图
###########
message("绘图")
output_dir = dirname(output)
setwd(output_dir)

# SNP manhanttan
snp_manhattan = result[, c('snps', 'snp_chr', 'snp_pos', 'pvalue')]
snp_manhattan = snp_manhattan[order(snp_manhattan$pvalue), ]  # 排序
snp_manhattan = snp_manhattan[!duplicated(snp_manhattan$snps), ]  # 去重复，保留Pvalue最小的结果
CMplot(snp_manhattan,plot.type="m",band=0,LOG10=TRUE)
file.rename("Rectangular-Manhattan.pvalue.jpg", "snp.most_sig_pvalue.manhattan.jpg")

# gene manhanttan
gene_manhattan = result[, c('gene', 'gene_chr', 'gene_left_pos', 'pvalue')]
gene_manhattan = gene_manhattan[order(gene_manhattan$pvalue), ]  # 排序
gene_manhattan = gene_manhattan[!duplicated(gene_manhattan$gene), ]  # 去重复，保留Pvalue最小的结果
CMplot(gene_manhattan,plot.type="m",band=0,LOG10=TRUE)
file.rename("Rectangular-Manhattan.pvalue.jpg", "gene.most_sig_pvalue.manhattan.jpg")

# 散点图
scatter_data = result[, c('snp_gene_distance', 'pvalue')]
scatter_data$pvalue = -log10(scatter_data$pvalue)

png("scatter.png",width = 1000, height = 1000) 
p <- ggscatter(scatter_data, x='snp_gene_distance', y='pvalue', color = 'black', 
          palette = "npg", #杂志nature的配色
          xlim = c(-cisDist_val, cisDist_val)
          )  + ylab('-log10(pvalue)') + xlab('Distance')  + theme(plot.title = element_text(hjust = 0.5))  # 去掉legend,标题居中  
print(p)
dev.off()
