#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: methylKit_site.r  -o <dir> --case_group_name <string> --control_group_name <string> --case_sample <string> --control_sample <string> --case_file <string> --control_file <string> [--difference <int> --qvalue <float> --min_cov <int> --min_per_group <int> --thread <int>  --Rlib <dir> --refgene <string> --cpgi <string>]
Options:
   -o, --output_dir <dir>                    结果输出目录，要提前建立
   --case_group_name <string>                case组名称
   --control_group_name <string>             control组名称
   --case_sample <string>                    case组样本列表，样本间用‘逗号’分隔
   --control_sample <string>                 control组样本列表，样本间用‘逗号’分隔
   --case_file <string>                      case组样本甲基化文件列表，样本间用‘逗号’分隔
   --control_file <string>                   control组样本甲基化文件列表，样本间用‘逗号’分隔 
   --min_cov <int>                           最低测序深度 [default: 10] 
   --min_per_group <int>                     差异分析时，每一组最低样本数量，样本数量过少，则跳过分析 [default: 1] 
   --qvalue <float>                          差异分析结果qvalue过滤阈值, <qvalue [default: 0.01] 
   --difference <int>                        差异分析结果组间平均甲基化程度差异阈值，百分制, |dif| > difference [default: 25] 
   --thread <int>                            差异分析时，并行线程数 [default: 10] 
   --refgene <string>                        refgene注释数据库，用于展示差异位点分布 [default: /home/genesky/database/ucsc/hg19/gene/hg19_refGene.bed12] 
   --cpgi <string>                           CPGI注释数据库，用于展示差异位点分布 [default: /home/genesky/database/ucsc/hg19/rrbs/cpgi.bed] 
   --Rlib <dir>                              methylKit R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts                   <- docopt(doc, version = '基于methylKit的RRBS数据差异分析,位点水平 \n')
output_dir             <- opts$output_dir
case_group_name        <- opts$case_group_name
control_group_name     <- opts$control_group_name
case_sample            <- opts$case_sample
control_sample         <- opts$control_sample
case_file              <- opts$case_file
control_file           <- opts$control_file
min_cov                <- opts$min_cov
min_per_group          <- opts$min_per_group
qvalue                 <- opts$qvalue
difference             <- opts$difference
thread                 <- opts$thread
Rlib                   <- opts$Rlib
refgene                <- opts$refgene
cpgi                   <- opts$cpgi

# output_dir         <- './output'
# case_group_name    <- 'G1'
# control_group_name <- 'G2'
# case_sample        <- 'M-P10-1,M-P11-1,M-P12-1,P3-1'
# control_sample     <- 'M-P10-2,M-P11-2,M-P12-2,P3-2'
# case_file          <- '/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/M-P10-1_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/M-P11-1_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/M-P12-1_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/P3-1_pe.bismark.cov.chr.gz'
# control_file       <- '/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/M-P10-2_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/M-P11-2_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/M-P12-2_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/P3-2_pe.bismark.cov.chr.gz'
# min_cov            <- 10
# min_per_group      <- 1
# thread             <- 10
# qvalue             <- 0.01
# difference         <- 25
# refgene            <- '/home/genesky/database/ucsc/hg19/gene/hg19_refGene.bed12'
# cpgi               <- '/home/genesky/database/ucsc/hg19/rrbs/cpgi.bed'
# Rlib               <- '/home/genesky/software/r/3.5.1/lib64/R/library'
# other_sample <- 'C1,C2,C3,C4'
# other_file <- '/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/C1_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/C2_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/C3_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/C4_pe.bismark.cov.chr.gz'

options(warn=-1)
.libPaths(Rlib)
library(methylKit)
library(genomation)

################
# (1) 数据读入
################
# hg19
cat('Read RRBS data', '\n')
min_cov             <- as.numeric(min_cov)
file_list           <- c(unlist(strsplit(case_file, ',')), unlist(strsplit(control_file, ',')))
case_sample_list    <- unlist(strsplit(case_sample, ','))
control_sample_list <- unlist(strsplit(control_sample, ','))
sample_list         <- c(case_sample_list, control_sample_list)
treatment           <- c(rep(1, length(case_sample_list)), rep(0, length(control_sample_list)))
sample_count        <- length(sample_list)

myobj <- methRead( location = as.list(file_list), 
				sample.id = as.list(sample_list), 
				assembly = "hg19",
				treatment = treatment,
				context = "CpG",
				pipeline = "bismarkCoverage", 
				mincov = min_cov
	)

# 数据合并
cat('Unite RRBS data', '\n')
meth <- unite(myobj, destrand=FALSE)

################
# (2) 数据初期质控
################

# 样本相关性图
cat('Quality Control : Sample Correlation', '\n')
pdf(paste0(output_dir, '/sample_correlation.pdf'), width = 2.5 * sample_count, height = 2.5 * sample_count)
pdf_tmp <-dev.cur()
png(paste0(output_dir, '/sample_correlation.png'), width = 150 * sample_count, height = 150 * sample_count)
dev.control("enable")
waste_message <- capture.output(getCorrelation(meth, plot=TRUE))
dev.copy(which=pdf_tmp)  #复制来自png设备的图片到pdf
dev.off()
dev.off()

# 样本聚类图
cat('Quality Control : Sample Cluster', '\n')
pdf(paste0(output_dir, '/sample_cluster.pdf'), width = 14, height = 8)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE) 
dev.off() 	

# 样本pca图
cat('Quality Control : Sample PCA', '\n')
pdf(paste0(output_dir, '/sample_pca.pdf'), width = 14, height = 10)
pca_index <- PCASamples(meth, obj.return = TRUE) 
write.csv(pca_index$x, file = paste0(output_dir, '/sample_pca.csv'))
dev.off() 

# 样本甲基化统计图
cat('Quality Control : Sample methylation stats', '\n')
pdf(paste0(output_dir, '/Sample_methylation_stats.pdf'))
for(i in 1:sample_count)
{
	getMethylationStats(myobj[[i]], plot = TRUE, both.strands = FALSE)
}
dev.off()

# 样本覆盖深度图
cat('Quality Control : Sample coverage stats', '\n')
pdf(paste0(output_dir, '/Sample_coverage_stats.pdf'))
for(i in 1:sample_count)
{
	getCoverageStats(myobj[[i]], plot = TRUE, both.strands = FALSE)
}
dev.off()

################
# (3) 差异分析
################
# 重新汇总，设定每组最小样本数量
cat('Diff Analysis', '\n')
meth <- unite(myobj, destrand=FALSE, min.per.group = as.integer(min_per_group) )
myDiff <- calculateDiffMeth(meth, mc.cores = as.integer(thread)) ##num.cores设定计算核心数加快计算速度 ###差异计算，根据样本大小，自动选择Fisher’s exact or logistic regression计算P值，并使用 SLIM方法校正得到Q-values
myDiff25p <- getMethylDiff(myDiff, difference=as.integer(difference), qvalue=as.numeric(qvalue) ) # 按照要求提取差异结果

cat('Prepare output Diff Result', '\n')
# 数据输出准备
all_meth           <- getData(meth)        # 总甲基化信息      
all_meth_diff      <- getData(myDiff)      # 总甲基化差异信息
sig_meth_diff      <- getData(myDiff25p)   # 显著差异甲基化信息

# 样本名替换
colnames(all_meth) = c('chr', 'start', 'end', 'strand', unlist(lapply(sample_list, function(x){ paste(x, c('coverage', 'numCs', 'numTs'), sep = '_')})) )
# 计算样本甲基化程度
for(sample in sample_list)
{	
	sample_freq_title = paste(sample, 'freq', sep = '_')
	sample_cov_title  = paste(sample, 'coverage', sep = '_')
	sample_c_title    = paste(sample, 'numCs', sep = '_')
	all_meth[[sample_freq_title]] = all_meth[[sample_c_title]] / all_meth[[sample_cov_title]]
}
# 提取差异信息
all_meth$pvalue     <- myDiff$pvalue
all_meth$qvalue     <- myDiff$qvalue
all_meth$meth.diff  <- myDiff$meth.diff

# 获取显著差异的详细信息
is_sig_meth_diff = paste(all_meth$chr, all_meth$start, sep = '_') %in% paste(sig_meth_diff$chr, sig_meth_diff$start, sep = '_')
sig_meth = all_meth[is_sig_meth_diff, ]


write.table(all_meth, file = paste0(output_dir, '/meth_all_result.txt'),  append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
write.table(sig_meth, file = paste0(output_dir, '/meth_sig_result.txt'),  append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
# 差异统计
sig_count = nrow(sig_meth)
all_count = nrow(all_meth)
sig_perc  = sig_count / all_count

stats = matrix(c(all_count, sig_count, sig_perc), nrow = 1, ncol = 3) 
colnames(stats) = c('ALL_METH_COUNT', 'SIG_DIFFERENT_METH_COUNT', 'SIG_PERCENT')
write.table(stats, file = paste0(output_dir, '/meth_all_stats.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# 输出样本分组信息
sample_info = data.frame(Group = c(case_group_name, control_group_name), Sample = c(case_sample, control_sample))
write.table(sample_info, file = paste0(output_dir, '/meth_all_group_sample.txt'),  append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))

################
# (4) 差异分布
################

#差异位点染色体上分布
cat('Diff Distribution in chr', '\n')
pdf(paste0(output_dir, '/diff_meth_chr_distribution.pdf'))
result = tryCatch( {diffMethPerChr(myDiff, plot=TRUE, qvalue.cutoff=as.numeric(qvalue), meth.cutoff=as.integer(difference))}, 
	                error = function(e){ cat("ERROR :",conditionMessage(e), "\n")} 
	                ) 
dev.off()

# 差异位点基因注释
if(file.exists(refgene) && nrow(myDiff25p) > 0)
{	
	cat('Diff Distribution in gene', '\n')
	pdf(paste0(output_dir, '/diff_meth_gene_distribution.pdf'))
	gene.obj <- readTranscriptFeatures(refgene)
	diff_gene_ano  <- annotateWithGeneParts(as(myDiff25p, "GRanges"), gene.obj)
	result = tryCatch( { plotTargetAnnotation(diff_gene_ano, precedence=TRUE, main="differential methylation annotation by genes") },
		                error = function(e){cat("ERROR :",conditionMessage(e),"\n")}
		)
	dev.off()
}

# 差异位点CPGI注释
if(file.exists(cpgi) && nrow(myDiff25p) > 0)
{	
	cat('Diff Distribution in CPGI', '\n')
	pdf(paste0(output_dir, '/diff_meth_CPGI_distribution.pdf'))
	cpg.obj = readFeatureFlank(cpgi, feature.flank.name = c("CpGi", "shores"))
	diff_cpgi_ano=annotateWithFeatureFlank(as(myDiff25p, "GRanges"), cpg.obj$CpGi, cpg.obj$shores, feature.name = "CpGi", flank.name = "shores")
	result = tryCatch({ plotTargetAnnotation(diff_cpgi_ano, col=c("green","gray","white"), main="differential methylation annotation CpG islands") },
		                error=function(e){cat("ERROR :",conditionMessage(e),"\n")}
		)
	dev.off()
}
  
 

