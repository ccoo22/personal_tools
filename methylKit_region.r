#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: methylKit_region.r  -o <dir> --case_group_name <string> --control_group_name <string> --case_sample <string> --control_sample <string> --case_file <string> --control_file <string> REGIONS...[--difference <int> --qvalue <float> --min_cov <int> --min_per_group <int> --thread <int>  --Rlib <dir>] 
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
   --Rlib <dir>                              methylKit R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]
Arguments:
	REGIONS                                  区域文件列表,格式：区域1名称 区域1文件 区域2名称 区域2文件...。注意:文件格式固定，表头至少要包含：#ID，Chr，Start，End" -> doc

opts                   <- docopt(doc, version = '基于methylKit的RRBS数据差异分析,区域水平 \n')
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
REGIONS                <- opts$REGIONS
 

# output_dir         <- './output'
# case_group_name    <- 'G1'
# control_group_name <- 'G2'
# case_sample        <- 'M-P10-1,M-P11-1,M-P12-1,P3-1,P9-1'
# control_sample     <- 'M-P10-2,M-P11-2,M-P12-2,P3-2,P9-2'
# case_file          <- '/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/M-P10-1_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/M-P11-1_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/M-P12-1_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/P3-1_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/P9-1_pe.bismark.cov.chr.gz'
# control_file       <- '/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/M-P10-2_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/M-P11-2_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/M-P12-2_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/P3-2_pe.bismark.cov.chr.gz,/home/wuj/project/Methylation_RRBS/16B0824B/report/bismark/P9-2_pe.bismark.cov.chr.gz'
# min_cov            <- 10
# min_per_group      <- 1
# thread             <- 10
# qvalue             <- 0.01
# difference         <- 25
# refgene            <- '/home/genesky/database/ucsc/hg19/gene/hg19_refGene.bed12'
# REGIONS <- c('Gene', '/home/genesky/database/ucsc/hg19/rrbs/hg19_RRBS_GENE.txt',
#              'Promoter', '/home/genesky/database/ucsc/hg19/rrbs/hg19_RRBS_Promoter500.txt',
#              'CPGI', '/home/genesky/database/ucsc/hg19/rrbs/hg19_RRBS_CpgIslandExt.txt',
#              'Target_Region', '/home/genesky/database/ucsc/hg19/rrbs/hg19_RRBS_Target_Region.txt'
#              )

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


################
# (2) 每个区域处理
################
region_count <- length(REGIONS) / 2
for( count in 1:region_count)
{	
	# 区域数据准备
	region_name <- REGIONS[2*count - 1]
	region_file <- REGIONS[2*count]
	cat('Process Region:', region_name, '\n')
	cat('    Generate Region Meth Data', '\n')
	target_region <- read.table(region_file, header = TRUE, comment.char="", check.names = F) 
    target_GRanges <- GRanges(seqnames = target_region$Chr, ranges = IRanges(start = target_region$Start, end = target_region$End))  # 构建区域信息
	target_obj  <- regionCounts(myobj, target_GRanges) # 区域甲基化汇总
	target_meth <- unite(target_obj,  destrand=FALSE, min.per.group = as.integer(min_per_group))

	# 差异分析
	cat('    Start Region Diff Analysis', '\n')
	target_diff = calculateDiffMeth(target_meth, mc.cores = as.integer(thread))
	myDiff25p <- getMethylDiff(target_diff, difference=as.integer(difference), qvalue=as.numeric(qvalue) ) # 按照要求提取差异结果

	# 数据输出准备
	cat('    Prepare output Diff Result', '\n')
	all_meth           <- getData(target_meth)        # 总甲基化信息      
	all_meth_diff      <- getData(target_diff)      # 总甲基化差异信息
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
	all_meth$pvalue     <- target_diff$pvalue
	all_meth$qvalue     <- target_diff$qvalue
	all_meth$meth.diff  <- target_diff$meth.diff

	# 获取显著差异的详细信息
	is_sig_meth_diff = paste(all_meth$chr, all_meth$start, all_meth$end, sep = '_') %in% paste(sig_meth_diff$chr, sig_meth_diff$start, sig_meth_diff$end, sep = '_')
	sig_meth = all_meth[is_sig_meth_diff, ]


	write.table(all_meth, file = paste0(output_dir, '/region_', region_name,'_meth_all_result.txt'),  append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
	write.table(sig_meth, file = paste0(output_dir, '/region_', region_name,'_meth_sig_result.txt'),  append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))

	# 差异统计
	sig_count = nrow(sig_meth)
	all_count = nrow(all_meth)
	sig_perc  = sig_count / all_count

	stats = matrix(c(all_count, sig_count, sig_perc), nrow = 1, ncol = 3) 
	colnames(stats) = c('ALL_METH_COUNT', 'SIG_DIFFERENT_METH_COUNT', 'SIG_PERCENT')
	write.table(stats, file = paste0(output_dir, '/region_', region_name,'_meth_all_stats.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}