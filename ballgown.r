#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: ballgown.r -i <dir> -o <dir> --case_sample_list <string> [--control_sample_list <string> --rlib <dir>]

Options:
    -i, --input_dir <dir>           ballgown格式输入路径
    -o, --output_dir <dir>          结果输出目录
    --case_sample_list <string>     case组样本列表，用“逗号”分隔,例如：C1,C2,C3,C4
    --control_sample_list <string>  control组样本列表，用“逗号”分隔，例如：S1,S2,S3,S5
                                    注意：如果没有提供control_sample_list参数，则不进行差异分析，只给出FPKM信息
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，ballgown 定量FPKM与差异分析\n')
input_dir              <- opts$input_dir
output_dir             <- opts$output_dir
case_sample_list       <- opts$case_sample_list
control_sample_list    <- opts$control_sample_list
rlib                   <- opts$rlib
#case_sample_list    = 'T26_WTS,T10_WTS,T12_WTS,T18_WTS,T19_WTS,T23_WTS'
#control_sample_list = 'N19_WTS,N1_WTS,N10_WTS,N12_WTS,N18_WTS'

.libPaths(rlib)

# 加载ballgown
library(ballgown)

# 确认样本
case_samples        = unlist(strsplit(case_sample_list, ','))  # case分组样本 
control_samples     = c()  # control分组样本 
file_dir = paste(input_dir, case_samples, sep = '/')  # 样本数据路径
groups   = rep(1, length(case_samples)) # 分组信息

# 是否存在control样本
if(! is.null(control_sample_list))
{   
    control_samples = unlist(strsplit(control_sample_list, ','))
    file_dir = c(file_dir, paste(input_dir, control_samples, sep = '/'))
    groups   = c(groups, rep(0, length(control_samples)))
}

####################################################################################
message("读入ballgown数据")
bg = ballgown(samples=file_dir, meas = 'FPKM')
pData(bg) = data.frame(id=sampleNames(bg), group=groups)  # 添加分组信息

message("提取fpkm数据")
transcript_fpkm = texpr(bg, 'all')  # 提取FPKM数据，返回data.frame格式，其中，前10列包含了转录本的多个注释信息：t_id,chr,strand,start,end,t_name,num_exons,length,gene_id,gene_name
colnames(transcript_fpkm)   = gsub("FPKM.", '', colnames(transcript_fpkm)) # 上面的结果中，会额外添加FPKM前缀，这里去除
transcript_fpkm = transcript_fpkm[, -1] # 删除第一列id,没用



# 存在control样本，计算差异
if(! is.null(control_sample_list))
{   
    message("差异分析")
    stat_results                = stattest(bg, feature='transcript', meas='FPKM', covariate='group')  # 结果行数与transcript_fpkm的行数、id一一对应
    transcript_fpkm$pval        = stat_results$pval
    transcript_fpkm$qval        = stat_results$qval
    transcript_fpkm$caseMean    = apply(transcript_fpkm[, case_samples], 1, mean) # 计算case组均值
    transcript_fpkm$controlMean = apply(transcript_fpkm[, control_samples], 1, mean) # 计算control组均值
    transcript_fpkm$log2FC      = log2(transcript_fpkm$caseMean / transcript_fpkm$controlMean) # 计算log2(FoldChange)
}


output_file <- ifelse(is.null(control_sample_list), paste0(output_dir, "/transcript_fpkm.txt", ""), paste0(output_dir, "/transcript_fpkm.DE.txt", ""))
message("结果输出：", output_file)
write.table(transcript_fpkm, output_file, quote = F, sep = "\t", row.names = FALSE)
 
