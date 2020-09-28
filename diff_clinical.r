#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: diff_clinical.r -i <file> -o <dir> --sample_group <file>  --case_group_name <string> --control_group_name <string> [--rlib <dir>]

Options:
    -i, --input <file>              临床文件矩阵，每一行对应一个样本，每一列对应一个特征，第一列是样本名，第一行是特征名称，第二列及之后的信息为所有要分析的临床信息，如果有缺失，空着即可。
                                    临床信息可以是连续型数值（ttest）、离散数值（卡方, 当uniq数值数量 <= 4 认为是离散，否则认为是连续）、离散型字符（卡方）
    --sample_group <file>           样本分组文件，两列数据，第一列样本名，第二列样本分组，有表头。
    --case_group_name <string>      分组文件中，case分组名称
    --control_group_name <string>   分组文件中，control分组名称
    -o, --output_file <dir>         输出文件
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，做临床样本做两组之间的差异分析，支持：连续变量（ttest）、离散变量(卡方)、字符变量（卡方）\n')
input               <- opts$input
sample_group        <- opts$sample_group
case_group_name     <- opts$case_group_name
control_group_name  <- opts$control_group_name
output_file         <- opts$output_file

rlib                <- opts$rlib
.libPaths(rlib)

###################################################################### 主程序
message("read input")
data_input = read.table(input, header = T, row.names = 1, check.names = F, sep = "\t", stringsAsFactors=F)

# 读取分组信息
message("read group")
sampleinfo <- read.table(sample_group, header = TRUE, sep = "\t", colClasses = 'character') # 读取分组信息
colnames(sampleinfo) <- c("sample", "group")  
rownames(sampleinfo) <- sampleinfo[, 1] 

# 错误检测
message("sample/group check")
if(!case_group_name %in% sampleinfo$group)
{
    message("[Error] 你输入的case分组名称不在 sample_group中")
    q()
}
if(!control_group_name %in% sampleinfo$group)
{
    message("[Error] 你输入的control分组名称不在 sample_group中")
    q()
}
if(sum(!sampleinfo$sample %in% rownames(data_input)) > 0)
{
    losts = sampleinfo$sample[!sampleinfo$sample %in% rownames(data_input)]
    message("[Error] sample_group中的样本在input文件中缺失:", losts)
    q()
}

# 分析准备
message("prepare analysis")
case_samples      = sampleinfo$sample[which(sampleinfo$group == case_group_name)]
control_samples   = sampleinfo$sample[which(sampleinfo$group == control_group_name)]

result = matrix(0,0,4)
colnames(result) = c('Variables', 'Case cohort No. (%)', 'Control cohort No. (%)', 'P-value')
result = rbind(result, c('No. of patients' , length(case_samples), length(control_samples), ''))

# 分析每一个临床特征
for(variable in colnames(data_input))
{   
    message('process variable: ', variable)
    case_data = data_input[case_samples, variable]
    control_data = data_input[control_samples, variable]

    # 去掉缺失、空值
    case_data = case_data[!is.na(case_data) & case_data!='']
    control_data = control_data[!is.na(control_data) & control_data!='']

    # uniq数量
    uniq_count = unique(c(case_data, control_data))
    # 数据类型
    data_class = class(data_input[, variable])


    # 字符变量、数值变量小于等于4个的，都用卡方
    if(data_class == "character" || uniq_count <= 4)
    {
        data_tmp = data.frame(value = c(case_data, control_data),
                               class=c(rep('case', length(case_data)), rep('control', length(control_data)))
                               )  # value class 顺序不要变
        data_table = table(data_tmp) # 数量矩阵， 每一行是一个特征，列分别是case/control
        data_prop  = round(prop.table(data_table, 2) * 100, 2) # 比例矩阵 按列计算，双精度 %
        pvalue = chisq.test(data_table, correct = FALSE)$p.value
        result = rbind(result, c(variable, '', '', pvalue))

        # 把每一个值依次并入结果
        for(row in 1:nrow(data_table))
        {
            item  = rownames(data_table)[row]
            count = data_table[item, ]  # 数量： 615     286
            prob  = data_prop[item, ]  # 比例 ： 97.00   98.96
            count_prob = paste0(paste(count, prob, sep='('), ')')  # 联合 ： "615(97)"    "286(98.96)"
            # 结果并入
            result = rbind(result, c(item, count_prob, ''))
        }
    }else  # ttest
    {
        pvalue = t.test(case_data, control_data)$p.value
        case_mean    = round(mean(case_data), 2)
        case_sd      = round(sd(case_data), 2)
        control_mean = round(mean(control_data), 2)
        control_sd   = round(sd(control_data), 2)
        result = rbind(result, c(variable, paste(case_mean, case_sd, sep = ' \u00b1 ') , paste(control_mean, control_sd, sep = ' \u00b1 '), pvalue ))
    }
}


write.table(result, output_file, quote = FALSE, row.names = FALSE, sep = '\t')


