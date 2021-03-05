#!/usr/bin/env Rscript

library(docopt)
"Usage: differ_analysis_normal.r  -i <file> --sample_group_file <file> -o <file> -k <column num>
Options:
    -i , --input <file>                   the input file, each row is a feature, each column is sample
    -o , --output <file>                  the output file
    -k , --keep_column <column num>       keep input column to output, example: -k 1,2
    --sample_group_file <file>            the input group file, without header and only two column[1: sample, 2: group]" -> doc

opts                      <- docopt(doc, version = 'Program : metabolite diff analysis anova v1.0 \n          朱喜 320\n')
inputfile                 <- opts$input
output                    <- opts$output
groupfile                 <- opts$sample_group_file
keep_column_list          <- opts$k

# 测试用参数
# inputfile                 <- '/home/zhuxi/reserch/metabolite/MetaboAnalystR/script/normalized_data.txt'
# groupfile                 <- '/home/zhuxi/reserch/metabolite/MetaboAnalystR/script/group_anova.txt'
# output                    <- '/home/zhuxi/reserch/metabolite/MetaboAnalystR/script/tmp2/test.txt'
# keep_column_list          <- '1,2'


options(warn = -1)
library(tidyr, quietly = TRUE)   # gather 函数
library(fBasics, quietly = TRUE) # 正态分布检测

message("start metabolite_anova_analysis.r")

# 正态分布检测
run_normal_test <- function(data)
{
    # Shapiro-Wilk's正态检验,数据必须在3-5000,并且数据不能全相等
    if(length(data) < 3 || length(data) > 5000 || diff(range(data)) == 0)
    {
        pvalue <- NA
        wvalue <- NA
    }else
    {
        result <- normalTest(data, 'sw')
        pvalue <- result@test$p.value[[1]]
        wvalue <- result@test$statistic[[1]]
    } 
    return(c(wvalue, pvalue))
}

# 读入数据,data.frame格式
data          <- read.table(inputfile, header = T, sep = "\t", comment.char = "", check.names = F, stringsAsFactors = FALSE, quote = "")
sample.groups <- read.table(groupfile, header = F, sep = "\t", quote = "", stringsAsFactors = FALSE, colClasses = 'character')

sample_need   <- as.character(sample.groups[, 1])
sample_group  <- as.character(sample.groups[, 2])

group         <- unique(sample_group)
group_count   <- length(group)        # 分组数量
combine_count <- group_count * (group_count - 1) / 2  # 分组组合数量

###分组文件检测
if (length(sample_need[!(sample_need %in% colnames(data))]) != 0)
{
    cat('[Error] 分组文件中，部分样本编号在丰度文件中不存在! 请仔细检查分组文件！\n')
    q()
}
if(group_count < 3)
{
    cat('[Error] group < 3, please check groupfile !\n')
    q()
}
if('' %in% group || NA %in% group)
{
    cat('[Error] exist missing value, please check groupfile !\n')
    q()
}

# 获取需要保留的列编号
keep_column   = as.numeric(unlist(strsplit(keep_column_list, split = ",")))

# 数据整理，按照输入的分组信息提取数据，并整理
data_re_order <- matrix(nrow = nrow(data), ncol = 0) # 新数据
sample_group  <- c() # 对应新数据的样本分组
sample_name   <- c()
group_levels  <- c() # 组名称水平

for(count in 1:group_count)
{
    group_name        <- group[count]  # 组名
    group_sample_list <- as.character(sample.groups[, 1][sample.groups[, 2] == group_name])       # 组内样本

    group_data        <- data[, group_sample_list, drop = FALSE]
    group_data        <- as.matrix(group_data)
    # 如果数据只有一行，会返回的是向量，而不再是矩阵，所以不能用下面的方式
    # group_data        <- sapply(unlist(strsplit(group_sample_list, split=",")), function(x){data[[x]] })
    data_re_order     <- cbind(data_re_order, group_data)
    sample_group      <- c(sample_group, rep(group_name, ncol(group_data)))
    group_levels      <- c(group_levels, group_name)
    sample_name       <- c(sample_name, group_sample_list)
}

# 分组定义为因子类型，且levels中，control在前，case在后
sample_group <- factor(sample_group, levels = rev(group_levels))

# anova给出的分组组合方式， 后面用于表头
combine_pair <- t(combn(group_levels, 2))
combine_pair <- paste(combine_pair[ ,1], combine_pair[ ,2], sep = '-')

# 结果列表，空矩阵
result <- matrix(, nrow = nrow(data), ncol = length(keep_column) + 2 + length(sample_name)) # 保留表头、P值/FDR

# # 设置表头
# col_names <- c()

# # （1）描述信息
# col_names <- c(col_names, colnames(data)[keep_column])

# # （2）anova检测
# col_names <- c(col_names, 'P_value', 'FDR P_value')

# # (4) 分组配对数据统计
# pair_result_names <- unlist(lapply(combine_pair, function(x){ paste(x, c('diff', 'lwr', 'upr', 'p.adj'), sep = '_')}  ))
# col_names         <- c(col_names, pair_result_names)
colnames(result)  <- c(colnames(data)[keep_column], 'P_value', 'FDR P_value', sample_name)

################################################## 计算主体 ##########################################

for(row in 1:nrow(data))
{
    x_tmp         <- sample_group
    y_tmp         <- data_re_order[row, ]
    analysis_data <- data.frame(x = x_tmp, y = y_tmp)
    analysis_data <- analysis_data[complete.cases(analysis_data), ]
    x             <- analysis_data$x
    y             <- analysis_data$y

    # 去掉缺失值后，只有1组数据了，不能计算
    if(length(unique(x)) <= 1)
    {
        pvalue = NA
    }else
    {
        fit              <- aov(y ~ x)
        result_all       <- summary(fit)  # 提取p
        # result_each_pair <- TukeyHSD(fit) # 提取组组之间的分析结果
        pvalue           <- result_all[[1]]['x', 'Pr(>F)']

        # result_each_pair = as.data.frame(result_each_pair$x) # 转化为数据框
        # result_each_pair = data.frame(pair = row.names(result_each_pair), result_each_pair)
        # result_each_pair = gather(result_each_pair, type, statistic, -pair) # 矩阵转换形式，以TukeyHSD原始数据行、列名称重新给出数据

        # 按照列名填充数据，原因：去掉缺失样本后，分组的数量可能会减少。故通过这种方式更具有通用性。
        # result[row, paste(result_each_pair[, 1], result_each_pair[, 2], sep = '_')] = result_each_pair[, 3]
    }

    # 数据量统计
    annotation       <- as.matrix(data[row, keep_column])[1, ]

    result[row, ] = c(annotation, pvalue, NA, y_tmp)

}
# FDR矫正
result[, 'FDR P_value']       <- p.adjust(result[, 'P_value'], method = 'fdr')
write.table(result, output , sep = "\t", quote = F, row.names = F, col.names = TRUE )

# 差异ID
diff_id = result[,1][result[, 'P_value'] < 0.05]
write.table(diff_id, paste0(output, ".diff_id.txt"), sep = "\t", quote = F, row.names = F, col.names = F)


message("finish metabolite_anova_analysis.r")
