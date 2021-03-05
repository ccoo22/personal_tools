#!/usr/bin/env Rscript

library(docopt)

"Usage: diff.R  -i <file>  --sample_group_file <file>  -o <file> --volcano <pdf>   -k <column num>  [-m <ttest model> --LOG2FC <num>]
Options:
   -i , --input <file>                   the input file, with header, each row is a point, each column is sample
   -o , --output <file>                  the output file
   --volcano <pdf>                       volcano plot pdf file
   -m , --model <ttest model>            independent and dependent [default: independent]
   -k , --keep_column <column num>       keep input column to output, example: -k 1,2
   --sample_group_file <file>            the input group file, without header and only two column(1: sample, 2: group)
   --LOG2FC <num>                        the log2 fold change [default: 1]" -> doc

opts                      <- docopt(doc, version = 'Program : metabolite diff analysis v1.0 \n          朱喜 320\n')
inputfile                 <- opts$input
groupfile                 <- opts$sample_group_file
output                    <- opts$output
volcano                   <- opts$volcano
model                     <- opts$m
keep_column_list          <- opts$k
LOG2FC                    <- as.numeric(opts$LOG2FC)

# 测试用参数
# inputfile                 <- '/home/zhuxi/reserch/metabolite/MetaboAnalystR/script/normalized_data.txt'
# groupfile                 <- '/home/zhuxi/reserch/metabolite/MetaboAnalystR/script/groupfile.txt'
# output                    <- '/home/zhuxi/reserch/metabolite/MetaboAnalystR/script/tmp'
# model                     <- 'independent'
# keep_column_list          <- '1,2'
# LOG2FC                    <- as.numeric('1')

message("start metabolite_diff_analysis.r")
 
options(warn=-1)
library(ggplot2)
library(RColorBrewer) # 颜色包
library(gplots)

# TTEST子函数
run_ttest <- function(data1, data2, model)
{
    # t检验， 返回 Pvalue Tvalue
    is_good_run = 1
    if(length(data1) <= 1 || length(data2) <= 1 || ( length(unique(data1)) == 1 & length(unique(data2)) == 1) ) is_good_run = 0 # 两组数据数量都要大于1，且两组数据中的唯一数据不能同时只有一个

    pvalue = NA
    tvalue = NA
    if(is_good_run == 1)
    {
        if(model == "dependent"){
            result = tryCatch( {t.test(data1, data2, alternative = "two.sided", paired = TRUE)},
                    error = function(e){ cat("ERROR :", conditionMessage(e), "\n")}
                    )
        }else{
            result = tryCatch( {t.test(data1, data2, alternative = "two.sided")},
                    error = function(e){ cat("ERROR :", conditionMessage(e), "\n")}
                    )
        } 

        if(!is.null(result)) # 如果报错了，result是NULL
        {
            pvalue <- result$p.value
            tvalue <- result$statistic
        }
    }

    return(c(pvalue, NA, tvalue)) # 保留NA因为要加FDR矫正值
}

# FoldChange子函数
run_foldchange <- function(data1, data2)
{
    # 差异倍数计算， 返回 foldchange
    if(sum(data1) != 0 && sum(data2) != 0){
        foldchange_1 <- mean(data1)
        foldchange_2 <- mean(data2)
        foldchange   <- log2(foldchange_1/foldchange_2)
    }
    else{
        cat('[Error] please check input file !\n')
        q()
    }
    return(c(foldchange))
}

run_type <- function(P_value, log2FC)
{
    type <- 'Not DEG'
    if(P_value <= 0.05 & log2FC >= LOG2FC){
        type <- "Up"
    }
    if(P_value <= 0.05 & log2FC <= -(LOG2FC)){
        type <- "Down"
    }
    return(c(type))
}

################################################# 计算主体 #################################################

# 读取文件
data          <- read.table(inputfile, header = T, sep = "\t", comment.char = "", check.names = F, quote = "") # 列为样本，行为特征值
sample.groups <- read.table(groupfile, header = F, sep = "\t", quote = "", colClasses = 'character')                         # 读取需要分析样本以及分组名称

# 获取分组文件数据
sample_need  = as.character(sample.groups[, 1])
sample_group = as.character(sample.groups[, 2])

# 分组文件检测
if(sum(sample_need %in% colnames(data)) != length(sample_need) )
{
    cat('[Error] 分组文件中，部分样本编号在丰度文件中不存在! 请仔细检查分组文件！\n')
    q()
}

data_need    = data[, sample_need, drop = FALSE]
group        = unique(sample_group)

if(length(group) > 2 || '' %in% group || NA %in% group)
{
    cat('[Error] group num > 2 or exist missing value, please check groupfile !\n')
    q()
}

# 获取原始文件中相关分组数据
case_data    = data[, as.character(sample.groups[, 1][sample.groups[, 2] == group[1]]), drop = FALSE]
case_data    = as.matrix(case_data)
control_data = data[, as.character(sample.groups[, 1][sample.groups[, 2] == group[2]]), drop = FALSE]
control_data = as.matrix(control_data)

# 获取需要保留的列编号
keep_column = as.numeric(unlist(strsplit(keep_column_list, split = ",")))

# 结果列表，空矩阵
result_anno         <- matrix(, nrow = nrow(data), ncol = length(keep_column))  # 保留的特征值等注释
result_sample       <- matrix(, nrow = nrow(data), ncol = length(data_need))    # 保留的分析样本表达值
result_ttest        <- matrix(, nrow = nrow(data), ncol = 3)                    # ttest结果汇总
result_foldchange   <- matrix(, nrow = nrow(data), ncol = 1)                    # foldchange结果
result_type         <- matrix(, nrow = nrow(data), ncol = 1)                    # 上下调结果

# 添加表头

colnames(result_anno)         <- colnames(data)[keep_column]
colnames(result_ttest)        <- c("P_value", "FDR P_value", "t_value")
colnames(result_foldchange)   <- c("log2FoldChange")
colnames(result_type)         <- c("type")

# 开始计算
for(row in 1:nrow(data))
{ 
    case_data_1    = case_data[row, ]
    control_data_1 = control_data[row, ]
    # Ttest 配对模式
    if(model == "dependent")
    {
        analysis_data        <- data.frame(case = case_data_1, control = control_data_1)
        analysis_data        <- analysis_data[complete.cases(analysis_data), ] # 去掉缺失值
        case_data_1          <- analysis_data$case
        control_data_1       <- analysis_data$control
    }else{
        case_data_1          <- case_data_1[complete.cases(case_data_1)]
        control_data_1       <- control_data_1[complete.cases(control_data_1)]
    }

    # 数据分析
    result_anno[row, ]       <- as.matrix(data[row, keep_column])[1, ]

    result_ttest[row, ]      <- run_ttest(case_data_1, control_data_1, model)

    # #foldchange统计
    result_foldchange[row, ] <- c(run_foldchange(case_data_1, control_data_1))

    #上下调统计
    result_type[row, ]       <- c(run_type(result_ttest[row, 1], result_foldchange[row, ]))

}

#获取所有样本的列数据
result_sample = data[, sample_need, drop = F]
#FDR矫正
result_ttest[, 'FDR P_value']  <- p.adjust(result_ttest[, 'P_value'], method = 'fdr')

# 结果矩阵汇总输出
result_final = cbind(result_anno, result_ttest, result_foldchange, result_type, result_sample)

write.table(result_final, output , sep = "\t", quote = F, row.names = F, col.names = T)

# 差异ID
diff_id = result_final[,1][result_final$P_value < 0.05]
write.table(diff_id, paste0(output, ".diff_id.txt"), sep = "\t", quote = F, row.names = F, col.names = F)


################################################ 绘图 #################################################

result_final$type = factor(result_final$type, levels = c('Up', 'Down', 'Not DEG'), order = T)

data_volcano <- data.frame(logFC = result_final$log2FoldChange, padj = -1*log10(result_final$P_value), type = result_final$type)

# 绘制火山图
x_lim <- max(5, -5)   # 选最大值作为xlim的上下边界
pdf(file = volcano, width = 10, height = 10)
theme_set(theme_bw())

p <- ggplot(data_volcano, aes(logFC, padj, color = type)) +
    geom_point() +
    xlim(-x_lim, x_lim) +
    ylim(0, 15) +
    labs(x = "log2(FoldChange)", y = "-log10(P_value)", title = "Volcanoplot") +
    scale_color_manual(values = c('Up' = "#BC3C28", 'Not DEG' = "grey", 'Down' = "#0072B5")) +  # 上下调颜色
    geom_hline(yintercept = -log10(0.05), linetype = 4) +                                       # 阈值虚线
    geom_vline(xintercept = c(-LOG2FC, LOG2FC), linetype = 4) +
    theme(panel.grid = element_blank(), legend.title = element_blank()) +                       #不显示legend的title
    theme(plot.title = element_text(size = 20, face = 'bold', vjust = 0.5, hjust = 0.5), axis.line = element_line(size = 0)) +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))
p
dev.off()

message("finish metabolite_diff_analysis.r")
