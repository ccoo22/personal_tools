#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)
"Usage: differ_analysis_normal.r -i <file> -o <file>  --case_group_name <string> --control_group_name <string> --case_group_sample_file <file> --control_group_sample_file <file>  [--paired  -k <column num>]
Options:
   -p , --paired                         TTEST 模型，样本是否配对
   -k , --keep_column <column num>       keep input column to output, example: -k 1,2,3  [default: 1]
   -i , --input <file>                   the input file, each row is a point, each column is sample
   -o , --output <file>                  the output file
   --case_group_name <string>            the case group name
   --control_group_name <string>         the control group name
   --case_group_sample_file <file>       case组样本列表文件，一列数据，没有表头，第一列就是样本名，一行一个样本
   --control_group_sample_file <file>    control组样本列表文件，一列数据，没有表头，第一列就是样本名，一行一个样本" -> doc

opts                      <- docopt(doc, version = 'Program : methylation diff analysis Ttest/Utest/Logistic v1.1 \n          甘斌 129\n')
input                     <- opts$input
case_group_name           <- opts$case_group_name
control_group_name        <- opts$control_group_name
case_group_sample_file    <- opts$case_group_sample_file
control_group_sample_file <- opts$control_group_sample_file
# case_group_sample_list    <- opts$case_group_sample_list
# control_group_sample_list <- opts$control_group_sample_list
output                    <- opts$output
paired                    <- opts$paired
keep_column_list          <- opts$k
 
case_sample_name <- read.table(case_group_sample_file, header=F, sep="\t", comment.char="", check.names=F, stringsAsFactors=F)
case_group_sample_list = paste(case_sample_name[,1], collapse = ',')
control_sample_name <- read.table(control_group_sample_file, header=F, sep="\t", comment.char="", check.names=F, stringsAsFactors=F)
control_group_sample_list = paste(control_sample_name[,1], collapse = ',')
# 检测用
# cat(input, "\n")
# cat(case_group_name, "\n")
# cat(control_group_name, "\n")
# cat(case_group_sample_list, "\n")
# cat(control_group_sample_list, "\n")
# cat(output, "\n")
# cat(model, "\n")
# cat(keep_column_list, "\n")
# q()


# 测试用参数
# input                     <- '/home/ganb/work/research/methyltarget_diff_program/19B0314A/report/methylation_result/2.methyl.site.txt'
# case_group_name           <- 'Treat'
# control_group_name        <- 'Control'
# case_group_sample_list    <- '016_24h,30_24h,41_24h,47_24h,61_24h,46_24h'
# control_group_sample_list <- '016_0h,30_0h,41_0h,47_0h,61_0h,46_0h'
# output                    <- './test.txt'
# model                     <- 'independent'
# keep_column_list          <- '1,2,6'
 
options(warn=-1)
library(pROC, quietly = TRUE) # roc要用
library(fBasics, quietly = TRUE) # 正态分布检测要用
# TTEST子函数
run_ttest <- function(data1, data2, paired)
{
    # t检验， 返回 Pvalue Tvalue
    is_good_run = 1
    if(length(data1) <= 1 || length(data2) <= 1 || ( length(unique(data1)) == 1 & length(unique(data2)) == 1) ) is_good_run = 0 # 两组数据数量都要大于1，且两组数据中的唯一数据不能同时只有一个

    pvalue = NA
    tvalue = NA
    if(is_good_run == 1)
    {

        result = tryCatch( {t.test(data1, data2, alternative="two.sided", paired=paired)}, 
                    error = function(e){ cat("ERROR :",conditionMessage(e), "\n")} 
                    ) 

        if(!is.null(result)) # 如果报错了，result是NULL
        {
            pvalue <- result$p.value
            tvalue <- result$statistic
        }
           
    }

    return(c(pvalue, NA, tvalue))  # 保留的NA是为了后面加FDR矫正值
}

# UTEST子函数
run_utest <- function(data1, data2, paired)
{
    # u检验， 返回 Pvalue Zvalue
    pvalue <- NA
    zvalue <- NA
    if(sum(data1) != 0 && sum(data2) != 0){
        result <- wilcox.test(data1, data2, paired=paired, exact = NULL, correct = TRUE)
        pvalue <- result$p.value
        zvalue <- qnorm(result$p.value)
    }

    return(c(pvalue, NA, zvalue))  # 保留的NA是为了后面加FDR矫正值
}

# 逻辑回归子函数
run_logistic <- function(data1, data2)
{
    # 逻辑回归 data1 = case, data2 = control
    # 返回：P值、or95%CI、敏感性、特异性、准确性、AUC
    if(length(data1) == 0 || length(data2) == 0 || length(c(data1, data2)) == 2 ) return (c(NA, NA, NA, NA, NA, NA, NA, NA)) 

    x         <- c(data1 * 100, data2 * 100)
    y         <- c(rep(1, length(data1)), rep(0, length(data2)))
    data_logi <- data.frame(x=x, y=y)
    glm.fit   <- glm(y~x, data=data_logi, family=binomial(link='logit'))
    
    result_glm <- summary(glm.fit)$coefficient
    result_glm <- cbind(result_glm, or = exp(result_glm[, 1])) # 添加or值
    
    # 可能会出现部分数据无法计算的情况，为了通用性，这里补充完整
    if( ("(Intercept)" %in% rownames(result_glm) ) == FALSE) result_glm <- rbind(result_glm, "(Intercept)"=NA) 
    if( ("x" %in% rownames(result_glm)) == FALSE ) result_glm <- rbind(result_glm, "x"=NA) 
    
    # y完全一样不能做roc
    if(diff(range(y)) == 0) 
    {
        AUC  <- NA
        spec <- NA
        sens <- NA
        accu <- NA
    }else{
        modelroc <- roc(y, predict(glm.fit, type='response'), quiet = TRUE)
        best <- coords(modelroc, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy")) # 返回最佳阈值结果，注意可能会有多个结果，后面默认只取第一个
        best <- as.matrix(best)
        AUC  <- modelroc$auc[[1]]
        spec <- best[1, "specificity"]
        sens <- best[1, "sensitivity"]
        accu <- best[1, "accuracy"]
    }
    
    pvalue_logi <- result_glm["x", "Pr(>|z|)"]
    est_logi    <- result_glm["x", "Estimate"]
    se_logi     <- result_glm["x", "Std. Error"]
    or_logi     <- result_glm["x", "or"]            
    l95_logi    <- exp(est_logi - 1.96 * se_logi)
    u95_logi    <- exp(est_logi + 1.96 * se_logi)
    or_ci_logi  <- paste0(round(or_logi, 4), "(", round(l95_logi, 4), "-", round(u95_logi, 4), ")")

    return(c(pvalue_logi, NA, est_logi, or_ci_logi, sens, spec, accu, AUC))  # 保留的NA是为了后面加FDR矫正值

}

# 数据量统计
run_vector_desc <- function(data)
{
    # 数据量统计
    count  <- length(data)
    mean   <- mean (data)
    median <- median (data)
    sd     <- sd (data)
    min    <- min (data)
    max    <- max (data)
    q25     <- quantile (data, .25)[[1]]
    q75     <- quantile (data, .75)[[1]]   

    return(c(count, mean, median, sd, min, max, q25, q75)) 
}

# 正态分布检测
run_normal_test <- function(data)
{
    
    # Shapiro-Wilk's正态检验,数据量必须在3-5000,并且数据不能全相等
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

################################################# 计算主体 #############################################################

# 读入数据,data.frame格式
data <- read.table(input, header=T, sep="\t", comment.char="", check.names=F)

# 获取分组数据
case_data    = data[, unlist(strsplit(case_group_sample_list, split=",")), drop = FALSE]
case_data    = as.matrix(case_data)
control_data = data[, unlist(strsplit(control_group_sample_list, split=",")), drop = FALSE]
control_data = as.matrix(control_data)

# 如果数据只有一行，会出问题
# case_data    <- sapply(unlist(strsplit(case_group_sample_list, split=",")), function(x){data[[x]] })
# control_data <- sapply(unlist(strsplit(control_group_sample_list, split=",")), function(x){data[[x]] })

# 获取需要保留的列编号
keep_column = as.numeric(unlist(strsplit(keep_column_list, split=",")))


# 结果列表，空矩阵
result_anno      <- matrix(, nrow=nrow(data), ncol=length(keep_column))  # 保留的注释
result_ttest     <- matrix(, nrow=nrow(data), ncol=3)  # ttest结果汇总
result_utest     <- matrix(, nrow=nrow(data), ncol=3)  # utest结果汇总
result_logistic  <- matrix(, nrow=nrow(data), ncol=8)  # 逻辑回归结果汇总
result_other     <- matrix(, nrow=nrow(data), ncol=21) # 其他数据描述信息汇总

# 添加表头
colnames(result_anno)      <- colnames(data)[keep_column]
colnames(result_ttest)     <- c("P-value(Ttest)", "FDR P-value(Ttest)", "t-value(Ttest)")
colnames(result_utest)     <- c("P-value(Utest)", "FDR P-value(Utest)", "Z-value(Utest)")
colnames(result_logistic)  <- c("P-value(Logistic)", "FDR P-value(Logistic)", "Estimate(Logistic)", "OR(L95-U95)(Logistic)", "Sensitivity(Logistic)", "Specificity(Logistic)", "Accuracy(Logistic)", "AUC(Logistic)")

other_anno_names = c('Count', 'Mean', 'Median', 'StdDev', 'Min', 'Max', 'Quantile_25%', 'Quantile_75%', 'Shapiro-Wilk_Wvalue', 'Shapiro-Wilk_Pvalue')
colnames(result_other)        <- c('GroupDiff', paste(case_group_name, other_anno_names, sep = '_'), paste(control_group_name, other_anno_names, sep = '_'));
 
# 开始计算
for(row in 1:nrow(data))
{ 
    case_data_1 = case_data[row, ]
    control_data_1 = control_data[row, ]
 
    # Ttest 配对模式
    if(paired)
    {   
        analysis_data <- data.frame(case=case_data_1, control=control_data_1)
        analysis_data <- analysis_data[complete.cases(analysis_data), ] # 去掉缺失值
        case_data_1_dep <- analysis_data$case
        control_data_1_dep <- analysis_data$control
    }

    # 去掉缺失值
    case_data_1 <- case_data_1[complete.cases(case_data_1)]
    control_data_1 <- control_data_1[complete.cases(control_data_1)]
 
    # 数据分析
 
    result_anno[row, ] <- as.matrix(data[row, keep_column])[1, ]
    if(!paired)    result_ttest[row, ]     <- run_ttest(case_data_1,     control_data_1,     paired) # TTEST 非配对模式
    if(paired)   result_ttest[row, ]     <- run_ttest(case_data_1_dep, control_data_1_dep,  paired) # TTEST 配对模式
    
    if(!paired)    result_utest[row, ]     <- run_utest(case_data_1,     control_data_1,     paired) # Utest 非配对模式
    if(paired)   result_utest[row, ]     <- run_utest(case_data_1_dep, control_data_1_dep, paired)# Utest 配对模式

    result_logistic[row, ]  <- run_logistic(case_data_1, control_data_1)     # 逻辑回归
 

    # 其他统计
    case_desc            <- run_vector_desc(case_data_1)                  # 数据量统计
    control_desc         <- run_vector_desc(control_data_1)               # 数据量统计
    meth_diff            <- case_desc[2] - control_desc[2]                # 组间平均甲基化差异
    case_normal_desc     <- run_normal_test(case_data_1)                  # 正态分布检验
    control_normal_desc  <- run_normal_test(control_data_1)          # 正态分布检验
    result_other[row, ] <- c(meth_diff, case_desc, case_normal_desc, control_desc, control_normal_desc)
}

# FDR矫正
result_ttest[, 'FDR P-value(Ttest)']       <- p.adjust(result_ttest[, 'P-value(Ttest)'], method = 'fdr')
result_utest[, 'FDR P-value(Utest)']       <- p.adjust(result_utest[, 'P-value(Utest)'], method = 'fdr')
result_logistic[, 'FDR P-value(Logistic)'] <- p.adjust(result_logistic[, 'P-value(Logistic)'], method = 'fdr')

# 结果矩阵汇总输出
result_final = cbind(result_anno, result_ttest, result_utest, result_logistic, result_other)

write.table(result_final, output , sep = "\t", quote = F, row.names = F, col.names = TRUE )
