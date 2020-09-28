#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)
"Usage: differ_analysis_normal.r -i <file> -o <file> --sample_group <file> --case_group_name <string> --control_group_name <string>  [--paired  -k <column num>]
Options:
   -i , --input <file>                   输入文件，每一行是一个特征，每一列是一个样本，允许缺失(脚本会自动剔除)。第一列是特征名，第一行是样本名
   -o , --output <file>                  输出文件
   --sample_group <file>                 样本分组文件，两列数据，第一列样本名，第二列样本分组，有表头。注意，如果是paired模式，务必保证该文件中的case/control样本先后顺序一致。
   --case_group_name <string>            case分组名称
   --control_group_name <string>         control分组名称

   -p , --paired                         TTEST 模型，样本是否配对
   -k , --keep_column <column num>       输入文件中的哪几列放入结果里？即追加input中的注释内容到结果文件里，多列用逗号分隔  example: -k 1,2,3  [default: 1]
" -> doc

opts                      <- docopt(doc, version = '差异分析 ttest/utest/logistic \n          甘斌 129\n')
input                     <- opts$input
output                    <- opts$output
sample_group              <- opts$sample_group
case_group_name           <- opts$case_group_name
control_group_name        <- opts$control_group_name
paired                    <- opts$paired
keep_column_list          <- opts$k


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

    x         <- c(data1, data2)
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
data_input <- read.table(input, header=T, sep="\t", comment.char="", check.names=F)

# 读入分组文件
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
if(sum(!sampleinfo$sample %in% colnames(data_input)) > 0)
{
    losts = sampleinfo$sample[!sampleinfo$sample %in% colnames(data_input)]
    message("[Error] sample_group中的样本在input文件中缺失:", losts)
    q()
}

case_samples      = sampleinfo$sample[which(sampleinfo$group == case_group_name)]
control_samples   = sampleinfo$sample[which(sampleinfo$group == control_group_name)]

if(paired & length(case_samples) != length(control_samples))
{
    message("[Error] 分析要求为paired模式，但是case组样本数量与control组样本数量不一致")
    q()
}

# 获取分组数据
case_data    = data_input[, case_samples, drop = FALSE]
case_data    = as.matrix(case_data)
control_data = data_input[, control_samples, drop = FALSE]
control_data = as.matrix(control_data)

# 获取需要保留的列编号
keep_column = as.numeric(unlist(strsplit(keep_column_list, split=",")))


# 结果列表，空矩阵
result_anno      <- matrix(, nrow=nrow(data_input), ncol=length(keep_column))  # 保留的注释
result_ttest     <- matrix(, nrow=nrow(data_input), ncol=3)  # ttest结果汇总
result_utest     <- matrix(, nrow=nrow(data_input), ncol=3)  # utest结果汇总
result_logistic  <- matrix(, nrow=nrow(data_input), ncol=8)  # 逻辑回归结果汇总
result_other     <- matrix(, nrow=nrow(data_input), ncol=21) # 其他数据描述信息汇总

# 添加表头
colnames(result_anno)      <- colnames(data_input)[keep_column]
colnames(result_ttest)     <- c("P-value(Ttest)", "FDR P-value(Ttest)", "t-value(Ttest)")
colnames(result_utest)     <- c("P-value(Utest)", "FDR P-value(Utest)", "Z-value(Utest)")
colnames(result_logistic)  <- c("P-value(Logistic)", "FDR P-value(Logistic)", "Estimate(Logistic)", "OR(L95-U95)(Logistic)", "Sensitivity(Logistic)", "Specificity(Logistic)", "Accuracy(Logistic)", "AUC(Logistic)")

other_anno_names = c('Count', 'Mean', 'Median', 'StdDev', 'Min', 'Max', 'Quantile_25%', 'Quantile_75%', 'Shapiro-Wilk_Wvalue', 'Shapiro-Wilk_Pvalue')
colnames(result_other)        <- c('GroupDiff', paste(case_group_name, other_anno_names, sep = '_'), paste(control_group_name, other_anno_names, sep = '_'));
 
# 开始计算
for(row in 1:nrow(data_input))
{ 
    case_data_1 = case_data[row, ]
    control_data_1 = control_data[row, ]
 
    # 数据配对处理
    if(paired) 
    {   
        analysis_data <- data.frame(case=case_data_1, control=control_data_1)
        analysis_data <- analysis_data[complete.cases(analysis_data), ] # 去掉缺失值
        case_data_1   <- analysis_data$case
        control_data_1 <- analysis_data$control
    }else
    {
        case_data_1 <- case_data_1[complete.cases(case_data_1)]
        control_data_1 <- control_data_1[complete.cases(control_data_1)]
    }

    # 数据分析
    result_anno[row, ] <- as.matrix(data_input[row, keep_column])[1, ]
    result_ttest[row, ]     <- run_ttest(case_data_1,     control_data_1, paired) # TTEST 非配对模式
    result_utest[row, ]     <- run_utest(case_data_1,     control_data_1, paired) # Utest 非配对模式
    result_logistic[row, ]  <- run_logistic(case_data_1, control_data_1)     # 逻辑回归
 

    # 其他统计
    case_desc            <- run_vector_desc(case_data_1)                  # 数据量统计
    control_desc         <- run_vector_desc(control_data_1)               # 数据量统计
    diff                 <- case_desc[2] - control_desc[2]                # 组间平均甲基化差异
    case_normal_desc     <- run_normal_test(case_data_1)                  # 正态分布检验
    control_normal_desc  <- run_normal_test(control_data_1)               # 正态分布检验
    result_other[row, ] <- c(diff, case_desc, case_normal_desc, control_desc, control_normal_desc)
}

# FDR矫正
result_ttest[, 'FDR P-value(Ttest)']       <- p.adjust(result_ttest[, 'P-value(Ttest)'], method = 'fdr')
result_utest[, 'FDR P-value(Utest)']       <- p.adjust(result_utest[, 'P-value(Utest)'], method = 'fdr')
result_logistic[, 'FDR P-value(Logistic)'] <- p.adjust(result_logistic[, 'P-value(Logistic)'], method = 'fdr')

# 结果矩阵汇总输出
result_final = cbind(result_anno, result_ttest, result_utest, result_logistic, result_other)

write.table(result_final, output , sep = "\t", quote = F, row.names = F, col.names = TRUE )
