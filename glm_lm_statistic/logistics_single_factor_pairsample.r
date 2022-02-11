#!/home/genesky/software/r/4.0.3/bin/Rscript

library(docopt)

"Usage: logistics_single_factor.r  -i <file> -p <prefix> [ -c <file> --rlib <dir>]

Options:
    -i, --input <file>      输入文件，第一列为配对样本信息，配对样本编号相同，从1开始编号，依次类推，缺失数据所在的样本行会被删除；第二列为样本名，第三列为样本分组（case=1，control=0），往后为表型列，缺失数据所在的样本行会被删除。该软件对每一个基因单独做逻辑回归分析。 列名不能包含空格、“-”等特殊字符
    -c, --cov <file>        协变量输入文件，用于矫正。第一列样本名，后面所有列都是协变量。样本顺序没有要求。默认不用输入。 务必保证cov中包含input的所有样本，cov文件中的所有列都用于矫正。列名不能包含空格、“-”等特殊字符
    -p, --prefix <prefix>   输出文件前缀， 生成: (1) ROC曲线图 prefix.ROC.pdf  (2) 逻辑回归模型结果 prefix.model.txt
    --rlib <dir>            R包路径 [default: /home/genesky/software/r/4.0.3/lib64/R/library]
                            注意：有时会报错 Error in fitter(X, Y, istrat, offset, init, control, weights = weights, 外接函数调用时不能有NA/NaN/Inf(arg5)，这表明模型不合适，忽略即可。
    " -> doc

opts   <- docopt(doc, version='条件逻辑回归、ROC \n')
input       <- opts$input
cov         <- opts$cov
prefix      <- opts$prefix
rlib        <- opts$rlib

if(!is.null(rlib)) .libPaths(rlib)

 
library(ROCR)
library(pROC)
library(Hmisc)
library(survival)

set.seed(91)

##测试数据
# input = "/home/dongxj/work/research/logistics/v3/input_gene.txt"
# cov   = "/home/dongxj/work/research/logistics/v3/cov_2.txt"
# prefix = "update"

# 读入 input
data_input = read.table(input, head = T, check.names = F, sep="\t")
colnames(data_input)[1:3] = c('matset', 'sample', 'group')
genes <- colnames(data_input)[4:ncol(data_input)]

if(sum(grepl(" ", colnames(data_input))) > 0 |  sum(grepl("-", colnames(data_input))))
{
    message("[Error]  输入的input文件表头包含空格或者-\n")
    q()
}

# 读入cov
data_cov = NULL
cov_names = c()
if(! is.null(cov))
{
    data_cov = read.table(cov, head = T, row.names = 1, check.names = F, sep="\t")
    if(sum(grepl(" ", colnames(data_cov))) > 0 |  sum(grepl("-", colnames(data_cov))))
    {
        message("[Error]  输入的cov文件表头包含空格或者-\n")
        q()
    }
    if(sum(data_input[,"sample"] %in% row.names(data_cov)) != nrow(data_input) )
    {
        message("[Error]  输入的cov文件没有包含所有的input样本\n")
        q()  
    }
    cov_names = colnames(data_cov)
}

# 保留模型值
result <- matrix(ncol = 37 + 2 * length(cov_names), nrow = length(genes))
col_names = c('Var', 'NMISS case', 'NMISS control', "Pvalue", "Pvalue_FDR", "Likelihood_ratio_test_Pvalue", "Wald_test_Pvalue", "Score_(logrank)_test_Pvalue", "coef", "exp(coef)",  "se(coef)", "z", "lower .95", "upper .95")  # 必要结果
for(cov_name in cov_names)  # 协变量结果
{
    col_names = c(col_names, paste0(c("Pvalue", "exp(coef)"), "(", cov_name, ")") )
}
col_names = c(col_names, "AUC", "AUC_CI95", 'threshold', "specificity", "sensitivity", "accuracy")  # 模型结果

for (i in c('Case', 'Control')){
    col_names = c(col_names, paste0(i,c('_Count', '_Mean', '_Median', '_StdDev', '_Min', '_Max', '_Quantile_25%', '_Quantile_75%'))) ##各统计量结果
}

colnames(result) = c(col_names, 'group_diff')

pdf_file = paste0(prefix, '.ROC.pdf')
pdf(pdf_file)

for( col in 1:length(genes))
{   
    gene = genes[col]
    vars = c(gene)  # 变量列表

    message("process : ", gene)
    # 准备输入数据、cov数据
    data_tmp = data_input[, c('matset', 'sample', 'group', gene)]
    if(! is.null(cov))
    {
        data_tmp = cbind(data_tmp, data_cov[data_tmp[,'sample'], ])
        vars = c(vars, cov_names)
    }
    data_tmp = data_tmp[complete.cases(data_tmp), ]  # 去掉缺失值
    matset_count = as.data.frame(table(data_tmp$matset)) ## table()函数统计每个配对编号下样本的数量
    data_tmp = data_tmp[data_tmp$matset %in%  as.numeric(as.character(matset_count[which(matset_count$Freq==2),'Var1'])),]  ## 检查matset列，去掉不配对的样本 


    # 检查数据是否异常
    nmiss_case    = sum(data_tmp$group == 1)
    nmiss_control = sum(data_tmp$group == 0)
    if(nmiss_case == 0 | nmiss_control == 0 | sum(!duplicated(data_tmp[, 4])) == 1)
    {   
        message("[Error] 数据错误，存在的缺失太多或者x值只有一种数据，跳过, skip. case = ", nmiss_case, ' . control = ', nmiss_control)
        result[col, 1:3] = c(gene, nmiss_case, nmiss_control) 
        next
    }

    #计算case组、control组的均值、标准差、group_diff等统计量
    # 数据量统计函数
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
    ## 数据准备、计算
    case_data = data_tmp[data_tmp$group == 1, gene]
    control_data = data_tmp[data_tmp$group == 0, gene]
    case_statistics = run_vector_desc(case_data)
    control_statistics = run_vector_desc(control_data)
    group_diff = case_statistics[2] - control_statistics[2]

    # 逻辑回归
    formula    <- as.formula(paste('group ~ ', paste(vars, collapse = ' + '), '+ strata(matset)' ) )
    clogit_fit = tryCatch( {clogit(formula, data = data_tmp, method = "exact") }, 
        error = function(e){ cat("ERROR :",conditionMessage(e), "\n")} 
        ) 
    if(!is.null(clogit_fit)) # 如果报错了，result是NULL
    {
        clogit_summary <- summary(clogit_fit)  ##str(clogit_summary) 查看结果
        clogit_result  <- clogit_summary$coefficients
        clogit_conf    <- clogit_summary$conf.int
        modelroc = roc(data_tmp$group, predict(clogit_fit, type='lp'), transpose = FALSE, quiet = TRUE) ##?predict.coxph查看帮助文档，type类型：“lp”, “risk”, “expected”, “terms”, “survival”
        best = coords(modelroc, "best", ret = c("threshold", "specificity", "sensitivity", "accuracy"))
        best = as.matrix(best)
        AUC  = round(modelroc$auc[[1]], 4)
        thr  = round(best[1, "threshold"], 4)
        spe  = round(best[1, "specificity"], 4)
        sen  = round(best[1, "sensitivity"], 4)
        acc  = round(best[1, "accuracy"], 4)
        CI   = signif(ci.auc(modelroc, conf.level=0.95, method=c("delong", "bootstrap")), 4)
        L95 = CI[1]
        U95 = CI[3]
        AUC_CI95 = paste(L95, "-", U95)

        # 绘图
        prob <- predict(clogit_fit, newdata = data_tmp, type = "lp")
        pred <- prediction(prob, data_tmp$group)
        perf = performance(pred, "tpr", "fpr")
        par(family="Times")  # times new roman 字体
        plot(perf@x.values[[1]], perf@y.values[[1]], lwd = 3, col = "#2171B5", type = "l", bty = "n", axes = FALSE, ann = FALSE)
        abline(a = 0, b = 1, lty = 5, lwd = 3, col = "#D94801" )
        par(tck = 0.03, col.axis = "#084594",lwd = 1, font.axis = 2, cex.axis = 1 )
        axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), col =  "#8C2D04", lwd = 2)
        axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), col =  "#8C2D04", lwd = 2)
        minor.tick(nx = 4, ny = 4, tick.ratio = 0.01)
        box(which = "plot", col = "#8C2D04", lwd = 3)
        title(xlab = 'False positive rate', ylab = 'True positive rate',  main = gene, col.lab = "#084594", col.main = "#084594")
        text (0.4, 0.3, c(paste("AUC:", AUC)), adj = 0, font = 2)
        text (0.4, 0.2, c(paste("95%CI:", L95, "-", U95)), adj = 0, font = 2)
    
        # 结果记录
        test_pvalues = c(clogit_summary$logtest[3], clogit_summary$waldtest[3], clogit_summary$sctest[3])
        result_value = c(gene, nmiss_case, nmiss_control, clogit_result[gene, 'Pr(>|z|)'], NA, test_pvalues, clogit_result[gene, c('coef', 'exp(coef)', 'se(coef)', 'z')], clogit_conf[gene, c("lower .95", "upper .95")])
        for(cov_name in cov_names)  # 协变量结果
        {   
            result_value = c(result_value, clogit_result[cov_name, c('Pr(>|z|)', 'exp(coef)')]) 
        }
        result_value = c(result_value, AUC, AUC_CI95, thr, spe, sen, acc)  # 模型结果
        
        result[col, ] <- c(result_value, case_statistics, control_statistics, group_diff)
    }else{
        result[col, 1:3] = c(gene, nmiss_case, nmiss_control)
    } 
}
dev.off()
result[,'Pvalue_FDR'] =  p.adjust(result[,'Pvalue'], method = "fdr", n=nrow(result))


# 结果输出
model_file = paste0(prefix, '.model.txt')
write.table(result, model_file, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )

message("ROC:              ", pdf_file)
message("logistic model:   ", model_file)

