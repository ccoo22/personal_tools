#!/home/genesky/software/r/4.0.3/bin/Rscript

library(docopt)

"Usage: logistics_single_factor.r  -i <file> -p <prefix>  -c <file> --interaction <string> [ --rlib <dir>]

Options:
    -i, --input <file>      输入文件，第一列为配对样本信息，配对样本编号相同，从1开始编号，依次类推，缺失数据所在的样本行会被删除；第二列为样本名，第三列为样本分组（case=1，control=0），往后为表型列，缺失数据所在的样本行会被删除。该软件对每一个基因单独做逻辑回归分析。 列名不能包含空格、“-”等特殊字符
    -c, --cov <file>        协变量输入文件，用于矫正。第一列样本名，后面所有列都是协变量。样本顺序没有要求。默认不用输入。 务必保证cov中包含input的所有样本，cov文件中的所有列都用于矫正。列名不能包含空格、“-”等特殊字符
    --interaction <string>  进行交互作用分析的协变量名称，即--cov文件中某一列矫正变量的列名称
    -p, --prefix <prefix>   输出文件前缀， 生成: (1) ROC曲线图 prefix.ROC.pdf  (2) 逻辑回归模型结果 prefix.model.txt
    --rlib <dir>            R包路径 [default: /home/genesky/software/r/4.0.3/lib64/R/library]" -> doc

opts   <- docopt(doc, version='条件逻辑回归、ROC \n')
input       <- opts$input
cov         <- opts$cov
interaction <- opts$interaction
prefix      <- opts$prefix
rlib        <- opts$rlib

if(!is.null(rlib)) .libPaths(rlib)

 
library(ROCR)
library(pROC)
library(Hmisc)
library(survival)

set.seed(91)

##测试数据
# input = "/home/xuy/work_new/19B0819A-3_new/data/smoke/input_smoke.txt"
# cov   = "/home/xuy/work_new/19B0819A-3_new/data/smoke/cov/cov_smoke.txt"
# prefix = "update2"
# interaction = 'smoke'

# 读入 input
data_input = read.table(input, head = T, check.names = F, sep="\t")
colnames(data_input)[1:3] = c('matset', 'sample', 'group')
genes <- colnames(data_input)[4:ncol(data_input)]

if(sum(grepl(" ", colnames(data_input))) > 0 |  sum(grepl("-", colnames(data_input))))
{
    message("[Error]  输入的input文件表头包含空格或者-\n")
    q()
}

# 读入cov(协变量)
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
col_names = c('Var', 'NMISS case', 'NMISS control', "Pvalue", "Pvalue_FDR", "coef", "exp(coef)",  "se(coef)", "z", "lower .95", "upper .95")  # 必要结果
#### cov_header_names 存储交互作用分析的协变量、交互作用分析项
##showvars 存储最后要展示的变量结果：gene、交互作用分析的协变量、gene:交互作用分析的协变量
cov_header_names = c(interaction, paste("interaction with", interaction)) 
result <- matrix(ncol = 20 + 8 * length(cov_header_names), nrow = length(genes))

for(cov_name in cov_header_names)  # 交互作用分析的协变量、交互作用分析项结果
{
    col_names = c(col_names, paste0(c("Pvalue", "Pvalue_FDR", "coef", "exp(coef)",  "se(coef)", "z", "lower .95", "upper .95"), "(", cov_name, ")") )
}
col_names = c(col_names, "Likelihood_ratio_test_Pvalue", "Wald_test_Pvalue", "Score_(logrank)_test_Pvalue", "AUC", "AUC_CI95", 'threshold', "specificity", "sensitivity", "accuracy")  # 模型结果

colnames(result) = col_names

pdf_file = paste0(prefix, '.ROC.pdf')
pdf(pdf_file)

for( col in 1:length(genes))
{   
    gene = genes[col]
    vars = c(gene)  # 变量列表
    
    interaction_model_name = paste(gene, ":", interaction, sep='')
    message("process : ", gene)
    # 准备输入数据、cov数据
    data_tmp = data_input[, c('matset', 'sample', 'group', gene)]
    if(! is.null(cov))
    {
        data_tmp = cbind(data_tmp, data_cov[data_tmp[,'sample'], ])
        ##showvars 存储最后要展示的变量结果：gene、交互作用分析的协变量、gene:交互作用分析的协变量
        showvars = c(vars, interaction, interaction_model_name)
        vars = c(vars, cov_names) 
    }
    data_tmp = data_tmp[complete.cases(data_tmp), ]  # 去掉缺失值
    colnames(data_tmp) = c('matset', 'sample', 'group', gene, cov_names)

    # 检查数据是否异常
    nmiss_case    = sum(data_tmp$group == 1)
    nmiss_control = sum(data_tmp$group == 0)
    if(nmiss_case == 0 | nmiss_control == 0 | sum(!duplicated(data_tmp[, 4])) == 1)
    {   
        message("[Error] 数据错误，存在的缺失太多或者x值只有一种数据，跳过, skip. case = ", nmiss_case, ' . control = ', nmiss_control)
        result[col, 1:3] = c(gene, nmiss_case, nmiss_control) 
        next
    }   

    # 逻辑回归  
    formula    <- as.formula(paste('group ~ ', paste(vars, collapse = ' + '), '+ strata(matset)', '+', paste(vars[1], interaction, sep = '*'))) ## Y ~ x1 + x2 + x1*x2 + x3 + ...
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
        result_value = c(gene, nmiss_case, nmiss_control) 

        for(cov_name in showvars)  # 协变量结果
        {   
            result_value = c(result_value, clogit_result[cov_name, 'Pr(>|z|)'], NA, clogit_result[cov_name, c('coef', 'exp(coef)', 'se(coef)', 'z')], clogit_conf[cov_name, c("lower .95", "upper .95")]) 
        }
        result_value = c(result_value, test_pvalues, AUC, AUC_CI95, thr, spe, sen, acc)  # 模型结果
        
        result[col, ] <- result_value
    }else{
        result[col, 1:3] = c(gene, nmiss_case, nmiss_control)
    } 
}
dev.off()
result[,'Pvalue_FDR'] =  p.adjust(result[,'Pvalue'], method = "fdr", n=nrow(result))
for(cov_name in cov_header_names)  # 协变量结果
{   
    pvalue_header = paste0("Pvalue", "(", cov_name, ")")
    pvalue_fdr_header = paste0("Pvalue_FDR", "(", cov_name, ")")  
    result[,pvalue_fdr_header] =  p.adjust(result[,pvalue_header], method = "fdr", n=nrow(result))
}


# 结果输出
model_file = paste0(prefix, '.model.txt')
write.table(result, model_file, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )

message("ROC:              ", pdf_file)
message("logistic model:   ", model_file)

