#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: logistics_single_factor.r  -i <file> -p <prefix> [ -c <file> --maxit <integer> --rlib <dir>]

Options:
    -i, --input <file>      输入文件，第一列为样本名，第二列为样本分组（case=1，control=0），往后为基因列。该软件对每一个基因单独做逻辑回归分析。 列名不能包含空格、“-”等特殊字符
    -c, --cov <file>        协变量输入文件，用于矫正。第一列样本名，后面所有列都是协变量。样本顺序没有要求。默认不用输入。 务必保证cov中包含input的所有样本，cov文件中的所有列都用于矫正。列名不能包含空格、“-”等特殊字符
    -p, --prefix <prefix>   输出文件前缀， 生成: (1) ROC曲线图 prefix.ROC.pdf  (2) 逻辑回归模型结果 prefix.model.txt
    --maxit <integer>       逻辑回归最大迭代次数, 当保警告Warning: glm.fit: algorithm did not converge时，建议增加迭代次数 [default: 25]
    --rlib <dir>            R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，单因素逻辑回归、ROC \n')
input       <- opts$input
cov         <- opts$cov
prefix      <- opts$prefix
maxit       <- as.integer(opts$maxit)
rlib        <- opts$rlib

if(!is.null(rlib)) .libPaths(rlib)

 
library(ROCR)
library(pROC)
library(Hmisc)

set.seed(91)


# 读入 input
data_input = read.table(input, head = T, row.names = 1, check.names = F, sep="\t")
colnames(data_input)[1] <- 'y'
genes <- colnames(data_input)[2:ncol(data_input)]

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
    if(sum(row.names(data_input) %in% row.names(data_cov)) != nrow(data_input) )
    {
        message("[Error]  输入的cov文件没有包含所有的input样本\n")
        q()  
    }
    cov_names = colnames(data_cov)
}

# 保留模型值
result <- matrix(ncol = 13 + 2 * length(cov_names), nrow = length(genes))
col_names = c('Var', 'NMISS case', 'NMISS control', "Estimate", "Pvalue", "Pvalue_FDR", "Estimate(Intercept)", "Pvalue(Intercept)")  # 必要结果
for(cov_name in cov_names)  # 协变量结果
{
    col_names = c(col_names, paste0(c("Estimate", "Pvalue"), "(", cov_name, ")") )
}
col_names = c(col_names, "AUC", "AUC_CI95", 'threshold', "specificity", "sensitivity", "accuracy")  # 模型结果

colnames(result) = col_names

pdf_file = paste0(prefix, '.ROC.pdf')
pdf(pdf_file)

for( col in 1:length(genes))
{   
    gene = genes[col]
    vars = c(gene)  # 变量列表

    message("process : ", gene)

    # 准备输入数据、cov数据
    data_tmp = data_input[, c('y', gene)]
    if(! is.null(cov))
    {
        data_tmp = cbind(data_tmp, data_cov[row.names(data_tmp), ])
        vars = c(vars, cov_names)
    }

    data_tmp = data_tmp[complete.cases(data_tmp), ]  # 去掉缺失值

    # 检查数据是否异常
    nmiss_case    = sum(data_tmp$y == 1)
    nmiss_control = sum(data_tmp$y == 0)
    if(nmiss_case == 0 | nmiss_control == 0 | sum(!duplicated(data_tmp[, 2])) == 1)
    {   
        message("[Error] 数据错误，存在的缺失太多或者x值只有一种数据，跳过, skip. case = ", nmiss_case, ' . control = ', nmiss_control)
        result[col, 1:3] = c(gene, nmiss_case, nmiss_control) 
        next
    }   

    # 逻辑回归
    formula <- as.formula(paste('y ~ ', paste(vars, collapse = ' + ') ) )
    glm_fit <- glm(formula, data = data_tmp, family = binomial, control=list(maxit=maxit))
    fit_summary <- summary(glm_fit)
    glm_result  <- fit_summary$coefficients

    # ROC参数计算
    modelroc = roc(data_tmp$y, predict(glm_fit, type='response'), transpose = FALSE, quiet = TRUE)
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
    prob <- predict(glm_fit, newdata = data_tmp, type = "response")
    pred <- prediction(prob, data_tmp$y)
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
    result_value = c(gene, nmiss_case, nmiss_control, glm_result[gene, 'Estimate'], glm_result[gene, 'Pr(>|z|)'], NA, glm_result['(Intercept)', 'Estimate'], glm_result['(Intercept)', 'Pr(>|z|)'])
    for(cov_name in cov_names)  # 协变量结果
    {   
        result_value = c(result_value, glm_result[cov_name, 'Estimate'], glm_result[cov_name, 'Pr(>|z|)']) 
    }
    result_value = c(result_value, AUC, AUC_CI95, thr, spe, sen, acc)  # 模型结果
    
    result[col, ] <- result_value
}
dev.off()
result[,'Pvalue_FDR'] =  p.adjust(result[,'Pvalue'], method = "fdr", n=nrow(result))


# 结果输出
model_file = paste0(prefix, '.model.txt')

write.table(result, model_file, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )


message("ROC:              ", pdf_file)
message("logistic model:   ", model_file)

