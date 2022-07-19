#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: logistics_multi_factor.r  -i <file>  -p <prefix> [-g <string> -c <file> --maxit <integer>  --rlib <dir>]

Options:
    -i, --input <file>      输入文件，第一列为样本名，第二列为样本分组（case=1，control=0），往后为基因列。
    -g, --gene <string>     输入要纳入多元逻辑回归分析的基因名称，多个基因用逗号分隔， 默认是所有基因
    -c, --cov <file>        协变量输入文件，用于矫正。第一列样本名，后面所有列都是协变量。样本顺序没有要求。默认不用输入。 务必保证cov中包含input的所有样本，cov文件中的所有列都用于矫正。 
    -p, --prefix <prefix>   输出文件前缀， 生成: (1) ROC曲线图 prefix.ROC.pdf  (2) 逻辑回归模型结果 prefix.model.txt (3) 基于逻辑回归模型对输入样本做预测 prefix.predict.txt  (4) 最佳阈值结果 prefix.best.txt
    --maxit <integer>       逻辑回归最大迭代次数, 当保警告Warning: glm.fit: algorithm did not converge时，建议增加迭代次数 [default: 25]
    --rlib <dir>            R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，多因素逻辑回归、ROC \n')
input         <- opts$input
cov           <- opts$cov
prefix        <- opts$prefix
gene          <- opts$gene
maxit         <- as.integer(opts$maxit)
rlib          <- opts$rlib
genes <- c()
if(!is.null(gene))  genes <- unlist(strsplit(gene, ","))

if(!is.null(rlib)) .libPaths(rlib)
 
 
library(ROCR)
library(pROC)
library(Hmisc)

set.seed(91)


# 读入
data_input = read.table(input, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "")
genes_replace = data.frame(from=colnames(data_input)[2:ncol(data_input)], newname=paste('gene', 2:ncol(data_input), sep=''), stringsAsFactors=F)  # 基因名字替换，方式有特殊字符而无法分析
rownames(genes_replace) = genes_replace[, 'newname']
colnames(data_input) <- c('y', genes_replace[,'newname'])


if(is.null(gene)) 
{
    genes <- colnames(data_input)[2:ncol(data_input)]
}else{
    # 使用指定的基因名
    tmp = genes_replace
    rownames(tmp) = tmp[, 'from']
    genes = tmp[genes, 'newname']
}

 

# 读入cov
data_cov = NULL
cov_names = c()
cov_names_replace = data.frame()
if(! is.null(cov))
{
    data_cov = read.table(cov, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "")
    cov_names_replace = data.frame(from=colnames(data_cov), newname=paste('cov', 1:ncol(data_cov), sep=''), stringsAsFactors=F)   
    rownames(cov_names_replace) = cov_names_replace[,'newname']
    colnames(data_cov) <- cov_names_replace[,'newname']
    cov_names = colnames(data_cov)

    if(sum(row.names(data_input) %in% row.names(data_cov)) != nrow(data_input) )
    {
        message("[Error]  输入的cov文件没有包含所有的input样本\n")
        q()  
    }
}

# 命名汇总
intercept = data.frame(from='(Intercept)', newname='(Intercept)')
rownames(intercept) = '(Intercept)'
names_replace = rbind(genes_replace, cov_names_replace, intercept)



# 开始分析
vars = genes  # 所有变量

data_tmp = data_input[, c('y', genes)]
if(! is.null(cov))
{
    data_tmp = cbind(data_tmp, data_cov[row.names(data_tmp), ])
    vars = c(vars, cov_names)
}

data_tmp = data_tmp[complete.cases(data_tmp), ]  # 去掉缺失值

# 检查数据是否异常
nmiss_case    = sum(data_tmp$y == 1)
nmiss_control = sum(data_tmp$y == 0)
if(nmiss_case == 0 | nmiss_control == 0)
{   
    message("[Error]  data error, skip. case = ", nmiss_case, ' . control = ', nmiss_control)
    next
} 

# 逻辑回归
message("build model")
f <- as.formula(paste('y ~', paste(vars, collapse = ' + ') ) )
glm_fit <- glm(f, data = data_tmp, family = binomial, control=list(maxit=maxit))
fit_summary <- summary(glm_fit)
glm_result  <- fit_summary$coefficients
predict_value = predict(glm_fit, type='response')
predict_value = predict_value[rownames(data_tmp)]  # 保证顺序一致

result_predict = data.frame(sample = rownames(data_tmp), group= data_tmp$y, predicted = predict_value)  # 预测值用于输出

# ROC参数计算
message("ROC plot")
modelroc = roc(data_tmp$y, predict_value, transpose = FALSE, quiet = TRUE)
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
pdf_file = paste0(prefix, '.ROC.pdf')
pdf(pdf_file)
par(family="Times")  # times new roman 字体
prob <- predict(glm_fit, newdata = data_tmp, type = "response")
pred <- prediction(prob, data_tmp$y)
perf = performance(pred, "tpr", "fpr")
plot(perf@x.values[[1]], perf@y.values[[1]], lwd = 3, col = "#2171B5", type = "l", bty = "n", axes = FALSE, ann = FALSE)
abline(a = 0, b = 1, lty = 5, lwd = 3, col = "#D94801" )
par(tck = 0.03, col.axis = "#084594",lwd = 1, font.axis = 2, cex.axis = 1 )
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), col =  "#8C2D04", lwd = 2)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), col =  "#8C2D04", lwd = 2)
minor.tick(nx = 4, ny = 4, tick.ratio = 0.01)
box(which = "plot", col = "#8C2D04", lwd = 3)
title(xlab = 'False positive rate', ylab = 'True positive rate',  main = "multi factor predict", col.lab = "#084594", col.main = "#084594")
text (0.4, 0.3, c(paste("AUC:", AUC)), adj = 0, font = 2)
text (0.4, 0.2, c(paste("95%CI:", L95, "-", U95)), adj = 0, font = 2)


# 结果输出
result_model = data.frame(Var = names_replace[rownames(glm_result), 'from'], glm_result, check.names = F)
result_best  = matrix(c(nmiss_case, nmiss_control, AUC, AUC_CI95, thr, spe, sen, acc), 1, 8)
colnames(result_best) = c('NMISS case', 'NMISS control', "AUC", "AUC_CI95", 'threshold', "specificity", "sensitivity", "accuracy")

file_model    = paste0(prefix, '.model.txt')
file_best     = paste0(prefix, '.best.txt')
file_predict  = paste0(prefix, '.predict.txt')
 
write.table(result_model, file_model, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )
write.table(result_best, file_best, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )
write.table(result_predict, file_predict, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )


message("ROC:              ", pdf_file)
message("logistic model:   ", file_model)
message("logistic best:    ", file_best)
message("logistic predict: ", file_predict)

