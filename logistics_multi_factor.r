#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: logistics_multi_factor.r  -i <file> -g <string> -p <prefix> [ --rlib <dir>]

Options:
    -i, --input <file>      输入文件，第一列为样本名，第二列为样本分组（case=1，control=0），往后为基因列。
    -g, --gene <string>     输入要纳入多元逻辑回归分析的基因名称，多个基因用逗号分隔
    -p, --prefix <prefix>   输出文件前缀， 生成: (1) ROC曲线图 prefix.ROC.pdf  (2) 逻辑回归模型结果 prefix.model.txt (3) 基于逻辑回归模型对输入样本做预测 prefix.predict.txt  (4) 最佳阈值结果 prefix.best.txt
    --rlib <dir>            R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，多因素逻辑回归、ROC \n')
input         <- opts$input
prefix        <- opts$prefix
gene          <- opts$gene
rlib          <- opts$rlib
genes <- unlist(strsplit(gene, ","))

if(!is.null(rlib)) .libPaths(rlib)
 
 
library(ROCR)
library(pROC)
library(Hmisc)

set.seed(91)


# 读入
data = read.table(input, head = T, row.names = 1, check.names = F, sep="\t")
colnames(data)[1] <- 'y'

if(sum(genes %in% colnames(data)) != length(genes) )
{
    message("[Error] 输入的基因列表有问题，部分基因在输入文件中缺失")
    q()
}


# 开始分析
data_tmp = data[, c('y', genes)]
data_tmp = data_tmp[complete.cases(data_tmp), ]  # 去掉缺失值

# 检查数据是否异常
nmiss_case    = sum(data_tmp$y == 1)
nmiss_control = sum(data_tmp$y == 0)
if(nmiss_case == 0 | nmiss_control == 0)
{   
    message("    data error, skip. case = ", nmiss_case, ' . control = ', nmiss_control)
    next
} 

# 逻辑回归
message("build model")
f <- as.formula(paste('y ~', paste(genes, collapse = ' + ') ) )
glm_fit <- glm(f, data = data_tmp, family = binomial)
fit_summary <- summary(glm_fit)
glm_result  <- fit_summary$coefficients


# ROC参数计算
message("ROC plot")
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
result_pred  = data.frame(sample = rownames(data_tmp), group = data_tmp$y, prob = prob[rownames(data_tmp)], check.names = F) 
result_model = data.frame(gene = rownames(glm_result), glm_result, check.names = F)
result_best  = matrix(c(nmiss_case, nmiss_control, AUC, AUC_CI95, thr, spe, sen, acc), 1, 8)
colnames(result_best) = c('NMISS case', 'NMISS control', "AUC", "AUC_CI95", 'threshold', "specificity", "sensitivity", "accuracy")

file_predict  = paste0(prefix, '.predict.txt')
file_model    = paste0(prefix, '.model.txt')
file_best     = paste0(prefix, '.best.txt')
 
write.table(result_pred, file_predict, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )
write.table(result_model, file_model, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )
write.table(result_best, file_best, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )


message("ROC:              ", pdf_file)
message("logistic model:   ", file_model)
message("logistic predict: ", file_predict)
message("logistic best:    ", file_best)

