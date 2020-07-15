#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: logistics_single_factor.r  -i <file> -p <prefix> [ --rlib <dir>]

Options:
    -i, --input <file>      输入文件，第一列为样本名，第二列为样本分组（case=1，control=0），往后为基因列。该软件对每一个基因单独做逻辑回归分析。
    -p, --prefix <prefix>   输出文件前缀， 生成: (1) ROC曲线图 prefix.ROC.pdf  (2) 逻辑回归模型结果 prefix.model.txt (3) 基于逻辑回归模型对输入样本做预测 prefix.predict.txt
    --rlib <dir>            R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，单因素逻辑回归、ROC \n')
input         <- opts$input
prefix        <- opts$prefix
rlib          <- opts$rlib

if(!is.null(rlib)) .libPaths(rlib)

 
library(ROCR)
library(pROC)
library(Hmisc)

set.seed(91)


# 读入
data = read.table(input, head = T, row.names = 1, check.names = F, sep="\t")
colnames(data)[1] <- 'y'
genes <- colnames(data)[2:ncol(data)]



# 保留预测值
result_pred <- matrix(ncol = length(genes) + 2, nrow = nrow(data))  
rownames(result_pred) = rownames(data)
colnames(result_pred) = c('sample', 'group', genes)
result_pred[, 'sample'] = rownames(data)
result_pred[, 'group']  = data[,1]

# 保留模型值
result <- matrix(ncol = 17, nrow = length(genes))
colnames(result) = c('Gene', 'NMISS case', 'NMISS control', "Estimate(argument)", "Std. Error(argument)", "z value(argument)", "Pr(>|z|)(argument)", "Estimate(Intercept)", "Std. Error(Intercept)", "z value(Intercept)", "Pr(>|z|)(Intercept)", "AUC", "AUC_CI95", 'threshold', "specificity", "sensitivity", "accuracy")

pdf_file = paste0(prefix, '.ROC.pdf')
pdf(pdf_file, family = "GB1")

for( col in 1:length(genes))
{   
    gene = genes[col]

    message("process : ", gene)

    data_tmp = data[, c('y', gene)]
    data_tmp = data_tmp[complete.cases(data_tmp), ]  # 去掉缺失值
    colnames(data_tmp)[2] = 'gene'  #临时替换列名，防止存在“-”符号，导致逻辑回归出错

    # 检查数据是否异常
    nmiss_case    = sum(data_tmp$y == 1)
    nmiss_control = sum(data_tmp$y == 0)
    if(nmiss_case == 0 | nmiss_control == 0)
    {   
        message("    data error, skip. case = ", nmiss_case, ' . control = ', nmiss_control)
        result[col, 1:3] = c(gene, nmiss_case, nmiss_control) 
        next
    }   

    # 逻辑回归
    glm_fit <- glm(y ~ gene, data = data_tmp, family = binomial)
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
    result[col, ] <- c(gene, nmiss_case, nmiss_control, glm_result['gene', ], glm_result['(Intercept)', ], AUC, AUC_CI95, thr, spe, sen, acc)
    result_pred[names(prob), gene] = prob
}
dev.off()


# 结果输出
model_file = paste0(prefix, '.model.txt')
pred_file = paste0(prefix, '.predict.txt')

write.table(result, model_file, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )
write.table(result_pred, pred_file, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )


message("ROC:              ", pdf_file)
message("logistic model:   ", model_file)
message("logistic predict: ", pred_file)

