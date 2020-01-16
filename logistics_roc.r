#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: logistic_roc.r  -i <file> -o <dir> [ --no_single_result --no_summary_result --rlib <dir>]

Options:
    -i, --input <file>   输入文件，第一列为样本名，第二列为样本分组（case=1，control=0），往后为基因列
    -o, --output <dir>   输出文件路径
    --no_single_result   不做每个基因独立建模分析 
    --no_summary_result  不做所有基因汇总建模分析
    --rlib <dir>         R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，逻辑回归、ROC\n支持单因素与多因素逻辑回归\n')
input             <- opts$input
output_dir        <- opts$output
no_single_result  <- opts$no_single_result
no_summary_result <- opts$no_summary_result
rlib              <- opts$rlib

if(!is.null(rlib)) .libPaths(rlib)

 
library(ROCR)
library(pROC)
library(Hmisc)

set.seed(91)


# 输出逻辑拟合结果

 
roc_result <- function(data, prefix, genes, main_title){

    # 逻辑拟合结果
	
    f <- as.formula(paste('y ~',paste(genes, collapse = ' + ') ) )
    glm_fit <- glm(f, data = data, family = binomial)
    fit_summary <- summary(glm_fit)
    glm_result  <- fit_summary$coefficients
    if(main_title == 'SUMMARY')
    {   
        glm_result = cbind(variables = rownames(glm_result), glm_result)
        write.table(glm_result, paste0(prefix, '.logistic.txt') , sep = "\t", quote = F, row.names = F , col.names = T)
        
    }

    # ROC参数计算

    modelroc = roc(data$y, predict(glm_fit, type='response'), transpose = FALSE)
    best = coords(modelroc, "best", ret = c("threshold", "specificity", "sensitivity", "accuracy"))
    best = as.matrix(best)
    AUC  = round(modelroc$auc[[1]], 4)
    spe  = round(best["specificity", 1], 4)
    sen  = round(best["sensitivity", 1], 4)
    acc  = round(best["accuracy", 1], 4)
    CI   = signif(ci.auc(modelroc, conf.level=0.95, method=c("delong", "bootstrap")), 4)
    L95 = CI[1]
    U95 = CI[3]
	AUC_CI95 = paste(L95, "-", U95)
    NMISS = nrow(data)
	

   
    # 绘图
    pdf(paste0(prefix, '.logistic.ROC.pdf'), family = "GB1")
    prob <- predict(glm_fit, newdata = data, type = "response")
    pred <- prediction(prob, data$y)
    perf = performance(pred, "tpr", "fpr")
    plot(perf@x.values[[1]], perf@y.values[[1]], lwd = 3, col = "#2171B5", type = "l", bty = "n", axes = FALSE, ann = FALSE)
    abline(a = 0, b = 1, lty = 5, lwd = 3, col = "#D94801" )
    par(tck = 0.03, col.axis = "#084594",lwd = 1, font.axis = 2, cex.axis = 1 )
    axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), col =  "#8C2D04", lwd = 2)
    axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), col =  "#8C2D04", lwd = 2)
    minor.tick(nx = 4, ny = 4, tick.ratio = 0.01)
    box(which = "plot", col = "#8C2D04", lwd = 3)
    title(xlab = 'False positive rate', ylab = 'True positive rate',  main = main_title, col.lab = "#084594", col.main = "#084594")
    text (0.4, 0.3, c(paste("AUC:", AUC)), adj = 0, font = 2)
    text (0.4, 0.2, c(paste("95%CI:", L95, "-", U95)), adj = 0, font = 2)
    dev.off()

    # 逻辑拟合结果与ROC参数输出

    result = c(genes, NMISS, glm_result[2, ], AUC, AUC_CI95, spe, sen, acc)
    return(result)
	
}


data = read.table(input, head = T, row.names = 1, check.names = F)
colnames(data)[1] <- 'y'
genes <- colnames(data)[2:ncol(data)]

# 基因汇总
if (no_summary_result == FALSE){
    message("summary process")

    prefix = paste0(output_dir, "/summary")
    data_tmp = data[complete.cases(data), ]
    tmp <- roc_result(data_tmp, prefix, genes, "SUMMARY") #因为返回的逻辑结果只有一个基因，所以这一项返回值不用
}

# 单个基因
if (no_single_result == FALSE){

    message("single process")
    result <- matrix(ncol = 11, nrow = length(genes))
    colnames(result) = c('Gene', 'NMISS', "Estimate", "Std. Error", "z value", "Pr(>|z|)", "AUC", "AUC_CI95", "specificity", "sensitivity", "accuracy")

    for (col in 1:length(genes))
    {   
        gene = genes[col]
		prefix = paste0(output_dir, '/single.', gene)
        data_tmp = data[, c('y', gene)]
        data_tmp = data[complete.cases(data), ]
        result[col,] <- roc_result(data, prefix, gene, gene)
    }
    write.table(result, paste(output_dir, 'single.logistic.txt', sep = '/'), quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )
}
