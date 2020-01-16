.libPaths("/home/ganb/soft/R-3.4.0/lib64/R/library/")

library(docopt)

"Usage: lasso.r  [options]  INPUT1 INPUT2 NAME1 NAME2 OUTPUT 
Options:
   -s --srgan=srgan    the organism, the value in <hsa, mmu, rno> [default: hsa]

Arguments:
   INPUT1       input file1
   INPUT2       input file2
   NAME1        the name of input file1
   NAME2        the name of input file2
   OUTPUT       the output dir name" -> doc

opts    <- docopt(doc)
input1  <- opts$INPUT1
input2  <- opts$INPUT2
name1   <- opts$NAME1
name2   <- opts$NAME2
prefix  <- opts$OUTPUT 

# 输入文件数据格式，第一行表头，至少三列，第一列样本名，第二列分组（取值0、1），第三列往后是每一个基因,

library(ROCR)
library(pROC)
set.seed(91)

# 数据 1
table = read.table(input1, head = T, row.names = 1, check.names = F)
colnames(table)[1] <- 'y'
table_colnames <- colnames(table)
genes <- table_colnames[1:length(table_colnames)]

# 拟合 1
f <- as.formula(paste('y ~',paste(genes , collapse = ' + ') ) )
glm_fit <- glm(f, data = table, family = binomial)
fit_summary <- summary(glm_fit)
glm_result  <- fit_summary$coefficients
write.table(glm_result, paste(prefix, '.', name1, '.logistic.txt', sep = ''), quote = FALSE, row.names = TRUE, sep = '\t')


# 数据 2
table2 = read.table(input2, head = T, row.names = 1, check.names = F)
colnames(table2)[1] <- 'y'
table2_colnames <- colnames(table2)
genes2 <- table2_colnames[1:length(table2_colnames)]

# 拟合 2
f <- as.formula(paste('y ~',paste(genes2 , collapse = ' + ') ) )
glm_fit2 <- glm(f, data = table2, family = binomial)
fit_summary2 <- summary(glm_fit2)
glm_result2  <- fit_summary2$coefficients
write.table(glm_result2, paste(prefix, '.', name2, '.logistic.txt', sep = ''), quote = FALSE, row.names = TRUE, sep = '\t')


# ROC 1
prob <- predict(glm_fit, newdata = table, type = "response")
modelroc = roc(table$y, predict(glm_fit, type='response'))
best = coords(modelroc, "best", ret = c("threshold", "specificity", "sensitivity", "accuracy"))
best = as.matrix(best)
AUC  = round(modelroc$auc[[1]], 4)
spe  = round(best["specificity",1], 4)
sen  = round(best["sensitivity",1], 4)
acc  = round(best["accuracy",1], 4)

# ROC 2
prob2 <- predict(glm_fit2, newdata = table2, type = "response")
modelroc2 = roc(table2$y, predict(glm_fit2, type='response'))
best2 = coords(modelroc2, "best", ret = c("threshold", "specificity", "sensitivity", "accuracy"))
best2 = as.matrix(best2)
AUC2  = round(modelroc2$auc[[1]], 4)
spe2  = round(best2["specificity",1], 4)
sen2  = round(best2["sensitivity",1], 4)
acc2  = round(best2["accuracy",1], 4)

# 绘图
pdf(paste(prefix, '.logistic.ROC.pdf', sep = ''), family="GB1")
pred <- prediction(prob, table$y)
pred2 <- prediction(prob2, table$y)
roc_test = roc.test(modelroc, modelroc2)
roc_pvalue = round(roc_test$p.value, 4)

plot(performance(pred, "tpr", "fpr") , lwd = 3, main = "ROC Curves", col = 'red')
plot(performance(pred2, "tpr", "fpr") , lwd = 3, main = "ROC Curves", col = 'blue', add = T)
abline(a = 0, b = 1, lty = 2, lwd = 3, col = "black")

legend_info = c(paste("pvalue(ROC compare) = ", roc_pvalue, sep = ''))
# legend_info = c(paste("ROC pvalue(", name1, " vs ", name2, ") = ", roc_pvalue, sep = ''))
legend_info = c(legend_info, paste("AUC(", name1, ") = ", AUC, sep = ''), paste("specificity(", name1, ") = ", spe, sep = ''), paste("sensitivity(", name1, ") = ", sen, sep = ''), paste("accuracy(", name1, ") = ", acc, sep = '') )
legend_info = c(legend_info, paste("AUC(", name2, ") = ", AUC2, sep = ''), paste("specificity(", name2, ") = ", spe2, sep = ''), paste("sensitivity(", name2, ") = ", sen2, sep = ''), paste("accuracy(", name2, ") = ", acc2, sep = '') )
legend("bottomright", legend = legend_info, text.col = c('black', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue') )  
dev.off()


