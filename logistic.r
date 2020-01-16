.libPaths("/home/ganb/soft/R-3.4.0/lib64/R/library/")

library(docopt)

"Usage: lasso.r  [options]  INPUT OUTPUT 
Options:
   -s --srgan=srgan    the organism, the value in <hsa, mmu, rno> [default: hsa]

Arguments:
   INPUT       the name of input files
   OUTPUT      the output dir name prefix" -> doc

opts   <- docopt(doc)
input  <- opts$INPUT
prefix <- opts$OUTPUT 


# 输入文件数据格式，第一行表头，至少三列，第一列样本名，第二列分组（取值0、1），第三列往后是每一个基因,

library(ROCR)
library(pROC)
set.seed(91)

table = read.table(input, head = T, row.names = 1, check.names = F, sep = "\t")
colnames(table)[1] <- 'y'
table_colnames <- colnames(table)
genes <- table_colnames[2:length(table_colnames)]

# 拟合
f <- as.formula(paste('y ~',paste(genes , collapse = ' + ') ) )
glm_fit <- glm(f, data = table, family = binomial)
fit_summary <- summary(glm_fit)
glm_result  <- fit_summary$coefficients
glm_result <- cbind(type = rownames(glm_result), glm_result)
write.table(glm_result, paste(prefix, '.logistic.txt', sep = ''), quote = FALSE, row.names = FALSE, sep = '\t')

# ROC
pdf(paste(prefix, '.logistic.ROC.pdf', sep = ''), family="GB1")
prob <- predict(glm_fit, newdata = table, type = "response")
modelroc = roc(table$y, predict(glm_fit, type='response'))
best = coords(modelroc, "best", ret = c("threshold", "specificity", "sensitivity", "accuracy"))
best = as.matrix(best)
AUC  = round(modelroc$auc[[1]], 4)
spe  = round(best["specificity",1], 4)
sen  = round(best["sensitivity",1], 4)
acc  = round(best["accuracy",1], 4)

pred <- prediction(prob, table$y)
plot(performance(pred, "tpr", "fpr") , lwd = 3, main = "ROC Curves")
abline(a = 0, b = 1, lty = 2, lwd = 3, col = "black")
legend("bottomright", legend = c(paste("AUC =", AUC), paste("specificity =", spe), paste("sensitivity =", sen), paste("accuracy =", acc) ) ) 
dev.off()

# predict
predict_result = data.frame(sample=rownames(table), class=table$y, predict=prob)
write.table(predict_result, paste(prefix, '.logistic.predict.txt', sep = ''), quote = FALSE, row.names = FALSE, sep = '\t')

