.libPaths("/home/ganb/soft/R-3.4.0/lib64/R/library/")

library(docopt)

"Usage: lasso.r  [options]  INPUT OUTPUTDIR LAMBDA
Options:
   -s --srgan=srgan    the organism, the value in <hsa, mmu, rno> [default: hsa]

Arguments:
   INPUT       the name of input files
   OUTPUTDIR   the output dir name
   LAMBDA      lambda type min/1se" -> doc

opts   <- docopt(doc)
input  <- opts$INPUT
prefix <- opts$OUTPUTDIR 
lambda <- opts$LAMBDA 

# 输入文件数据格式，第一行表头，至少三列，第一列样本名，第二列分组（取值0、1），第三列往后是每一个基因,

library(glmnet)
library(ROCR)
library(pROC)
set.seed(91)

table = read.table(input, head = T, row.names = 1, check.names = F)
table = data.matrix(table)

y = table[, 1]
x = table[, 2:dim(table)[2]]
y[y == 1] = 2
y[y == 0] = 1
genes = colnames(x)
geneRename = paste0('Gene', seq(1,length(genes))) # 基因重命名，因为名字有特殊符号，影响后面分析
colnames(x) = geneRename
colnames(table) = c('y', geneRename)

# 系数
cvfit  = cv.glmnet(x, y, family = "binomial", type.measure = "class")
lambda_value = ifelse(lambda == 'min', cvfit$lambda.min, cvfit$lambda.1se)
lambda_file  = ifelse(lambda == 'min', paste(prefix, '.coef.lambda.min.txt', sep = ''), paste(prefix, '.coef.lambda.1se.txt', sep = '')) 
names(lambda_value) = ifelse(lambda == 'min', 'lambda.min', 'lambda.1se')

coef   = coef(cvfit, s=c(lambda_value))
coef   = cbind(Gene = c('', genes), GeneRename = rownames(coef), Coef = coef[, 1])
write.table(coef, paste(prefix, '.coef.txt', sep = ''), quote = FALSE, row.names = FALSE, sep = '\t')
write.table(lambda_value, lambda_file, quote = FALSE, row.names = FALSE, sep = '\t')
    
# 绘图
pdf(paste(prefix, '.coef.pdf', sep = ''), family="GB1")
plot(cvfit)
dev.off()

# 绘图 2
pdf(paste(prefix, '.dif.pdf', sep = ''), family="GB1")
fit<-glmnet(x, y, family = "binomial")
plot(fit, xvar="lambda")
dev.off()

# 逻辑回归
coef = coef[-1, ] # 去掉第一行常数项
table = as.data.frame(table)
choose_gene       = coef[ coef[, 3] != "0", 1] # 要分析的基因
choose_geneRename = coef[ coef[, 3] != "0", 2] # 要分析的基因重命名

if(length(choose_gene) == 0)
{
    stop("no gene select, please change lambda\n")
} 
 
f <- as.formula(paste('y ~',paste(choose_geneRename , collapse = ' + ') ) )
glm_fit = glm(f, data = table, family = binomial)
fit_summary = summary(glm_fit)
glm_result  = fit_summary$coefficients
glm_result = cbind(Gene = c('(Intercept)',choose_gene), glm_result)
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
# coef_name = paste(file, 'fit.txt', sep = '.')
# write.table(fit_summary,coef_name)

