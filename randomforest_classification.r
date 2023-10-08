#!/home/genesky/software/r/4.0.3/bin/Rscript
library("argparse")

# 参数解析
parse_args <- function(){
    
    parser <- ArgumentParser(description='基于随机森林，对样本做二分类。基于 Gini 系数计算节点不纯度', formatter_class="argparse.RawTextHelpFormatter", epilog=" 结果文件包括： \\n prefix.train_set.importance.pdf 训练集下，随机森林中，每个特征的重要性绘图 \\n prefix.train_set.importance.txt 训练集下，随机森林中，每个特征的重要性。有4个统计量（MeanDecreaseGini最重要），值越大越重要。\\n prefix.train_set.predict.txt    训练集下，随机森林分类结果（默认参数）\\n prefix.train_set.ROC.pdf        训练集下，ROC 绘图\\n prefix.train_set.best.txt       训练集下，ROC 最佳阈值\\n prefix.test_set.predict.txt     测试集下，训练模型预测结果\\n prefix.test_set.statistic.txt   测试集下，训练模型预测结果统计\\n")

    parser$add_argument('-i', '--input', type='character', metavar='file', required=TRUE,
        help = "输入文件，第一列为样本名，第二列为样本分组（取值0或1，其中case=1，control=0），往后为基因列。")

    parser$add_argument('-g', '--gene', type='character', metavar='gene_list', required=FALSE, default='所有基因',
        help = "输入要纳入随机森林分析的基因名称，多个基因用逗号分隔， [default: %(default)s]")

    parser$add_argument('-c', '--cov', type='character', metavar='cov', required=FALSE,
        help = "协变量输入文件，用于矫正。第一列样本名，后面所有列都是协变量。样本顺序没有要求。默认不用输入。 务必保证cov中包含input的所有样本，cov文件中的所有列都用于矫正。列名不能包含空格、“-”等特殊字符")

    parser$add_argument('--test', type='character', metavar='file', required=FALSE,
        help = "测试集输入文件,格式与input一样。使用 input数据训练出来的模型，对测试集数据做分类检验，必须包含 input/cov 中的所有特征（第一列也必须是样本的分类信息）。如果不提供，则不做")

    parser$add_argument( '-o', '--output_dir', type='character', metavar='dir', default='./',
        help = "结果输出目录, 默认是当前路径。流程会自动创建。[default: %(default)s]。")

    parser$add_argument('-p', '--prefix', type='character', metavar='name', default="result",
        help = "输出的结果文件前缀。  [default: %(default)s] ")

    parser$add_argument('--ntree', type='integer', metavar='tree_count', default=100,
        help = "随机森林所包含的决策树数目 [default: %(default)s] ")

    args <- parser$parse_args()
    # 命令行测试参数
    # args <- parser$parse_args(c('--ntree','100', '-o','./','-p','result','-i', '/home/pub/output/research_and_customized_project/19B0925F-2/logistics_multi_factor/input/DMP.profile.txt', '-c', '/home/pub/output/research_and_customized_project/19B0925F-2/logistics_multi_factor/input/cov.txt', '-g', 'cg14827090,cg16287373,cg22113540,cg26654770,cg01201512,cg12746557,cg11070274,cg12102398,cg16154810,cg07234876,cg16474696,cg07909254,cg06715136,cg02095219,cg22325292,cg09863391,cg19944656,cg10786572,cg12690462,cg05840533,cg09417038,cg00036352,cg19097359'))
 
    return(args)
}

# 保存命令行信息，方便debug
OUT_COMMAND <- function(log_file){
    # 1. 时间
    time_now = format(Sys.time(), "%Y-%m-%d %H:%M:%OS")
    # 2. R路径
    rscript = gsub('/lib64/R/bin/exec/R$', '/bin/Rscript', commandArgs(trailingOnly = FALSE)[1], perl=TRUE)
    # 3. R脚本路径
    r_script = ""
    for(arg in commandArgs(trailingOnly = FALSE)){
        if(grepl('^--file=', arg, perl=TRUE)){
            r_script = gsub('^--file=', '', arg, perl=TRUE)
            break
        }
    }
    # 4. 参数组合
    parameter = c()
    for(arg in commandArgs(trailingOnly = TRUE)){
        if(grepl('^-', arg, perl=TRUE)){
            parameter = c(parameter, arg)
        }else{
            # 所有的值添加 单引号，防止有特殊字符导致传参异常
            parameter = c(parameter, paste0('\'', arg, '\''))
        }
    }

    # 5. 输出
    # 命令行
    command = paste(c(rscript, r_script, parameter), collapse = " ")
    # 测试用 argparse
    command_argparse = paste(commandArgs(trailingOnly = TRUE), collapse='\',\'')
    command_argparse = paste0('args <- parser$parse_args(c(\'', command_argparse, '\'))')
    
    write(paste0("[Command]: ", time_now, '\n', command, '\n', command_argparse, '\n\n'), file = log_file, append = TRUE)
}

args = parse_args()
if(!file.exists(args$output_dir)){
    dir.create(args$output_dir)
}

OUT_COMMAND(paste0(args$output_dir, '/command.info'))


# 加载R包
library(randomForest)
library(rfPermute)
library(pROC)
library(Hmisc)

set.seed(666)

# 读入
message("1. 读入input")
data_input = read.table(args$input, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "", fill=T)
colnames(data_input)[1] <- 'y'

# 要分析的基因
genes <- colnames(data_input)[2:ncol(data_input)]
if(args$gene != "所有基因"){
    genes = unlist(strsplit(args$gene, ","))
    lost_gene = genes[! genes %in% colnames(data_input)]
    if(length(lost_gene) > 0){
        stop("参数 gene 异常，基因在文件中丢失：", paste(lost_gene, collapse=','))
    }
}

# 可用数据
data_input_select = data_input[, c('y', genes)]
vars = genes 

# 读入cov
if(! is.null(args$cov))
{
    message("1.1 读入cov")
    data_cov = read.table(args$cov, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "", fill=T)
    if(sum(row.names(data_input_select) %in% row.names(data_cov)) != nrow(data_input) )
    {
        message("[Error]  输入的cov文件没有包含所有的input样本\n")
        q()  
    }
    vars = c(vars, colnames(data_cov))
    data_input_select = cbind(data_input_select, data_cov[row.names(data_input_select), ,drop=F])
}
# 去掉有缺失值的样本
message("2. 删除包含缺失值的样本")
data_clean = data_input_select[complete.cases(data_input_select), ]  # 去掉缺失值

# 检查数据是否异常
nmiss_case    = sum(data_clean$y == 1)
nmiss_control = sum(data_clean$y == 0)
if(nmiss_case == 0 | nmiss_control == 0)
{   
    message("[Error]  data error, skip. case = ", nmiss_case, ' . control = ', nmiss_control)
    q()
} 

# 建模
message("3. 随机森林建模")
f <- as.formula(paste('y ~', paste(vars, collapse = ' + ') ) )
data_clean_rf = data_clean
data_clean_rf$y = factor(as.character(data_clean_rf$y), levels=c('0', '1'))

modelRFP = rfPermute(f, data=data_clean_rf, ntree=args$ntree, proximity=TRUE, importance=TRUE)  # 带有置换检验、检查importance pvalue功能
modelRF = as.randomForest(modelRFP)  # 转换为 随机森林模型
# modelRF = randomForest(f, data=data_clean_rf, ntree=args$ntree, proximity=TRUE, importance=TRUE)
pred_value = predict(modelRF, data_clean, type='response')
pred_value = as.numeric(as.character(pred_value[rownames(data_clean)]))
result_predict = data.frame(Sample=rownames(data_clean), Class=data_clean$y, PredictedClass=pred_value)



if(! is.null(args$test)){
    message("3.1 基于训练模型，对测试集做分类")

    data_input_test = read.table(args$test, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "")
    colnames(data_input_test)[1] = 'y'
    # 测试集基因是否缺失
    lost_gene = vars[! vars %in% colnames(data_input_test)]
    if(length(lost_gene) > 0){
        stop("测试集缺失基因：", paste(lost_gene, collapse=','))
    }
    # 提取基因，并重命名
    data_input_test_select = data_input_test[, c('y', vars)]

    # 删除存在缺失值的样本
    data_clean_test = data_input_test_select[complete.cases(data_input_test_select), ]
    nmiss_case    = sum(data_clean_test$group == 1)
    nmiss_control = sum(data_clean_test$group == 0)

    data_clean_test_rf = data_clean_test
    data_clean_test_rf$y = factor(as.character(data_clean_test_rf$y), levels=c('0', '1'))

    pred_value_test = predict(modelRF, data_clean_test_rf, type='response')
    pred_value_test = as.numeric(as.character(pred_value_test[rownames(data_clean_test_rf)]))

    result = data.frame(Sample=rownames(data_clean_test), Class=data_clean_test$y, PredictClass=pred_value_test)

    # 分类效果
    modelroc = roc(data_clean_test$y, pred_value_test, transpose = FALSE, quiet = TRUE)
    # cutoff  系数阈值。
    choice_model = coords(modelroc, 0.5, ret = c("threshold", "specificity", "sensitivity", "accuracy"))
    choice_model = as.matrix(choice_model)
    spe  = round(choice_model[1, "specificity"], 4)
    sen  = round(choice_model[1, "sensitivity"], 4)
    acc  = round(choice_model[1, "accuracy"], 4)
    result_statistic = data.frame(specificity=spe, sensitivity=sen, accuracy=acc)

    # 预测结果输出 评价
    write.table( result, file.path(args$output_dir, paste0(args$prefix, '.test_set.predict.txt')), quote=F, row.names=F, col.names=T, sep='\t')
    write.table( result_statistic, file.path(args$output_dir, paste0(args$prefix, '.test_set.statistic.txt')), quote=F, row.names=F, col.names=T, sep='\t')
}

message("4. ROC")
modelroc = roc(data_clean$y, pred_value)
AUC  = round(modelroc$auc[[1]], 4)
CI   = signif(ci.auc(modelroc, conf.level=0.95, method=c("delong", "bootstrap")), 4)
L95 = CI[1]
U95 = CI[3]
AUC_CI95 = paste(L95, "-", U95)

best   = coords(modelroc, "best", ret=c("threshold","specificity","sensitivity","accuracy"))
best = as.matrix(best)
thr  = round(best[1, "threshold"], 4)
spe  = round(best[1, "specificity"], 4)
sen  = round(best[1, "sensitivity"], 4)
acc  = round(best[1, "accuracy"], 4)

# 绘图
pdf_file = file.path(args$output_dir, paste0(args$prefix, '.train_set.ROC.pdf'))  
pdf(pdf_file)
par(family="Times")
plot(1-modelroc$specificities,modelroc$sensitivities, lwd = 3, col = "#2171B5", type = "l", bty = "n", axes = FALSE, ann = FALSE)
abline(a = 0, b = 1, lty = 5, lwd = 3, col = "#D94801" )
par(tck = 0.03, col.axis = "#084594",lwd = 1, font.axis = 2, cex.axis = 1 )
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), col =  "#8C2D04", lwd = 2)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), col =  "#8C2D04", lwd = 2)
minor.tick(nx = 4, ny = 4, tick.ratio = 0.01)
box(which = "plot", col = "#8C2D04", lwd = 3)
title(xlab = 'False positive rate', ylab = 'True positive rate',  main = "randomForest predict", col.lab = "#084594", col.main = "#084594")
text (0.4, 0.3, c(paste("AUC:", AUC)), adj = 0, font = 2)
text (0.4, 0.2, c(paste("95%CI:", L95, "-", U95)), adj = 0, font = 2)
dev.off()



# 结果输出
message("5. 结果输出")

# ROC 结果文件
result_best  = matrix(c(nmiss_case, nmiss_control, AUC, AUC_CI95, thr, spe, sen, acc), 1, 8)
colnames(result_best) = c('NMISS case', 'NMISS control', "AUC", "AUC_CI95", 'threshold', "specificity", "sensitivity", "accuracy")
file_best     = file.path(args$output_dir, paste0(args$prefix, '.train_set.best.txt'))  
write.table(result_best, file_best, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )

# 预测结果文件
file_predict  = file.path(args$output_dir, paste0(args$prefix, '.train_set.predict.txt'))  
write.table(result_predict, file_predict, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )

# 特征重要性输出
file_importance = file.path(args$output_dir, paste0(args$prefix, '.train_set.importance.txt'))
importance_tmp = importance(modelRFP)
importance = data.frame(feature=rownames(importance_tmp), importance_tmp, check.names=F)
write.table(importance, file_importance, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )

# 特征重要性绘图
pdf_file_importance = file.path(args$output_dir, paste0(args$prefix, '.train_set.importance.pdf'))
height = 14
if(length(vars) > 28){
    height = length(vars) * 0.5
}

pdf(pdf_file_importance, width=14, height=height)
par(family="Times")
plotImportance(modelRFP)
dev.off()

message("ROC:              ", pdf_file)
message("randomForest best:    ", file_best)
message("randomForest predict:    ", file_predict)
message("randomForest importance:    ", file_importance)
message("randomForest importance:    ", pdf_file_importance)
