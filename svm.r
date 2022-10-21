#!/home/genesky/software/r/4.0.3/bin/Rscript
library("argparse")

# 参数解析
parse_args <- function(){
    
    parser <- ArgumentParser(description='基于SVM，对样本做二分类', formatter_class="argparse.RawTextHelpFormatter", epilog="\\n\\n结果生成4-7个文件： \\n\\nprefix.train_set_classify.best_cutoff_statistic.txt      训练样本集合的ROC最佳阈值的分类结果 - 分类敏感性、特异性、准确性。\\nprefix.train_set_classify.default0_cutoff_statistic.txt  训练样本集合的默认SVM阈值0的分类结果 - 分类敏感性、特异性、准确性。\\nprefix.train_set_classify.txt                            训练样本集合的默认SVM阈值0的分类结果。\\nprefix.train_set.ROC.pdf                                 训练样本集合的ROC曲线。\\nprefix.test_set_classify.default0_cutoff_statistic.txt   测试样本集合的默认SVM阈值0的分类结果 - 分类敏感性、特异性、准确性。\\nprefix.test_set_classify.txt                             测试样本集合的默认SVM阈值0的分类结果。\\nprefix.test_set.ROC.pdf                                  测试样本集合的ROC曲线。\\nprefix.svm.model.scale.txt/prefix.svm.model.txt          svm模型，客户不需要。\\n")

    parser$add_argument('-i', '--input', type='character', metavar='file', required=TRUE,
        help = "输入文件，第一列为样本名，第二列为样本分组（取值必须是 1或者0，其中 case=1，control=0），往后为基因列。")

    parser$add_argument('--test', type='character', metavar='file', required=FALSE,
        help = "测试集输入文件,格式与input一样。使用 input数据训练出来的模型，对测试集数据做分类检验。如果不提供，则不做")

    parser$add_argument('-g', '--gene', type='character', metavar='gene_list', required=FALSE, default='所有基因',
        help = "输入要纳入SVM分类器分析的基因名称，多个基因用逗号分隔， [default: %(default)s]")

    parser$add_argument( '-o', '--output_dir', type='character', metavar='dir', default='./',
        help = "结果输出目录, 默认是当前路径。流程会自动创建。[default: %(default)s]。")

    parser$add_argument('-p', '--prefix', type='character', metavar='name', default="svm",
        help = "输出的结果文件前缀。  [default: %(default)s] ")

    parser$add_argument('-k', '--kernel', type='character', metavar='kernel', default="radial",
        help = "设置svm核函数。没有最好的选择，需要全都试一遍才行,根据ROC曲线、分类准确率信息，自己选择。 可选参数有： linear / radial / polynomial / sigmoid [default: %(default)s] ")

    parser$add_argument('-t', '--tune', action="store_true", default=FALSE,
        help = "是否开启svm最优参数组合选取模式。非常耗时。 [default: %(default)s]")

    args <- parser$parse_args()
    # 命令行测试参数
    # args <- parser$parse_args(c('-p', '1', '-g', '2', '-o', '3'))

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
library(e1071)
library(pROC)
library(ROCR)
library(Hmisc)
set.seed(91)

do_roc <- function(raw_class, predict_value, cutoff, pdf_file, txt_file, if_plot){
    # raw_class: 样本的原始分类向量，用0、1 表示， 0 表示 control, 1表示对照
    # predict_value: 样本基于模型的预测值，是连续型变量
    #                要保证 raw_class 和 predict_value 的样本 一一对应
    # cutoff: 对 predict_value 选取的分类阈值。可选值： best 或 数值
    #         best: 表示自动选择最佳分类阈值
    #         数值：根据自定义数值，大于阈值的分类为 1、小于阈值的分类为0
    # pdf_file: roc绘图输出pdf文件路径
    # txt_file: roc结果文件
    # if_plot: 取值 TRUE、FALSE。 TRUE 表示绘图。FALSE 表示不绘图，同时txt结果也不给AUC信息
    modelroc = roc(raw_class, predict_value, transpose = FALSE, quiet = TRUE)
    # cutoff  系数阈值。
    choice_model = coords(modelroc, cutoff, ret = c("threshold", "specificity", "sensitivity", "accuracy"))
    choice_model = as.matrix(choice_model)
    AUC  = round(modelroc$auc[[1]], 4)
    thr  = round(choice_model[1, "threshold"], 4)
    spe  = round(choice_model[1, "specificity"], 4)
    sen  = round(choice_model[1, "sensitivity"], 4)
    acc  = round(choice_model[1, "accuracy"], 4)
    CI   = signif(ci.auc(modelroc, conf.level=0.95, method=c("delong", "bootstrap")), 4)
    L95 = CI[1]
    U95 = CI[3]
    AUC_CI95 = paste(L95, "-", U95)
    nmiss_case = sum(raw_class == 1)
    nmiss_control = sum(raw_class == 0)

    if(! if_plot){
        result = matrix(c(nmiss_case, nmiss_control, thr, spe, sen, acc), 1, 6)
        colnames(result) = c('NMISS case', 'NMISS control', 'threshold', "specificity", "sensitivity", "accuracy")
        write.table(result, txt_file, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )
    }else{
        result = matrix(c(nmiss_case, nmiss_control, AUC, AUC_CI95, thr, spe, sen, acc), 1, 8)
        colnames(result) = c('NMISS case', 'NMISS control', "AUC", "AUC_CI95", 'threshold', "specificity", "sensitivity", "accuracy")
        write.table(result, txt_file, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )
        # 绘图
        pdf(pdf_file)
        par(family="Times")  # times new roman 字体
        pred <- prediction(predict_value, raw_class)
        perf = performance(pred, "tpr", "fpr")
        plot(perf@x.values[[1]], perf@y.values[[1]], lwd = 3, col = "#2171B5", type = "l", bty = "n", axes = FALSE, ann = FALSE)
        abline(a = 0, b = 1, lty = 5, lwd = 3, col = "#D94801" )
        par(tck = 0.03, col.axis = "#084594",lwd = 1, font.axis = 2, cex.axis = 1 )
        axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), col =  "#8C2D04", lwd = 2)
        axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), col =  "#8C2D04", lwd = 2)
        minor.tick(nx = 4, ny = 4, tick.ratio = 0.01)
        box(which = "plot", col = "#8C2D04", lwd = 3)
        title(xlab = 'False positive rate', ylab = 'True positive rate',  main = "ROC", col.lab = "#084594", col.main = "#084594")
        text (0.4, 0.3, c(paste("AUC:", AUC)), adj = 0, font = 2)
        text (0.4, 0.2, c(paste("95%CI:", L95, "-", U95)), adj = 0, font = 2)
        dev.off()
    }
}

message("1. 输入数据预处理")

data_input = read.table(args$input, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "")
colnames(data_input)[1] = 'group'

# 挑选基因
genes = colnames(data_input)[2:ncol(data_input)]

if(args$gene != "所有基因"){
    genes = unlist(strsplit(args$gene, ","))
    lost_gene = genes[! genes %in% colnames(data_input)]
    if(length(lost_gene) > 0){
        stop("参数 gene 异常，基因在文件中丢失：", paste(lost_gene, collapse=','))
    }
    
}
data_input_select = data_input[, c('group', genes)]

# 替换基因名，防止特殊字符

genes_replace = data.frame(from=colnames(data_input_select)[2:ncol(data_input_select)], newname=paste('gene', 2:ncol(data_input_select), sep=''), stringsAsFactors=F)
rownames(genes_replace) = genes_replace[, 'newname']
colnames(data_input_select) <- c('group', genes_replace[,'newname'])


# 去掉存在缺失值的样本
data_clean = data_input_select[complete.cases(data_input_select), ]  

# 检查数据是否异常
nmiss_case    = sum(data_clean$group == 1)
nmiss_control = sum(data_clean$group == 0)
if(nmiss_case == 0 | nmiss_control == 0)
{   
    stop("[Error]  数据异常, case/control 样本缺失 不能分析. case = ", nmiss_case, ' . control = ', nmiss_control)
} 
# 鉴于svm的特殊要求，需要把分组设为factor
data_clean_bak = data_clean[,1,drop=F]  # 备份group信息，后续ROC要用到
data_clean$group = as.character(data_clean$group)
data_clean$group[data_clean$group == '1'] = 'case'
data_clean$group[data_clean$group == '0'] = 'control'
data_clean$group = factor(data_clean$group, levels=c('control', 'case'))
samples = rownames(data_clean)



if(args$tune){
    message("2.0 svm 选取最优参数组合（比较耗时，请耐心等待）")
    obj_tune_svm <- tune.svm(group ~ ., data=data_clean, kernel=args$kernel, gamma = 10 ^ (-6:1), degree = c(2:10), cost = c(1:100))
}
message("2. svm 分类")
if(args$tune){
    svm_model <- svm(group ~ ., data=data_clean, kernel=args$kernel, gamma = obj_tune_svm$best.parameters$gamma, degree = obj_tune_svm$best.parameters$degree, cost = obj_tune_svm$best.parameters$cost)  # 训练模型
}else{
    svm_model <- svm(group ~ ., data=data_clean, kernel=args$kernel)  # 训练模型
}
# 模型输出
write.svm(svm_model, svm.file = paste0(args$output_dir, "/", args$prefix, ".svm.model.txt"),  scale.file = paste0(args$output_dir, "/", args$prefix, ".svm.model.scale.txt"))

# 预测结果输出
svm_dist <- svm_model$decision.values  # 样本距离决策平面的距离。在二分类模型中，默认情况下，以0作为分界点。
svm.pred <- predict(svm_model, data_clean)  # 预测分类

write.table( data.frame(sample=samples, group=data_clean$group, svm_decision.values=svm_dist[samples,1], svm_predict=svm.pred[samples]),
             paste0(args$output_dir, "/", args$prefix, ".train_set_classify.txt"), quote=F, row.names=F, col.names=T, sep='\t')




message("3. ROC (基于样本到超平面的距离)")
# 最佳阈值分析
do_roc(data_clean_bak$group, svm_dist[samples,1], 'best', paste0(args$output_dir, "/", args$prefix, '.train_set.ROC.pdf'), paste0(args$output_dir, "/", args$prefix, ".train_set_classify.best_cutoff_statistic.txt"), TRUE)
# 默认阈值 0 分析
do_roc(data_clean_bak$group, svm_dist[samples,1], 0, '', paste0(args$output_dir, "/", args$prefix, ".train_set_classify.default0_cutoff_statistic.txt"), FALSE)


if(! is.null(args$test)){
    message("4. 基于训练模型，对测试集做分类")

    data_input_test = read.table(args$test, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "")
    colnames(data_input_test)[1] = 'group'
    # 测试集基因是否缺失
    lost_gene = genes[! genes %in% colnames(data_input_test)]
    if(length(lost_gene) > 0){
        stop("测试集缺失基因：", paste(lost_gene, collapse=','))
    }
    # 提取基因，并重命名
    data_input_test_select = data_input_test[, c('group', genes)]
    colnames(data_input_test_select) <- c('group', genes_replace[,'newname'])

    # 删除存在缺失值的样本
    data_clean_test = data_input_test_select[complete.cases(data_input_test_select), ]
    nmiss_case    = sum(data_clean_test$group == 1)
    nmiss_control = sum(data_clean_test$group == 0)

    # 分组信息重命名
    data_clean_test_bak = data_clean_test[,1,drop=F]  # 备份group信息，后续ROC要用到
    data_clean_test$group = as.character(data_clean_test$group)
    data_clean_test$group[data_clean_test$group == '1'] = 'case'
    data_clean_test$group[data_clean_test$group == '0'] = 'control'
    data_clean_test$group = factor(data_clean_test$group, levels=c('control', 'case'))
    samples_test = rownames(data_clean_test)

    # 预测
    svm.pred_test <- predict(svm_model, data_clean_test, decision.values = TRUE)  # 预测分类
    svm.pred_classify <- as.data.frame(svm.pred_test)
    svm.pred_decision_values <- attr(svm.pred_test, "decision.values")

    # 预测结果输出
    write.table( data.frame(sample=samples_test, group=data_clean_test$group, svm_decision.values=svm.pred_decision_values[samples_test,1], svm_predict=svm.pred_classify[samples_test,1]),
             paste0(args$output_dir, "/", args$prefix, ".test_set_classify.txt"), quote=F, row.names=F, col.names=T, sep='\t')
    # ROC
    do_roc(data_clean_test_bak$group, svm.pred_decision_values[samples_test,1], 0, paste0(args$output_dir, "/", args$prefix, '.test_set.ROC.pdf'), paste0(args$output_dir, "/", args$prefix, ".test_set_classify.default0_cutoff_statistic.txt"), TRUE)
}
