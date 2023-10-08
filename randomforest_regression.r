#!/home/genesky/software/r/4.0.3/bin/Rscript
library("argparse")

# 参数解析
parse_args <- function(){
    
    parser <- ArgumentParser(description='基于随机森林，对样本做回归分析。基于  residual sum of squares 计算节点不纯度', formatter_class="argparse.RawTextHelpFormatter", epilog=" 结果文件包括： \\n prefix.train_set.importance.pdf 训练集下，随机森林中，每个特征的重要性绘图 \\n prefix.train_set.importance.txt 训练集下，随机森林中，每个特征的重要性。有2个统计量（%IncMSE最重要），值越大越重要。\\n prefix.train_set.predict.txt    训练集下，随机森林分类结果（默认参数）\\n prefix.train_set.statistic.txt  训练集下，MSE\\n prefix.test_set.predict.txt     测试集下，训练模型预测结果\\n prefix.test_set.statistic.txt   测试集下，MSE \\n")

    parser$add_argument('-i', '--input', type='character', metavar='file', required=TRUE,
        help = "输入文件，第一列为样本名，第二列为因变量，往后为基因列。")

    parser$add_argument('-g', '--gene', type='character', metavar='gene_list', required=FALSE, default='所有基因',
        help = "输入要纳入随机森林分析的基因名称，多个基因用逗号分隔， [default: %(default)s]")

    parser$add_argument('-c', '--cov', type='character', metavar='cov', required=FALSE,
        help = "协变量输入文件，用于矫正。第一列样本名，后面所有列都是协变量。样本顺序没有要求。默认不用输入。 务必保证cov中包含input的所有样本，cov文件中的所有列都用于矫正。列名不能包含空格、“-”等特殊字符")

    parser$add_argument( '-o', '--output_dir', type='character', metavar='dir', default='./',
        help = "结果输出目录, 默认是当前路径。流程会自动创建。[default: %(default)s]。")

    parser$add_argument('--test', type='character', metavar='file', required=FALSE,
        help = "测试集输入文件,格式与input一样。使用 input数据训练出来的模型，对测试集数据做分类检验，必须包含 input/cov 中的所有特征（第一列也必须是样本的分类信息）。如果不提供，则不做")


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
library(A3)

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

# 建模
message("3. 随机森林建模")
f <- as.formula(paste('y ~', paste(vars, collapse = ' + ') ) )
 

modelRFP = rfPermute(f, data=data_clean, ntree=args$ntree, proximity=TRUE, importance=TRUE)  # 带有置换检验、检查importance pvalue功能
modelRF = as.randomForest(modelRFP)  # 转换为 随机森林模型
pred_value = predict(modelRF, data_clean, type='response')
pred_value = pred_value[rownames(data_clean)]
result_predict = data.frame(Sample=rownames(data_clean), Class=data_clean$y, PredictedClass=pred_value[rownames(data_clean)])
mse_train = sum((data_clean$y - pred_value)^2) / length(pred_value)
rsq_train = 1- (mse_train/ var(data_clean$y))

# A3 包补充分析

message("    基于A3包，计算R2/pvalue")
a3_result = a3(f, data_clean, randomForest, p.acc = 0.1,  n.folds = 10, features=FALSE)


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


    pred_value_test = predict(modelRF, data_clean_test, type='response')
    pred_value_test = as.numeric(as.character(pred_value_test[rownames(data_clean_test)]))

    result = data.frame(Sample=rownames(data_clean_test), Value=data_clean_test$y, PredictValue=pred_value_test)
    mse = sum((data_clean_test$y - pred_value_test)^2) / length(pred_value_test)
    rsq = 1- (mse/ var(data_clean_test$y))
    result_statistic = data.frame(MSE=mse, RSQ=rsq)
    # 预测结果输出 评价
    write.table( result, file.path(args$output_dir, paste0(args$prefix, '.test_set.predict.txt')), quote=F, row.names=F, col.names=T, sep='\t')
    write.table( result_statistic, file.path(args$output_dir, paste0(args$prefix, '.test_set.statistic.txt')), quote=F, row.names=F, col.names=T, sep='\t')
}


message("4. 结果输出")

# 结果输出
# 预测结果文件
file_predict  = file.path(args$output_dir, paste0(args$prefix, '.train_set.predict.txt'))  
write.table(result_predict, file_predict, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )

file_statistic = file.path(args$output_dir, paste0(args$prefix, '.train_set.statistic.txt'))  
result_statistic = data.frame(MSE=mse_train, RSQ=rsq_train, model.R2=a3_result$model.R2, model.p=a3_result$model.p)
write.table(result_statistic, file_statistic, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )

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

message("randomForest predict:    ", file_predict)
message("randomForest mse:        ", file_statistic)
message("randomForest importance:    ", file_importance)
message("randomForest importance:    ", pdf_file_importance)
