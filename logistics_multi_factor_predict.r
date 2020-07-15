#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: logistics_multi_factor_predict.r  -i <file> -m <file> -t <float> -p <prefix> [ --rlib <dir>]

Options:
    -i, --input <file>      输入文件，第一列为样本名，第二列为样本分组（case=1，control=0），往后为基因列。该软件 对 --model 里的基因做预测。 该文件要包含model文件中的所有基因。
    -m, --model <file>      逻辑回归模型数据，第一列必须是基因ID, 另外，必须包含 Estimate 列信息。常数项的基因名称为 “(Intercept)”，也必须存在
    -p, --prefix <prefix>   输出文件前缀， 生成:  (1) 逻辑回归分类结果 prefix.condition.txt (3) 基于逻辑回归模型对输入样本做预测的详细情况 prefix.predict.txt
    -t, --threshold <float> 分类阈值， [0-1], 通常选择多元逻辑回归模型下的best值对应的 threshold
    --rlib <dir>            R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，多元逻辑回归模型应用\n')
input         <- opts$input
model         <- opts$model
threshold     <- as.numeric(opts$threshold)
prefix        <- opts$prefix
rlib          <- opts$rlib
 
if(!is.null(rlib)) .libPaths(rlib)
 

# 读入
data = read.table(input, head = T, row.names = 1, check.names = F, sep="\t")
model = read.table(model, head = T, row.names = 1, check.names = F, sep="\t")
colnames(data)[1] = 'group'

# 模型用到的基因
genes = rownames(model)   
genes = genes[genes != "(Intercept)"]

# 检查关键列是否存在
if(! '(Intercept)' %in% rownames(model))
{
    message("[Error] model文件里缺少 (Intercept) 常数项")
    q()
}
if(! 'Estimate' %in% colnames(model))
{
    message("[Error] model文件里缺少 Estimate 列信息")
    q()
}

if(sum( genes %in% colnames(data)) != length(genes)  )
{
    message("[Error] model中的基因在input文件中缺失")
    q()
}


# 预测

data_tmp = data[, c('group', genes)]
data_tmp = data_tmp[complete.cases(data_tmp), ]  # 去掉缺失值

# 缺失检测
nmiss_case = sum(data_tmp$group == 1)
nmiss_control = sum(data_tmp$group == 0)
if(nmiss_case == 0 & nmiss_control == 0)
{
    message("    [warning] input文件里基因的表达量缺失")
    q()
}

a = model[genes, 'Estimate']
b = model["(Intercept)", 'Estimate']
y = apply(data_tmp[,2:ncol(data_tmp)], 1, function(x){ 1 / ( 1 + exp(- sum(a * x) - b ))   }  )


class = y
class[y > threshold] = 1
class[y < threshold] = 0

# 统计当前数据下的分类效果
spe = 'NA'
sen = 'NA'
acc = sum(class == data_tmp$group) / length(class)
if(sum(data_tmp$group == 0) > 0)
{
    spe = sum(class[data_tmp$group == 0] == 0) / sum(data_tmp$group == 0)
} 
if(sum(data_tmp$group == 1) > 0)
{
    sen = sum(class[data_tmp$group == 1] == 1) / sum(data_tmp$group == 1)
} 

# 输出
result_pred = data.frame(sample = rownames(data_tmp), group = data_tmp[,'group'], predict = y, classify = class)
result = matrix(c(nmiss_case, nmiss_control, spe, sen, acc), 1, 5)
colnames(result) = c('nmiss case', 'nmiss control', "specificity", "sensitivity", "accuracy")

condition_file = paste0(prefix, '.condition.txt')
pred_file = paste0(prefix, '.predict.txt')
write.table(result, condition_file, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )
write.table(result_pred, pred_file, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )


message("model condition: ", condition_file)
message("model predict:   ", pred_file)

