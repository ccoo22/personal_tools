#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: logistics_single_factor_predict.r  -i <file> -m <file>  -p <prefix> [ --rlib <dir>]

Options:
    -i, --input <file>      输入文件，第一列为样本名，第二列为样本分组（case=1，control=0），用于计算ROC，往后为基因列。该软件对 --model的每一个模型做预测。务必保证input的列名与model的行名一致
    -m, --model <file>      逻辑回归模型数据，第一列必须是基因ID, 另外，必须包含以下三列数据Estimate(argument)， Estimate(Intercept)，threshold， 用于预测样本分组
    -p, --prefix <prefix>   输出文件前缀， 生成:  (1) 逻辑回归分类结果 prefix.condition.txt (3) 基于逻辑回归模型对输入样本做预测的详细情况 prefix.predict.txt
    --rlib <dir>            R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，逻辑回归模型应用\n')
input         <- opts$input
model         <- opts$model
prefix        <- opts$prefix
rlib          <- opts$rlib

if(!is.null(rlib)) .libPaths(rlib)
 

# 读入
data = read.table(input, head = T, row.names = 1, check.names = F, sep="\t")
model = read.table(model, head = T, row.names = 1, check.names = F, sep="\t")
colnames(data)[1] = 'group'
genes = rownames(model)

# 检查关键列是否存在
if(sum( c('Estimate(argument)', 'Estimate(Intercept)', 'threshold') %in% colnames(model)) != 3  )
{
    message("[error] model文件必须含有以下三个表头： 'Estimate(argument)', 'Estimate(Intercept)', 'threshold'")
    q()
}

# result
result = matrix(ncol = 6, nrow = nrow(model))  
colnames(result) = c('gene', 'nmiss case', 'nmiss control', "specificity", "sensitivity", "accuracy")

# 保留预测值
result_pred <- matrix(ncol = length(genes) * 2 + 2, nrow = nrow(data))  
rownames(result_pred) = rownames(data)
colnames(result_pred) = c('sample', 'group', paste0(genes, "-predict"), paste0(genes, "-class"))
result_pred[, 'sample'] = rownames(data)
result_pred[, 'group']  = data[, 1]

for(row in 1:length(genes))
{
    gene = genes[row]
    message("process : ", gene)
    # 模型
    a         = model[gene, 'Estimate(argument)']
    b         = model[gene, 'Estimate(Intercept)']
    threshold = model[gene, 'threshold']

    if(is.na(a) | is.na(b) | is.na(threshold))
    {   
        message("    [warning] model 数据为NA, 无法分析： ", gene)
        result[row, 1:2] = c(gene, 'NA')
        next
    }

    if(sum(gene %in% colnames(data)) == 0)
    {   
        message("    [warning] input文件里缺失基因： ", gene)
        result[row, 1:2] = c(gene, 0)
        next
    }
    
    # 数据
    data_tmp = data[, c('group', gene)]
    data_tmp = data_tmp[complete.cases(data_tmp), ]  # 去掉缺失值

    nmiss_case = sum(data_tmp$group == 1)
    nmiss_control = sum(data_tmp$group == 0)
    if(nmiss_case == 0 & nmiss_control == 0)
    {
        message("    [warning] input文件里基因的表达量缺失： ", gene)
        result[row, 1:2] = c(gene, 0)
        next
    }
    
    x = data_tmp[, 2]
    names(x) = rownames(data_tmp)

    # 预测值
    y = 1 / ( 1 + exp(- a * x - b ))
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

    result[row, ] = c(gene, nmiss_case, nmiss_control, spe, sen, acc)
    result_pred[names(y), paste0(gene, '-predict')] = y
    result_pred[names(y), paste0(gene, '-class')] = class
}


# 结果输出
condition_file = paste0(prefix, '.condition.txt')
pred_file = paste0(prefix, '.predict.txt')

write.table(result, condition_file, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )
write.table(result_pred, pred_file, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )


message("model condition: ", condition_file)
message("model predict:   ", pred_file)

