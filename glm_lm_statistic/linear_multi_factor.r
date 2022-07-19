#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: linear_multi_factor.r  -i <file>  -o <file> [-g <string> -c <file> --rlib <dir>]

Options:
    -i, --input <file>      输入文件，第一列为样本名，第二列为因变量，往后为基因列。
    -g, --gene <string>     输入要纳入多元线性回归分析的基因名称，多个基因用逗号分隔， 默认是所有基因
    -c, --cov <file>        协变量输入文件，用于矫正。第一列样本名，后面所有列都是协变量。样本顺序没有要求。默认不用输入。 务必保证cov中包含input的所有样本，cov文件中的所有列都用于矫正。列名不能包含空格、“-”等特殊字符
    -o, --output <file>     输出文件前缀。 ./output
                            给出三个结果： output.txt 拟合模型
                                          output.pdf 因变量、拟合值散点图
                                          output.predicted.txt 因变量、拟合值结果
    --rlib <dir>            R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，多因素线性回归 \n')
input         <- opts$input
cov           <- opts$cov
output        <- opts$output
gene          <- opts$gene
rlib          <- opts$rlib
genes <- c()
if(!is.null(gene))  genes <- unlist(strsplit(gene, ","))

if(!is.null(rlib)) .libPaths(rlib)
 
set.seed(91)


# 读入
data_input = read.table(input, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "")
genes_replace = data.frame(from=colnames(data_input)[2:ncol(data_input)], newname=paste('gene', 2:ncol(data_input), sep=''), stringsAsFactors=F)  # 基因名字替换，方式有特殊字符而无法分析
rownames(genes_replace) = genes_replace[, 'newname']
colnames(data_input) <- c('y', genes_replace[,'newname'])


if(is.null(gene)) 
{
    genes <- colnames(data_input)[2:ncol(data_input)]
}else{
    # 使用指定的基因名
    tmp = genes_replace
    rownames(tmp) = tmp[, 'from']
    genes = tmp[genes, 'newname']
}

 

# 读入cov
data_cov = NULL
cov_names = c()
cov_names_replace = data.frame()
if(! is.null(cov))
{
    data_cov = read.table(cov, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "")
    cov_names_replace = data.frame(from=colnames(data_cov), newname=paste('cov', 1:ncol(data_cov), sep=''), stringsAsFactors=F)   
    rownames(cov_names_replace) = cov_names_replace[,'newname']
    colnames(data_cov) <- cov_names_replace[,'newname']
    cov_names = colnames(data_cov)

    if(sum(row.names(data_input) %in% row.names(data_cov)) != nrow(data_input) )
    {
        message("[Error]  输入的cov文件没有包含所有的input样本\n")
        q()  
    }
}

# 命名汇总
intercept = data.frame(from='(Intercept)', newname='(Intercept)')
rownames(intercept) = '(Intercept)'
names_replace = rbind(genes_replace, cov_names_replace, intercept)






# 开始分析
vars = genes  # 所有变量

data_tmp = data_input[, c('y', genes)]
if(! is.null(cov))
{
    data_tmp = cbind(data_tmp, data_cov[row.names(data_tmp), ])
    vars = c(vars, cov_names)
}

data_tmp = data_tmp[complete.cases(data_tmp), ]  # 去掉缺失值

# 检查数据是否异常
nmiss = nrow(data_tmp)
if(nmiss == 0 )
{   
    message("[Error] 数据错误，存在的缺失太多，跳过, skip. sample count = ", nmiss)
    next
}   

# 逻辑回归
f <- as.formula(paste('y ~', paste(vars, collapse = ' + ') ) )
lm_fit <- lm(f, data = data_tmp)
fit_summary <- summary(lm_fit)
lm_result  <- fit_summary$coefficients

# 拟合值
fit_y = fitted(lm_fit)
r_squared = round(fit_summary$r.squared, 2)

# 系数结果输出
result_model = data.frame(Var = names_replace[rownames(lm_result), 'from'], lm_result, check.names = F)
write.table(result_model, paste0(output, '.txt'), quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )

# 拟合值、实际值输出
result_y = data.frame(sample=rownames(data_tmp), y = data_tmp$y, predicted_value = fit_y)
write.table(result_y, paste0(output, '.predicted.txt'), quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )

# 绘图
pdf(paste0(output, '.pdf'))
plot(x = data_tmp$y, y = fit_y , pch = 16, cex = 0.4, xaxs = "i", yaxs = "i", xlab = 'dependent variable', ylab = 'predicted value')
abline(a = 0, b = 1) # 直线
legend("bottomright", legend = paste("R2 = ", r_squared))
dev.off()

