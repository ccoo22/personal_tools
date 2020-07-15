#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: linear_multi_factor.r  -i <file>  -o <file> [-g <string> -c <file> --rlib <dir>]

Options:
    -i, --input <file>      输入文件，第一列为样本名，第二列为因变量，往后为基因列。
    -g, --gene <string>     输入要纳入多元线性回归分析的基因名称，多个基因用逗号分隔， 默认是所有基因
    -c, --cov <file>        协变量输入文件，用于矫正。第一列样本名，后面所有列都是协变量。样本顺序没有要求。默认不用输入。 务必保证cov中包含input的所有样本，cov文件中的所有列都用于矫正。列名不能包含空格、“-”等特殊字符
    -o, --output <file>     输出文件。 ./output.txt
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
data_input = read.table(input, head = T, row.names = 1, check.names = F, sep="\t")
colnames(data_input)[1] <- 'y'

if(is.null(gene)) genes <- colnames(data_input)[2:ncol(data_input)]

if(sum(grepl(" ", colnames(data_input))) > 0 |  sum(grepl("-", colnames(data_input))))
{
    message("[Error]  输入的input文件表头包含空格或者-\n")
    q()
}

# 读入cov
data_cov = NULL
cov_names = c()
if(! is.null(cov))
{
    data_cov = read.table(cov, head = T, row.names = 1, check.names = F, sep="\t")
    if(sum(grepl(" ", colnames(data_cov))) > 0 |  sum(grepl("-", colnames(data_cov))))
    {
        message("[Error]  输入的cov文件表头包含空格或者-\n")
        q()
    }
    if(sum(row.names(data_input) %in% row.names(data_cov)) != nrow(data_input) )
    {
        message("[Error]  输入的cov文件没有包含所有的input样本\n")
        q()  
    }
    cov_names = colnames(data_cov)
}




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



# 结果输出
result_model = data.frame(Var = rownames(lm_result), lm_result, check.names = F)

write.table(result_model, output, quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )




