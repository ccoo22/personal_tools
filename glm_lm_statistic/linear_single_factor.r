#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: linear_single_factor.r  -i <file> -o <dir> [ -c <file> --rlib <dir> --if_dotplot <logic>]

Options:
    -i, --input <file>      输入文件，第一列为样本名，第二列为因变量，往后为基因列。该软件对每一个基因单独做 线性回归 分析。 列名不能包含空格、“-”等特殊字符
    -c, --cov <file>        协变量输入文件，用于矫正。第一列样本名，后面所有列都是协变量。样本顺序没有要求。默认不用输入。 务必保证cov中包含input的所有样本，cov文件中的所有列都用于矫正。列名不能包含空格、“-”等特殊字符
    -o, --output <file>     结果输出目录
    --if_dotplot <logic>    是否绘制散点图 TRUE/FALSE [default: FALSE]
    --rlib <dir>            R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，单因素线性回归 \n')
input       <- opts$input
cov         <- opts$cov
output      <- opts$output
if_dotplot  <- opts$if_dotplot
rlib        <- opts$rlib

if(!is.null(rlib)) .libPaths(rlib)

if(!file.exists(output)){
    dir.create(output)
}

set.seed(91)

# 读入 input
data_input = read.table(input, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "")
genes_replace = data.frame(gene=colnames(data_input)[2:ncol(data_input)], newname=paste('gene', 2:ncol(data_input), sep=''), stringsAsFactors=F)  # 基因名字替换，方式有特殊字符而无法分析
y_name = colnames(data_input)[1]
rownames(genes_replace) = genes_replace[, 'newname']
colnames(data_input) <- c('y', genes_replace[,'newname'])
genes <- colnames(data_input)[2:ncol(data_input)]

# 读入cov
data_cov = NULL
cov_names = c()
cov_names_replace = data.frame()

if(! is.null(cov))
{
    data_cov = read.table(cov, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "")
    cov_names_replace = data.frame(cov=colnames(data_cov), newname=paste('cov', 1:ncol(data_cov), sep=''), stringsAsFactors=F)   
    rownames(cov_names_replace) = cov_names_replace[,'newname']
    colnames(data_cov) <- cov_names_replace[,'newname']
    cov_names = colnames(data_cov)

    if(sum(row.names(data_input) %in% row.names(data_cov)) != nrow(data_input) )
    {
        message("[Error]  输入的cov文件没有包含所有的input样本\n")
        q()  
    }

}

# 保留模型值
result <- matrix(ncol = 7 + 2 * length(cov_names), nrow = length(genes))
col_names = c('Var', 'NMISS', "Estimate", "Pvalue", "Pvalue_FDR", "Estimate(Intercept)", "Pvalue(Intercept)")  # 必要结果
for(cov_name in cov_names)  # 协变量结果
{
    col_names = c(col_names, paste0(c("Estimate", "Pvalue"), "(", cov_names_replace[cov_name, 'cov'], ")") )
}

colnames(result) = col_names

if(if_dotplot == "TRUE") cairo_pdf(paste0(output,'/dotplot.pdf'), onefile=T)

for( col in 1:length(genes))
{   
    gene = genes[col]
    vars = c(gene)  # 变量列表

    message("process : ", genes_replace[gene, 'gene'])

    # 准备输入数据、cov数据
    data_tmp = data_input[, c('y', gene)]
    if(! is.null(cov))
    {
        data_tmp = cbind(data_tmp, data_cov[row.names(data_tmp), ])
        vars = c(vars, cov_names)
    }

    data_tmp = data_tmp[complete.cases(data_tmp), ]  # 去掉缺失值

    # 检查数据是否异常
    nmiss = nrow(data_tmp)
    if(nmiss == 0 | sum(!duplicated(data_tmp[,2])) == 1)
    {   
        message("[Error] 数据错误，存在的缺失太多或者x值只有一种数据，跳过, skip. sample count = ", nmiss)
        result[col, 1:2] = c(genes_replace[gene, 'gene'], nmiss) 
        next
    }   

    # 逻辑回归
    mouse.regression<-
    formula <- as.formula(paste('y ~ ', paste(vars, collapse = ' + ') ) )
    lm_fit <- lm(formula, data=data_tmp)
    fit_summary <- summary(lm_fit)
    lm_result  <- fit_summary$coefficients

    if(if_dotplot == "TRUE"){
        data_plot = data.frame(x = data_tmp[, gene], y = data_tmp[, "y"], fitted_y = fitted(lm_fit)[rownames(data_tmp)])
        data_plot = data_plot[order(data_plot[, "x"]), , drop = F]  # 按x值从小到大排序，避免系数为负值曲线重叠
        x = data_plot[, "x"]
        y = data_plot[, "y"]
        plot(x = x, y = y, xlab = genes_replace[gene, 'gene'], ylab = y_name , pch = 16, main = paste0(genes_replace[gene, 'gene'], " ", y_name) ) # 黑色散点
        lines(x, data_plot$fitted_y, col="red") # 拟合曲线
        
        text_formula <- paste("y = ", format(lm_fit$coefficients[gene], scientific=TRUE, digit=4), "x + ", format(lm_fit$coefficients['(Intercept)'], scientific=TRUE, digit=4) , sep="")
        if(! is.null(cov)){
            text_formula = paste0(text_formula, " + ...")
        }
        text_rsquared <- paste("R-squared = ", round(fit_summary$r.squared, 3), sep = "")
        legend("topright", legend = c(text_formula, text_rsquared))
    }

    # 结果记录
    result_value = c(genes_replace[gene, 'gene'], nmiss, lm_result[gene, 'Estimate'], lm_result[gene, 'Pr(>|t|)'], NA ,lm_result['(Intercept)', 'Estimate'], lm_result['(Intercept)', 'Pr(>|t|)'])
    for(cov_name in cov_names)  # 协变量结果
    {   
        result_value = c(result_value, lm_result[cov_name, 'Estimate'], lm_result[cov_name, 'Pr(>|t|)']) 
    }
    
    result[col, ] <- result_value
}
if(if_dotplot == "TRUE") dev.off()

result[,'Pvalue_FDR'] =  p.adjust(result[,'Pvalue'], method = "fdr", n=nrow(result))
# 结果输出
write.table(result, paste0(output, '/linear_single_factor_result.txt'), quote = FALSE, row.names = FALSE, sep = '\t', col.names = T )

