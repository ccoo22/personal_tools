#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: plot_risk.r -i <file> -o <pdf> --lasso_coef <file>  --exp <file> [ --rlib <dir>]

Options:
    -i, --input <file>              输入文件，必须含有表头 sample  risk    score   time_day        vital_status，第一列是样本，不允许缺失值
    --lasso_coef <file>             lasso系数，两列数据，第一列基因名，第二列系数
    --exp <file>                    样本的表达量文件，每一行一个基因，每一列一个样本，务必保证lasso_coef中的基因都存在，且input中的样本都存在
    -o, --output <pdf>              输出pdf文件 
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，lasso分析\n')
input               <- opts$input
lasso_coef          <- opts$lasso_coef
exp                 <- opts$exp
output              <- opts$output

rlib                <- opts$rlib
.libPaths(rlib)

# 导入 -> package
 library("pheatmap")
 library('gplots')
set.seed(91)

###################################################################### 主程序
message("read input")
data_input      = read.table(input, head = T, row.names = 1, check.names = F, sep = "\t", stringsAsFactors=F)
data_lasso_coef = read.table(lasso_coef, head = T, check.names = F, sep = "\t", stringsAsFactors=F)
data_exp        = read.table(exp, head = T, row.names = 1, check.names = F, sep = "\t")
data_input = data_input[order(data_input$score), ]  # 按照risk排序


# 绘图
pdf(output, width=14)

# (1) risk score 图
color = data_input$risk
color[color=='Low risk'] = 'green'
color[color=='High risk'] = 'red'
plot(1:nrow(data_input), data_input$score, xlab='', ylab='Risk score' , pch=19, col = color)

# (2) risk score 图
color = data_input$vital_status
color[color=='Alive'] = 'green'
color[color=='Dead'] = 'red'
plot(1:nrow(data_input), data_input$time_day, xlab='', ylab='Survival time (days)' , pch=19, col = color)
legend('topright', c('Dead', 'Alive'), col=c('red', 'green'), pch = 19, box.lwd = 0)

# (3) heatmap
gene = data_lasso_coef[data_lasso_coef$coef != 0, ]$feature
gene = gene[gene != '(Intercept)']
data_heatmap = data.matrix(data_exp[gene, rownames(data_input)])
pheatmap(data_heatmap, color=redgreen(75), show_colnames = FALSE, scale = 'row', cluster_cols=FALSE, cluster_rows=FALSE)

dev.off()


