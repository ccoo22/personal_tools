#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: plot_correlation.r.r  -i <file> -o <pdf file> [ --no0 --lg10 --xlim <string> --ylim <string>]

Options:
    -i, --input <file>              输入文件，第一列为基因名，第二列(x)、第三列(y)分别为两个需要做比较的样本。注：空值会被去掉
    -o, --output <pdf file>         输出pdf文件路径,示例：./a.pdf
    --no0                           去掉表达量为0的基因
    --lg10                          对数据取log10
    --xlim <string>                 x轴显示范围,默认不限制。 示例： 0,5
    --ylim <string>                 y轴显示范围,默认不限制。 示例： 0,5 " -> doc

opts   <- docopt(doc, version='甘斌，热图\n')
input              <- opts$input
output             <- opts$output
no0                <- opts$no0
lg10               <- opts$lg10
xlim               <- opts$xlim
ylim               <- opts$ylim

# 参数预处理
if(! is.null(xlim)) xlim = as.numeric(unlist(strsplit(xlim, ','))) 
if(! is.null(ylim)) ylim = as.numeric(unlist(strsplit(ylim, ','))) 


message("读入文件")
data = read.table(input, head = T , row.names = 1, sep = '\t', check.names=F);
data_clean = data[complete.cases(data), ]

# 去掉0
if(no0 == TRUE)
{   
    message("去除表达量为0的行 : ", appendLF = FALSE)
    no0_row = apply(data_clean, 1, function(x){! 0 %in% x}) # 不含有0的行
    is0_count = nrow(data_clean) - sum(no0_row)
    message(is0_count, "个/", nrow(data_clean))
    data_clean = data_clean[no0_row, ]
}

# 取 log10
if(lg10 == TRUE)
{   
    message("表达量取log10")
    # 存在负数，不能取log
    if(sum(data_clean < 0) > 0)
    {
        message("存在负数，不能取log10")
        q()
    }
    data_clean = log10(data_clean + 0.00000000001) # 加上很小的数值后，再取log, 防止存在0
}

# 开始分析
sample_name1 = colnames(data_clean)[1]  # 样本名1
sample_name2 = colnames(data_clean)[2]  # 样本名2

# 相关性计算
message("计算相关性")
cor_coefficient = round(cor(data_clean[, 1], data_clean[, 2], method = "pearson" ), 4)

# 线性拟合，sample1作为x, sample2作为y
message("线性拟合")
cor_fit <- lm(data_clean[, 2] ~ data_clean[, 1] + 1)

# 绘图
message("绘图")
pdf(output);
# 横纵坐标名称
x_name = ifelse(lg10 == TRUE, paste(sample_name1, "(log10)"), sample_name1)
y_name = ifelse(lg10 == TRUE, paste(sample_name2, "(log10)"), sample_name2)

plot(x = data_clean[, 1], y = data_clean[, 2], xlab = x_name, ylab = y_name , pch = 16, xlim = xlim, ylim = ylim, cex = 0.4, xaxs = "i", yaxs = "i")
lines(data_clean[, 1], fitted(cor_fit))#添加拟合值对x的散点图并连线
legend("topleft", legend = paste("r = ", cor_coefficient))
dev.off()
