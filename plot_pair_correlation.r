#!/home/genesky/software/r/4.1.2/bin/Rscript

library(docopt)

"Usage: plot_pair_correlation.r --file1 <file> --file2 <file> -o <file> --pdf <pdf>  [--log10_1 --log10_2 --log10_1_offset <numeric>  --log10_2_offset <numeric> --t1 --t2 --select_row <file> --select_file1 <file> --select_file2 <file> --width <int> --height <int> --rlib <dir>]

Options:
    --file1 <file>                  输入文件，第一列为样本名称，第二列及之后为每一个特征的值
    --file2 <file>                  输入文件，第一列为样本名称，第二列及之后为每一个特征的值
                                    注：file1 和 file2 的第一列行名一定要有重叠，以保证配对。允许两个文件有不同。
                                    注：默认使用 file1 file2 共有的行名作为配对结果
                                    注：默认使用 file1/file2 所有的列两两进行比较
                                    注：相关性计算采用 pearson。线性拟合采用 lm
    -o, --output <file>             输出text文件路径，例如 result.txt
    --pdf <pdf>                     输出pdf文件路径，例如 result.pdf
    --t1                            对file1 数据做转置
    --t2                            对file2 数据做转置
                                    如果 file1/file2 的格式不满足上面的描述要求，可以通过转置的方式使其满足
    --log10_1                       是否对file1的数据取 log10
    --log10_2                       是否对file2的数据取 log10
    --log10_1_offset <numeric>      对file1取log10时，为防止有0存在，要对数据加上一个偏差量, 例如 1 [default: 0]
    --log10_2_offset <numeric>      对file2取log10时，为防止有0存在，要对数据加上一个偏差量, 例如 1 [default: 0]
    --select_row <file>             通过文件指定选择共有的行纳入分析。默认使用 file1 file2 共有的行名作为配对结果
                                    注：文件只有一列，一行一个名称，不要有表头
    --select_file1 <file>           从file1中选择指定的列进行分析，而不是用所有的。文件只有一列，一行一个名称，不要有表头。
    --select_file2 <file>           从file2中选择指定的列进行分析，而不是用所有的。文件只有一列，一行一个名称，不要有表头。
    --width <int>                   pdf宽度 [default: 7]
    --height <int>                  pdf高度 [default: 7]
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/4.1.2/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，比较两组数据的相关性\n')
file1              <- opts$file1
file2              <- opts$file2
output             <- opts$output
pdf_file           <- opts$pdf
t1                 <- opts$t1
t2                 <- opts$t2
log10_1            <- opts$log10_1
log10_2            <- opts$log10_2
log10_1_offset     <- as.numeric(opts$log10_1_offset)
log10_2_offset     <- as.numeric(opts$log10_2_offset)
select_row         <- opts$select_row
select_file1       <- opts$select_file1
select_file2       <- opts$select_file2
width              <- as.integer(opts$width)
height             <- as.integer(opts$height)
rlib               <- opts$rlib


# 开始分析
message('(1) 读入 data1 data2')
data1_raw <- read.table(file1, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "")
data2_raw <- read.table(file2, header = T, sep = "\t" , row.names = 1, check.name = F, stringsAsFactors = F, quote = "", comment.char = "")

# 确认是否转置处理
if(t1)
{
    data1_raw = t(data1_raw)
}
if(t2)
{
    data2_raw = t(data2_raw)
}


message('(2) 确定要分析的共有的行名')
common_row = c()
if(is.null(select_row))
{   
    common_row = intersect(rownames(data1_raw), rownames(data2_raw))# 交集
    message("        取交集，共有的行名数量：", length(common_row))
}else{
    select_row_data <- read.table(select_row, header = F, sep = "\t" , check.name = F, stringsAsFactors = F, quote = "", comment.char = "", colClasses = 'character')
    common_row = select_row_data[,1]
    diff1 = setdiff(common_row, rownames(data1_raw))
    diff2 = setdiff(common_row, rownames(data2_raw))
    if(length(diff1) > 0)
    {
        message("        [Error] select_row 文件异常，部分名称在 file1 中没有找到： ", paste0(diff1, collapse=','))
        q()
    }
    if(length(diff2) > 0)
    {
        message("        [Error] select_row 文件异常，部分名称在 file2 中没有找到： ", paste0(diff2, collapse=','))
        q()
    }
    message("        你自定义纳入分析的共有的行名数量：", length(common_row))
}

message('(3) 确定要分析的列')

file1_column = c()
file2_column = c()
if(is.null(select_file1))
{
    file1_column = colnames(data1_raw)
    message("        file1 要分析的列数量：", length(file1_column))
}else{
    select_file1_data <- read.table(select_file1, header = F, sep = "\t" , check.name = F, stringsAsFactors = F, quote = "", comment.char = "", colClasses = 'character')
    file1_column = select_file1_data[,1]
    diff1 = setdiff(file1_column, colnames(data1_raw))
    if(length(diff1) > 0)
    {
        message("        [Error] select_file1 文件异常，部分名称在 file1 中没有找到： ", paste0(diff1, collapse=','))
        q()
    }
    message("        你自定义纳入file1分析的列数量：", length(file1_column))
}
 
if(is.null(select_file2))
{
    file2_column = colnames(data2_raw)
    message("        file1 要分析的列数量：", length(file2_column))
}else{
    select_file2_data <- read.table(select_file2, header = F, sep = "\t" , check.name = F, stringsAsFactors = F, quote = "", comment.char = "", colClasses = 'character')
    file2_column = select_file2_data[,1]
    diff2 = setdiff(file2_column, colnames(data2_raw))
    if(length(diff2) > 0)
    {
        message("        [Error] select_file2 文件异常，部分名称在 file2 中没有找到： ", paste0(diff2, collapse=','))
        q()
    }
    message("        你自定义纳入file2分析的列数量：", length(file2_column))
}
 
message('(4) 开始分析')
data1_clean = data1_raw[common_row, file1_column, drop=F]
data2_clean = data2_raw[common_row, file2_column, drop=F]
if(log10_1)
{
    data1_clean = log10(data1_clean + log10_1_offset)
}
if(log10_2)
{
    data2_clean = log10(data2_clean + log10_2_offset)
}


pdf(pdf_file, width=width, height=height)
result = data.frame(matrix(nrow=0, ncol=5))
for(name1 in colnames(data1_clean))
{
    for(name2 in colnames(data2_clean))
    {   
        message('        name1=',name1, '    /    name2=', name2)
        data_tmp = data.frame(x=data1_clean[,name1], y=data2_clean[,name2])
        # 删除缺失值
        data_tmp = data_tmp[complete.cases(data_tmp), ]
        if(nrow(data_tmp) < 3)
        {
            message('            缺失值太多，删除缺失值后，剩下的数据量不足3，不适合分析')
            next
        }

        # 线性回归
        model<- lm('y ~ x ', data=data_tmp) 
        summary_res = summary(model)

        # 相关性分析
        cor_res = cor.test(data_tmp[, 'x'], y = data_tmp[, 'y'], method = 'pearson')

        result = rbind(result, c(name1, name2, nrow(data_tmp), cor_res$estimate, cor_res$p.value))
        # 绘图
        plot(x = data_tmp[, 'x'], y = data_tmp[, 'y'], xlab = name1, ylab = name2 , pch = 16, cex = 0.4 )
        # plot(x = data_tmp[, 'x'], y = data_tmp[, 'y'], xlab = name1, ylab = name2 , pch = 16, cex = 0.4, xaxs = "i", yaxs = "i")
        lines(data_tmp[, 'x'], fitted(model)) # 添加拟合值对x的散点图并连线
        legend("topright", legend = paste("R = ", round(cor_res$estimate, 4), "\nPvalue = ", round(cor_res$p.value, 4)) )
         
    }
}
dev.off()

colnames(result) = c('name1', 'name2', 'sample_count', 'pearson_estimate', 'pearson_pvalue')
write.table(result, output, quote=F, row.names=F, col.names=T, sep='\t')

# inputfile  = normalizePath(ARGS[1])
# outputpdf  =  ARGS[2]

# data = read.table(inputfile, head = T , row.names = 1, sep = '\t');
# data = log10(data + 0.00000000001)
# sample1 = colnames(data)[1]
# sample2 = colnames(data)[2]

# cor_count = round(cor(data[, 1], data[, 2], method = "pearson" ), 4)
# z<- lm(data[, 2]~data[, 1]+1)


# pdf(outputpdf);
# # plot(x = data[, 1], y = data[, 2], xlab = paste(sample1, "(log10)"), ylab = paste(sample2, "(log10)") , pch = 16, cex = 0.4, xaxs = "i", yaxs = "i")
# plot(x = data[, 1], y = data[, 2], xlab = paste(sample1, "(log10)"), ylab = paste(sample2, "(log10)") , pch = 16, xlim = c(0, 5), ylim = c(0, 5), cex = 0.4, xaxs = "i", yaxs = "i")
# lines(data[, 1], fitted(z))#添加拟合值对x的散点图并连线
# legend("topleft", legend = paste("r = ", cor_count))
# dev.off()




