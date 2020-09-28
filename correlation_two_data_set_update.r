#!/home/genesky/software/r/3.5.1/bin/Rscript

.libPaths("/home/genesky/software/r/3.5.1/lib64/R/library/")
library(docopt)
"Usage: correlation_two_data_set.r --data1 <file> --data2 <file> --sample_file <file>  --output <file>  [--method <string>  --plot_top_data1 <int> --plot_top_data2 <int> --pdf_width <numeric>  --pdf_height <numeric>] 
Options:
    --data1 <file>           数据表1，第一列是特征名称。每一行是一个特征，每一列是一个样本，含有表头。必须保证data1/data2中的样本名与 --sample_file 参数一致。对样本顺序没有要求。
    --data2 <file>           数据表2，第一列是特征名称。每一行是一个特征，每一列是一个样本，含有表头
    --sample_file <file>     样本文件，一列数据，有表头。对包含的样本做相关性分析
    --output <file>          输出文件
                             同时会生成矩阵形式的.estimate_matrix.txt、.pvalue_matrix.txt文件，方便绘制热图
    --plot_top_data1 <int>   对相关性结果绘制热图，选择data1数据中pvalue最显著的前n个 [default: 100]
    --plot_top_data2 <int>   对相关性结果绘制热图，选择data1数据中pvalue最显著的前n个 [default: 100]
    --pdf_width <numeric>    pdf图宽度 [default: 7]
    --pdf_height <numeric>   pdf图高度 [default: 7]
    --method <string>        分析方法，仅支持 pearson/spearman [default: spearman]" -> doc

opts                <- docopt(doc, version = 'Program : correlation_two_data_set v1.0 \n          甘斌 129\n')
data1_file          <- opts$data1
data2_file          <- opts$data2
sample_file         <- opts$sample_file
output              <- opts$output
method              <- opts$method
plot_top_data1      <- as.numeric(opts$plot_top_data1)
plot_top_data2      <- as.numeric(opts$plot_top_data2)
pdf_width           <- as.numeric(opts$pdf_width)
pdf_height          <- as.numeric(opts$pdf_height)

library(corrplot)

message(method)
# 数据读入
message("loading data : ", data1_file)
data1   = read.table(data1_file, header = T, check.names = F, stringsAsFactors=F, sep = "\t")
rownames(data1) = data1[,1]

message("loading data : ", data2_file)
data2   = read.table(data2_file, header = T, check.names = F, stringsAsFactors=F, sep = "\t")
rownames(data2) = data2[,1]

message("loading data : ", sample_file)
data_sample = read.table(sample_file, stringsAsFactors=F, sep = "\t", header = T, colClasses = 'character')
samples = data_sample[,1]
 
# 检查样本、特征名称是否异常
message("check sample name and feature")
if( sum(samples %in% colnames(data1)) != length(samples) )
{
    message('[Error] 数据集1 中缺失部分样本,请仔细核对')
    print(samples[! samples %in% colnames(data1)])
    q()
}
if( sum(samples %in% colnames(data2)) != length(samples) )
{
    message('[Error] 数据集2 中缺失部分样本,请仔细核对')
    print(samples[! samples %in% colnames(data2)])
    q()
}

# 开始计算相关性
message("start correlation")
result = matrix(nrow = nrow(data1) * nrow(data2), ncol =  4);
colnames(result) = c('feature1', 'feature2', 'estimate', 'pvalue')
 
result_estimate = matrix(nrow = nrow(data1), ncol = nrow(data2))
result_pvalue   = matrix(nrow = nrow(data1), ncol = nrow(data2))
 
rownames(result_estimate) = rownames(data1)
colnames(result_estimate) = rownames(data2)
 
rownames(result_pvalue) = rownames(data1)
colnames(result_pvalue) = rownames(data2)



row = 0
for(name_d1 in rownames(data1))
{   
    value1 = as.numeric(data1[name_d1, samples])
    for(name_d2 in rownames(data2))
    {   
        row = row + 1
        value2 = as.numeric(data2[name_d2, samples])
 
        tmp <- cor.test(value1, 
                    value2, 
                    alternative = "two.sided", 
                    method = method, 
                    conf.level = 0.95)
        pvalue   <- tmp$p.value
        estimate <- tmp$estimate[[1]]

        result[row, ] = c(name_d1, name_d2, estimate, pvalue)
        result_estimate[name_d1, name_d2] = estimate
        result_pvalue[name_d1, name_d2] = pvalue
    }
}
write.table(result, output , sep = "\t", quote = F, row.names = F ,col.names = T)

result_estimate = data.frame(id = rownames(result_estimate), result_estimate, check.names=F)
result_pvalue = data.frame(id = rownames(result_pvalue), result_pvalue, check.names=F)
write.table(result_estimate, paste0(output, ".estimate_matrix.txt"), sep = "\t", quote = F, row.names = F ,col.names = T)
write.table(result_pvalue, paste0(output, ".pvalue_matrix.txt") , sep = "\t", quote = F, row.names = F ,col.names = T)

# 绘图准备
message("start plot")
result_df = data.frame(feature1=result[, 'feature1'], feature2=result[, 'feature2'], estimate = as.numeric(result[, 'estimate']), pvalue = as.numeric(result[, 'pvalue']), stringsAsFactors=F)
result_df = result_df[order(result_df$pvalue), ]  # pvalue排序

# 选择绘图特征
if(nrow(data1) < plot_top_data1) plot_top_data1 = nrow(data1)
if(nrow(data2) < plot_top_data2) plot_top_data2 = nrow(data2)

feature1_plot = unique(result_df$feature1)[1:plot_top_data1]
feature2_plot = unique(result_df$feature2)[1:plot_top_data2]
result_estimate_plot = result_estimate[feature1_plot, feature2_plot, drop = F]  # 绘图系数
result_pvalue_plot   = result_pvalue[feature1_plot, feature2_plot, drop = F]  # 绘图pvalue
 
# 颜色
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                                "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("blue", "white", "red"))
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                                "cyan", "#007FFF", "blue", "#00007F"))


# corrplot 绘图
 
pdf(file=paste0(output, ".top.pdf"), width=pdf_width, height=pdf_height)

corrplot(as.matrix(result_estimate_plot),
        p.mat = as.matrix(result_pvalue_plot),
        method = 'circle',
        type = 'full',
        title = ' ',
        tl.col = 'black',
        is.corr = FALSE,
        cl.pos  = TRUE,
        col = col3(200),
        insig = "label_sig",
        sig.level = c(.001, .01, .05),
        pch.cex = 0.9
        )
dev.off()


