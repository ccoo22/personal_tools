#!/home/genesky/software/r/3.5.1/bin/Rscript

.libPaths("/home/genesky/software/r/3.5.1/lib64/R/library/")
library(docopt)
"Usage: correlation_two_data_set.r --data1 <file> --data2 <file> --compare <file> --sample <string>  --output <file>  [--method <string> --keep1_name <string>  --keep2_name <string> ] 
Options:
   --data1 <file>          数据表1，第一列是特征名称。每一行是一个特征，每一列是一个样本，含有表头。必须保证data1/data2中的样本名与--sample参数一致。对样本顺序没有要求。
   --data2 <file>          数据表2，第一列是特征名称。每一行是一个特征，每一列是一个样本，含有表头
   --compare <file>        两列特征名称数据，表示需要做相关性分析的特征配对。没有表头，第一列是data1中的特征名字，第二列是data2中的特征名字。
   --sample <string>       样本列表，逗号分隔
   --output <file>         输出文件
   --keep1_name <string>   数据表1中的指定列输出到结果中, 多个列名用逗号分隔， 默认仅保留特征名称
   --keep2_name <string>   数据表2中的指定列输出到结果中, 多个列名用逗号分隔， 默认仅保留特征名称
   --method <string>       分析方法，仅支持 pearson/spearman [default: pearson]" -> doc

opts                <- docopt(doc, version = 'Program : correlation_two_data_set v1.0 \n          甘斌 129\n')
data1_file          <- opts$data1
data2_file          <- opts$data2
compare_file        <- opts$compare
sample              <- opts$sample
output              <- opts$output
method              <- opts$method 
keep1_name          <- opts$keep1_name 
keep2_name          <- opts$keep2_name 




# 数据读入
message("loading data : ", data1_file)
data1   = read.table(data1_file, head = T, check.names = F, stringsAsFactors=F, sep = "\t")
rownames(data1) = data1[,1]

message("loading data : ", data2_file)
data2   = read.table(data2_file, head = T, check.names = F, stringsAsFactors=F, sep = "\t")
rownames(data2) = data2[,1]

message("loading data : ", compare_file)
compare = read.table(compare_file, stringsAsFactors=F, sep = "\t")
samples = unlist(strsplit(sample, ','))

# 确定需要保存的列
if(is.null(keep1_name))
{
    keep1_name = colnames(data1)[1]
}else {
   keep1_name = unlist(strsplit(keep1_name, ','))
}
if(is.null(keep2_name))
{
    keep2_name = colnames(data2)[1]
}else {
   keep2_name = unlist(strsplit(keep2_name, ','))
}

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
# if( sum(compare[,1] %in% rownames(data1)) != nrow(compare) )
# {
#     message('[Error] 数据集1 中缺失部分特征,请仔细核对')
#     print(compare[,1][! compare[,1] %in% rownames(data1)])
#     q()
# }
# if( sum(compare[,2] %in% rownames(data2)) != nrow(compare) )
# {
#     message('[Error] 数据集2 中缺失部分特征,请仔细核对')
#     print(compare[,2][! compare[,2] %in% rownames(data2)])
#     q()
# }
if( sum(keep1_name %in% colnames(data1)) != length(keep1_name) )
{
    message('[Error] 数据集1 中缺失部分keep1_name,请仔细核对')
    print(keep1_name[! keep1_name %in% colnames(data1)])
    q()
}
if( sum(keep2_name %in% colnames(data2)) != length(keep2_name) )
{
    message('[Error] 数据集2 中缺失部分keep2_name,请仔细核对')
    print(keep2_name[! keep2_name %in% colnames(data2)])
    q()
}

# 开始计算相关性
message("start correlation")
result = matrix(nrow = nrow(compare), ncol = (length(keep1_name) + length(keep2_name) + 2) );
colnames(result) = c(keep1_name, keep2_name, 'estimate', 'pvalue')

# 数据存在与否
is_exists = cbind(compare[,1] %in% rownames(data1), compare[,2] %in% rownames(data2))

for(row in 1:nrow(compare))
{   
    feature1 = compare[row, 1]
    feature2 = compare[row, 2]

    keep1_data = rep('NA', length(keep1_name))
    keep2_data = rep('NA', length(keep2_name))
    pvalue     = 'NA'
    estimate   = 'NA'
    
    if(is_exists[row, 1] == TRUE & is_exists[row, 2] == TRUE)
    {
        x = as.numeric(data1[feature1, samples])
        y = as.numeric(data2[feature2, samples])
 
        tmp <- cor.test(x, 
                    y, 
                    alternative = "two.sided", 
                    method = method, 
                    conf.level = 0.95)
        pvalue <- tmp$p.value
        estimate <- tmp$estimate[[1]]
    }
    if(is_exists[row, 1] == TRUE) keep1_data = unlist(data1[feature1, keep1_name, drop = TRUE])
    if(is_exists[row, 2] == TRUE) keep2_data = unlist(data2[feature2, keep2_name, drop = TRUE])

    if(is_exists[row, 1] == FALSE) keep1_data[1] = feature1 
    if(is_exists[row, 2] == FALSE) keep2_data[1] = feature2

    result[row, ] = c(keep1_data, keep2_data, estimate, pvalue)
}

write.table(result, output , sep = "\t", quote = F, row.names = F ,col.names = T)

 


