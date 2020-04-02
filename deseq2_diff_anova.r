#!/home/genesky/software/r/3.5.1/bin/Rscript

.libPaths("/home/genesky/software/r/3.5.1/lib64/R/library/")
library(docopt)
"Usage: deseq2_diff_anova.R --expression <file> --group_file <file> --output_prefix <file> 
Options:
   --expression <file>          reads count 表达量数据，行是基因/转录本，列是样本，第一列必须是基因/转录本名称
   --group_file <file>          样本分组文件，两列，没有表头。第一列样本名，第二列样本分组
   --output_prefix <file>       输出结果的文件前缀" -> doc

opts                <- docopt(doc, version = 'Program : deseq2_diff_anova.R v1.0 \n          甘斌 129\n')
expression          <- opts$expression
group_file          <- opts$group_file
output_prefix       <- opts$output_prefix
 
library(DESeq2)
library(gplots)
library(tidyr) # gather 函数要用到
 
countdata = read.table(expression, header = T, sep = "\t", row.names = 1, check.names = F, comment.char = "", quote = "", fill = T)

SampleInfo = read.table(group_file, header = F, sep="\t")
rownames(SampleInfo) = SampleInfo[,1]
colnames(SampleInfo) = c("sample","Group")

# sort by groups
SampleInfo = SampleInfo[order(SampleInfo$Group),]
Groups = levels(SampleInfo$Group)

# prepare data
count = countdata[,rownames(SampleInfo)]
count = count[rowSums(count) != 0 , ]
# re.count = apply(count, 2, function(d) d/sum(as.numeric(d)))
# #筛选相对丰度大于0.001的gene
# index = which(rowSums(re.count) > 0.001)

#标准化
dds <- DESeqDataSetFromMatrix(countData = count, colData = SampleInfo, design = ~ Group)
dds <- DESeq(dds,betaPrior=FALSE)
cnt <- as.data.frame(counts(dds, normalized=TRUE))
sample_num <- ncol(cnt)

# 表达量数据
DataForDEA = as.matrix(cnt)


# 计算组间比较模式
combine_count = length(Groups) * (length(Groups) - 1) / 2 # 分组组合数量
# anova给出的分组组合方式， 后面用于表头
combine_pair = t(combn(Groups, 2))
combine_pair = paste(combine_pair[,2], combine_pair[,1], sep = '-') 

# P值、分组两两之间结果存储矩阵
result <- matrix(, nrow=nrow(DataForDEA), ncol = 1 + combine_count * 4) # 保留表头、P值/FDR、数据总体统计量、分组两两之间统计量

# 设置表头
col_names = c('Pvalue');
pair_result_names = unlist(lapply(combine_pair, function(x){ paste(x, c('diff', 'lwr', 'upr', 'p.adj'), sep = '_')}  ))  # 两两组间差异

col_names = c(col_names, pair_result_names);
colnames(result) = col_names

######### 甘斌
message("start anova analysis")
for(row in 1:nrow(DataForDEA))
{   
    x <- SampleInfo$Group
    y <- DataForDEA[row,]
    analysis_data <- data.frame(x=x, y=y)
    analysis_data <- analysis_data[complete.cases(analysis_data),]
    x <- analysis_data$x
    y <- analysis_data$y
 
    # 去掉缺失值后，只有1组数据了，不能计算
    if(length(unique(x)) <= 1)
    {
    	pvalue = NA
    }else
    {
    	fit <- aov(y~x)
    	result_all       <- summary(fit) # 提取p
    	result_each_pair <- TukeyHSD(fit) # 提取组组之间的分析结果
    	pvalue <- result_all[[1]]['x', 'Pr(>F)']	

    	result_each_pair = as.data.frame(result_each_pair$x) # 转化为数据库
    	result_each_pair = data.frame(pair = row.names(result_each_pair), result_each_pair)
    	result_each_pair = gather(result_each_pair, type, statistic, -pair) # 矩阵转换形式，以TukeyHSD原始数据行、列名称重新给出数据
 
    	# 按照列名填充数据，原因：去掉缺失样本后，分组的数量可能会减少。故通过这种方式更具有通用性。
    	result[row, paste(result_each_pair[, 1], result_each_pair[, 2], sep = '_')] = result_each_pair[, 3] 
    }
    result[row, 1] <-  pvalue  # 前面半部分结果
}

############

fdr = p.adjust(result[,1], method = "BH")
means = sapply(Groups, function(g) apply(DataForDEA[,as.character(SampleInfo[as.character(SampleInfo$Group) %in% g, 1])], 1, mean))
colnames(means) = paste("Mean.In.", Groups, sep = "")
	
result  = data.frame(Gene = rownames(DataForDEA), DataForDEA, means, Pvalue=result[,1], FDR = fdr , Type = "Not DEG", result[, 2:ncol(result)], check.names = F, stringsAsFactors =F)
result  = result[order(result$Pvalue),]
result$Type[result$Pvalue < 0.05] = "DEG"


anova_file = paste(output_prefix, ".ANOVA_Test_Result.xls", sep="")
message("output : ", anova_file)
write.table(result, anova_file, row.names = F, quote = F,sep = "\t")

# heatmap

myheatcol    = colorpanel(75, 'green','black','red')
data         = data.matrix(result[result$Type == "DEG", rownames(SampleInfo)])
heatmap_file = paste(output_prefix, ".heatmap.pdf", sep="")
message("output : ", heatmap_file)
pdf(heatmap_file , width = 12, height = 12)
heatmap.2( data,
           dendrogram   = 'row', 
           Colv         = FALSE,
           col          = myheatcol, 
           scale        = "row", 
           margins      = c(5,10),
           density.info = "none",
           trace        = "none",
           key          = TRUE,
           keysize      = 0.6,
           cexRow       = 0.01,
           cexCol       = 0.8,
           srtCol       = 90
)
dev.off()

#Top50
res = result[result$Type == 'DEG', rownames(SampleInfo)]
if( nrow(res) < 50 )
{
    data = data.matrix(res)          # 小于50,全画
}else{
    data = data.matrix(res[1:50, ])       #  相对丰度最高的50个
}

heatmap_file = paste(output_prefix, ".heatmap.top50.pdf", sep="")
message("output : ", heatmap_file)
pdf(heatmap_file , width = 12, height = 12)
heatmap.2( data,
           dendrogram   = 'both', 
           #Colv         = FALSE,
           col          = myheatcol, 
           scale        = "row", 
           margins      = c(5,10),
           density.info = "none",
           trace        = "none",
           key          = TRUE,
           keysize      = 0.6,
           cexRow       = 0.8,
           cexCol       = 0.8,
           srtCol       = 90
)
dev.off()
