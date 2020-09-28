#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: /deseq2_normalized.r   -i <file> -o <dir> --sample_file <file>  [ --anno_col <string> --Rlib <dir>]
Options:
    -i, --input <file>              输入文件，样本reads count原始表达矩阵。 第一列是基因、转录本的ID，后面的所有列可以是样本，也可以是注释。 样本的表达量不能有缺失，缺失应该设为0。矩阵分隔符是'\\t'
    -o, --output <dir>              输出文件
    -s, --sample_file <file>        样本列表，一列数据，有表头。
    --anno_col <string>             input列名。从input文件中挑选指定列，放到差异结果尾部。多个列用逗号分隔。
    --Rlib <dir>                    R包路径  [default: /home/genesky/software/r/3.5.1/lib64/R/library] " -> doc

opts                 <- docopt(doc, version = '使用deseq2对reads count数据做矫正，得到矫正后的表达矩阵 \n          甘斌 129\n')
input                <- opts$input
output               <- opts$output
sample_file          <- opts$sample_file
anno_col                 <- opts$anno_col
Rlib                 <- opts$Rlib
.libPaths(Rlib)

anno_columns = c()
if(! is.null(anno_col) ) anno_columns = unlist(strsplit(anno_col, split=','))

# 加载R包
library(DESeq2)

message("load reads count matrix")
data_input <- read.table(input, header=T, sep="\t", row.names=1, comment.char="", check.names=F)
data_sample <- read.table(sample_file, header=T, sep="\t", comment.char="", stringsAsFactors=F, colClasses = 'character')
samples = data_sample[,1]



# 检查是否有丢失样本
if(sum(!samples %in% colnames(data_input)) > 0 )
{   
    lost_samples = samples[!samples %in% colnames(data_input)]
    message("[Error] input文件中的样本名与 case_sample_list control_sample_list 中的样本名不一致！ 丢失的样本为：", paste(lost_samples, collapse =','))
    q()
}
if(length(anno_columns) > 0  & sum(! anno_columns %in% colnames(data_input)) > 0 )
{
    lost_columns = anno_columns[!anno_columns %in% colnames(data_input)]
    message("[Error] input文件 不存在anno_col中声明的列名，请仔细检查", paste(lost_columns, collapse =','))
    q()
}



# 提取需要的样本、去掉表达量为0的基因
data_raw = data_input[, samples]
data_clean = data_raw[apply(data_raw, 1, sum) > 0 , ]

# 构建差异分析模型
message("deseq2")


group = factor(c('a', rep('b', length(samples) - 1)  ), levels = c('a', 'b'))  # levles， control要放在前面，保证后续差异分析时，以control为对照
names(group) = samples
colData = data.frame(condition = group)
deseq_formula = as.formula('~ condition')

 
# 开始deseq2处理
dds <- DESeqDataSetFromMatrix(countData = data_clean, colData = colData, design = deseq_formula )  # 初始化
dds <- estimateSizeFactors(dds, type = "ratio", quiet = FALSE)

# 提取矫正后的表达量
cnt_norm <- as.data.frame(counts(dds, normalized=TRUE))
cnt_norm = data.frame(ID=rownames(cnt_norm), cnt_norm, stringsAsFactors=F, check.names = F)

# 补充注释信息
if(length(anno_columns) > 0) cnt_norm = cbind(cnt_norm, data_input[cnt_norm$ID, anno_columns, drop=FALSE] )

# norm 结果输出
write.table(cnt_norm, output, sep="\t", quote=F, col.names = T, row.names = F)