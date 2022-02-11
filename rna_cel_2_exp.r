#!/home/genesky/software/r/4.0.3/bin/Rscript

library(docopt)
"Usage: rna_cel_2_exp.r --input <file>  --output  <file>
Options:
    --input, -i <file>   CEL信息文件，有表头（列名随意），两列数据，第一列样本名（用做表达量文件的样本名），第二列CEL文件路径
    --output, -o <file>  输出文件， 例如 ： result.txt
" -> doc

opts                   <- docopt(doc, version = '把RNA芯片原始CEL文件转换为表达量矩阵文件，主要是针对 Affymetrix Human Genome U133 平台 \n          甘斌 129\n')
input            <- opts$input
output            <- opts$output

library(affy)

data_input = read.table(input, header=T, sep='\t', check.names=F, quote='', stringsAsFactors=F, colClasses = 'character')
colnames(data_input) = c('sample', 'cel')
message('[process] 读入样本 ', nrow(data_input), ' 个')

# cel 读入
message('[process] 读入cel文件')
Data <- ReadAffy(filenames = data_input$cel, sampleNames = data_input$sample)

message('[process] 基于rma函数做 background correction、 normalization、 pm correction、 summarization')
eset <- rma(Data)

message('[process] 获取表达量、输出')
exp = exprs(eset) 
write.table(data.frame(id=rownames(exp), exp), output, sep='\t', row.names=F, quote=F)


