#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: plot_transcript_structure.r  -i <file> -o <pdf file>  [--type <string>  --rlib <string>]

Options:
    -i, --input <file>              输入文件，一列数据，没有表头，转录本ID号（例如：ENST00000216129），或者基因名称（例如：TTLL12）
    -o, --output <pdf file>         输出pdf文件路径,示例：./a.pdf
    -t, --type <string>             输入数据类型，基因/转录本，g/t [default: t]
    --rlib <string>                 r包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，绘制转录本或基因在基因组上的区域结构图\n')
input              <- opts$input
output             <- opts$output
type             <- opts$type
rlib               <- opts$rlib
 
.libPaths(rlib)
library(ggbio)
library(EnsDb.Hsapiens.v86)
ensdb <- EnsDb.Hsapiens.v86 # GRCH38 基因注释数据库

message("读入数据")
data = read.table(input, head = F,, sep = "\t", stringsAsFactors=F)


pdf(output, width = 21, height = 10)
for(name in data[, 1])
{   
    message("绘制：",name)
    pic = NULL
    if(type == 't') pic = autoplot(ensdb, TxIdFilter(c(name)) , label = TRUE )
    if(type == 'g') pic = autoplot(ensdb, GenenameFilter(c(name)) , label = TRUE, main = name )  
    print(pic) 
}
dev.off()
 
