#!/home/genesky/software/r/3.5.1/bin/Rscript

.libPaths("/home/genesky/software/r/3.5.1/lib64/R/library/")
library(docopt)
"Usage: rna_diff_limma_chip.r --input <file> --case_sample <string> --control_sample <string> --output_prefix <file>  [--do_log2  --pvalue_cutoff <numeric> --case_group_name <string>  --control_group_name <string> --append_anno_col <string>] 
Options:
   --input <file>                  样本的表达量矩阵，第一列是基因/转录本ID, 第二列开始是样本的表达量信息。含有表头，tab分隔符。 注意：通常情况下，当前脚本仅适用于rna芯片的表达量数据分析。如果你有reads count数据，建议使用deseq2。
   --case_sample <string>          case组样本编号，样本之间用逗号分隔 例如： 'GSM2183539,GSM2183540,GSM2183541,GSM2183542,GSM2183543,GSM2183544,GSM2183545,GSM2183546,GSM2183547,GSM2183548'
   --control_sample <string>       control组样本编号，样本之间用逗号分隔 例如： 'GSM2183532,GSM2183533,GSM2183534,GSM2183535,GSM2183536,GSM2183537,GSM2183538'
   --do_log2                       是否对输入数据做log2处理， 默认不处理。 limma的分析要求必须对数据做log2转换，如果输入矩阵已经转换完毕，可以不用理会这个参数
   --case_group_name <string>      case组名 [default: case]
   --control_group_name <string>   control组名 [default: control]
   --output_prefix <file>          输出文件的前缀，例如： ./result  。 最终生成的结果为 result.diff.xls。
   --pvalue_cutoff <numeric>       pvalue 阈值，仅输出小于该阈值的结果 [default: 1]
   --append_anno_col <string>      差异结果中，添加input文件中指定列的注释内容， 例如： 2,3,4   把input文件中的2,3,4数据追加到输出文件尾部。 默认，不追加 。 第一列编号是 1    " -> doc

opts                   <- docopt(doc, version = 'Program : rna_diff_limma_chip v1.0 \n          甘斌 129\n')
input                  <- opts$input
case_sample_list       <- opts$case_sample
control_sample_list    <- opts$control_sample
output_prefix          <- opts$output_prefix
do_log2                <- opts$do_log2
pvalue_cutoff          <- as.numeric(opts$pvalue_cutoff)
output_expr            <- opts$output_expr
case_group_name        <- opts$case_group_name
control_group_name     <- opts$control_group_name
append_anno_col        <- opts$append_anno_col

case_samples    = unlist(strsplit(case_sample_list, ','))
control_samples = unlist(strsplit(control_sample_list, ','))
sample_group  <- c(rep(case_group_name, length(case_samples)), rep(control_group_name, length(control_samples)))


# https://www.ncbi.nlm.nih.gov/geo/info/geo2r.html
##################
# （1） 加载R包
##################
library(Biobase)
library(limma)
library(ggplot2)
library(gplots)


##################
# (1) 读入数据
##################
message("[process] loading data")
data_input = read.table(input, head = T, row.names = 1, check.names = F, sep = "\t", quote ="")

##################
# (2) 检查样本存在与否
##################
message("[process] check input sample ok")
if(sum(case_samples %in% colnames(data_input)) != length(case_samples))
{   
    error = case_samples[! case_samples %in% colnames(data_input)]
    message("[Error] 输入的case样本中，部分样本在输入矩阵中丢失： ", error)
    q()
}
if(sum(control_samples %in% colnames(data_input)) != length(control_samples))
{   
    error = control_samples[! control_samples %in% colnames(data_input)]
    message("[Error] 输入的control样本中，部分样本在输入矩阵中丢失： ", error)
    q()
}

##################
# (3) 提取分析的数据
##################
message("[process] build sample data matrix")
data_matrix <- data_input[, c(case_samples, control_samples)]  # 仅保留需要的样本
data_matrix <- as.matrix(data_matrix)

# log2处理
if(do_log2)
{
    data_matrix = log2(data_matrix)
    data_matrix[data_matrix == -Inf] = NA
}

# 建立数据集
gset <- ExpressionSet(assayData = data_matrix) 

##################
# (4) 差异分析
##################
message("[process] diff analysis")
sml <- paste("G", c(rep(1, length(case_samples)), rep(0, length(control_samples))), sep="")

fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="p", n=Inf, p.value=pvalue_cutoff)  # 最终差异分析结果


##################
# (5) 差异分析结果整理输出
##################
message("[process] write diff result")

result <- data.frame(probe = rownames(tT),
                 data_matrix[rownames(tT), ],
                 caseExpr = apply(data_matrix[rownames(tT), case_samples], 1, function(x){mean(x, na.rm =T)}),
                 controlExpr = apply(data_matrix[rownames(tT), control_samples], 1, function(x){mean(x, na.rm =T)}),
                 averageExpr = tT$AveExpr,
                 log2FC = tT$logFC,
                 pvalue = tT$P.Value,
                 fdr = tT$adj.P.Val
                 )
result$type <- "Not DEG"
result$type[result$pvalue < 0.05 & result$log2FC >= 1 ] <- "Up"
result$type[result$pvalue < 0.05 & result$log2FC <= -1] <- "Down"
result$type <- factor(result$type, levels = c("Up", "Down", "Not DEG"))

 
# 追加指定注释
if(!is.null(append_anno_col))
{ 
    need_col <- as.integer(unlist(strsplit(append_anno_col, ','))) - 1  # 输入文件的第一列用作行名了，故前移一位
    append_info <- data_input[rownames(tT), need_col, drop = F]
    result <- data.frame(result[, ], append_info[, ])
 
}
 
file_diff = paste0(output_prefix, ".diff.xls")
message("output diff file : ", file_diff)
write.table(result, file=file_diff, row.names=F, sep="\t", quote = FALSE)
 
 

##################
# （7）MA plot
##################
ma_file <- paste(output_prefix, ".MA.pdf", sep="")
message("[process] plot MA : ", ma_file)
pdf(ma_file)
ggplot(result, aes(x = averageExpr, y = log2FC , colour = type)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 10)) + 
  scale_y_continuous(limits = c(-20, 20)) + 
  theme(legend.title = element_blank()) + 
  labs(x = "averageExpr", y="log2(FC)", tilte="MA plot")
dev.off()

##################
# （8）valcano plot
##################
volcano_file <- paste(output_prefix, ".volcano.pdf", sep="")
message("[process] plot volcano : ", volcano_file)
pdf(volcano_file)
volcano_title = paste0("genes: ", case_group_name, control_group_name)
ggplot(result, aes(x = log2FC, y = -log10(result$pvalue), color =type)) + 
  geom_point() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = bquote(paste(log[2],"(fold change)",sep="")), y = bquote(paste(-log[10],"(p value)",sep="")), title = volcano_title) +
  scale_x_continuous(limits=c(-10,10)) 
dev.off()

##################
# （9）热图绘制
##################
result_diff <- result[result$type != 'Not DEG', ]

# Top50 heatmap plot
if( nrow(result_diff) < 50 )
{
               data <- data.matrix(result_diff[ , c(case_samples, control_samples)])          # 小于50,全画
}else{
               result_diff <- result_diff[order(result_diff$pvalue), ]
               data <- data.matrix(result_diff[1:50, c(case_samples, control_samples)])           #  相对丰度最高的50个
}
 
heatmap_file <- paste(output_prefix, ".heatmap.top50.pdf", sep="")
message("[process] plot heatmap top50 : ", heatmap_file)

pdf(heatmap_file , width=12, height = 12 )

# 样本颜色
my_plot_color <- sml
my_plot_color[my_plot_color == 'G1'] <- 'red'
my_plot_color[my_plot_color == 'G0'] <- 'blue'
myheatcol     <- colorpanel(75, 'green','black','red')
heatmap.2(
          data,
        dendrogram = 'both',
      #Colv         = FALSE,
      col          = myheatcol,
      ColSideColors = my_plot_color,
      scale="row",
      margins	   = c(8,10),
      density.info = "none",
      trace="none",
      key=TRUE,
      keysize=0.6,
	  cexRow=0.8,
      cexCol=0.8,
      srtCol=90
)


# heatmap plot
data <- data.matrix(result_diff[, c(case_samples, control_samples)])
 
heatmap_file <- paste(output_prefix, ".heatmap.pdf", sep="")
message("[process] plot heatmap all : ", heatmap_file)

pdf(heatmap_file , width=12, height = 12 )

# 样本颜色
my_plot_color <- sml
my_plot_color[my_plot_color == 'G1'] <- 'red'
my_plot_color[my_plot_color == 'G0'] <- 'blue'
myheatcol     <- colorpanel(75, 'green','black','red')
heatmap.2(
      data,
      dendrogram = 'row', 
      Colv         = FALSE,
      col          = myheatcol, 
      ColSideColors = my_plot_color,
      scale="row", 
      margins      = c(8,10),
      density.info = "none",
      trace="none",
      key=TRUE,
      keysize=0.6,
      cexRow=0.01,
  #cexRow=0.8,
      cexCol=0.8,
     srtCol=90
)
