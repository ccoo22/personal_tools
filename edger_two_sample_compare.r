#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: edger_two_sample_compare.r  -i <file> -k <file> -o <dir> --case_sample <string> --control_sample <string> [--anno_col <string> --pvalue_cutoff <numeric> --log2fc_cutoff <numeric> --Rlib <dir>]
Options:
    -i, --input <file>              输入文件，样本reads count原始表达矩阵。 第一列是基因、转录本的ID，后面的所有列可以是样本，也可以是注释。 样本的表达量不能有缺失，缺失应该设为0。矩阵分隔符是'\\t'
    -k, --keep_house_gene <file>    管家基因列表，一列数据，没有表头。 对两个样本做差异分析时，必须以管家基因做参照计算离散值。这些基因名称必须存在于input文件中
    -o, --output <dir>              结果输出目录
    --case_sample <string>          case 样本名
    --control_sample <string>       control 样本名
    --pvalue_cutoff <numeric>       标记显著性符号时，pvalue阈值。 [default: 0.05]
    --log2fc_cutoff <numeric>       标记显著性符号时，log2 foldchange阈值, 绝对值 [default: 1]
    --anno_col <string>             input列名。从input文件中挑选指定列，放到差异结果尾部。多个列用逗号分隔。
    --Rlib <dir>                    R包路径  [default: /home/genesky/software/r/3.5.1/lib64/R/library] " -> doc

opts                     <- docopt(doc, version = '对两个样本的表达量做差异分析 \n          甘斌 129\n')
input                    <- opts$input
keep_house_gene          <- opts$keep_house_gene
output_dir               <- opts$output
case_sample              <- opts$case_sample
control_sample           <- opts$control_sample
pvalue_cutoff            <- as.numeric(opts$pvalue_cutoff)
log2fc_cutoff            <- as.numeric(opts$log2fc_cutoff)
anno_col                 <- opts$anno_col
Rlib                     <- opts$Rlib
.libPaths(Rlib)

# 部分信息提前准备
all_samples = c(case_sample, control_sample)
anno_columns = c()
if(! is.null(anno_col) ) anno_columns = unlist(strsplit(anno_col, split=','))

# 加载R包
library(gplots)
library(ggplot2)
library(edgeR)

message("load reads count matrix")
data_input <- read.table(input, header=T, sep="\t", row.names=1, comment.char="", check.names=F)
data_keep_house <- read.table(keep_house_gene, header=F, sep="\t" , comment.char="", check.names=F, stringsAsFactors=F)[, 1]

# 检查是否有丢失样本
if(sum(!all_samples %in% colnames(data_input)) > 0 )
{   
    lost_samples = all_samples[!all_samples %in% colnames(data_input)]
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
data_raw = data_input[, all_samples]
data_clean = data_raw[apply(data_raw, 1, sum) > 0 , ]
housekeeping = data_keep_house[data_keep_house %in% rownames(data_clean)]  # 管家基因列表
if(length(housekeeping) == 0)
{
    message("[Error] 在data_clean 数据中，不存在管家基因")
    q()
}

# 构建差异分析模型
message("edger analysis")

group = factor(c(case_sample, control_sample), levels = c(control_sample, case_sample))  # levles， control要放在前面，保证后续差异分析时，以control为对照
names(group) = all_samples
colData = data.frame(condition = group)
deseq_formula = as.formula('~ condition')

# 开始edgeR处理
message("edgeR diff output")
y <- DGEList(counts=data_clean, group=group)  # 注意：此处的group不能省
y <- calcNormFactors(y)  # 默认的矫正方法是TMM

# 基于管家基因 分析dispersion 
y1 <- y
y1$samples$group <- 1
# Then estimate the common dispersion from the housekeeping genes and all the libraries
# as one group:
y0 <- estimateDisp(y1[housekeeping,], trend="none", tagwise=FALSE)
# Then insert this into the full data object and proceed:
y$common.dispersion <- y0$common.dispersion

# 差异
design <- model.matrix(deseq_formula, data=colData)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)

# 差异结果提取
tTag = topTags(lrt, n=nrow(lrt)) # logFC   logCPM       LR       PValue          FDR。 其中logCPM是矫正后的表达量取均值，然后取log2
res <- as.data.frame(tTag)  # 差异分析结果 
colnames(res) = c('log2FoldChange', 'baseMean', 'LR', 'pvalue', 'padj')  # 修改一下列名，与deseq2保持一致
res$baseMean = 2^res$baseMean 
res = res[, c('baseMean', 'log2FoldChange', 'LR', 'pvalue', 'padj')] # 调整一下顺序， 与deseq2一致

# 提取矫正后的表达量
cnt_norm <- as.data.frame(cpm(y))
cnt_norm = data.frame(ID=rownames(cnt_norm), cnt_norm, stringsAsFactors=F, check.names = F)

# norm 结果输出
norm_file <- paste0(output_dir, "/diff_norm.txt")
write.table(cnt_norm, norm_file, sep="\t", quote=F, col.names = T, row.names = F)

# 差异分析结果提取
# (1) 补充case/control表达量均值
cnt_norm$baseMeanA <- cnt_norm[, case_sample,] 
cnt_norm$baseMeanB <- cnt_norm[, control_sample] 

# (2) 差异结果提取、标记Up/Down/Not DEG
res$type <- "Not DEG"
res$type[res$pvalue < pvalue_cutoff & res$log2FoldChange >= log2fc_cutoff ]   <- "Up"
res$type[res$pvalue < pvalue_cutoff & res$log2FoldChange <= -(log2fc_cutoff)] <- "Down"
res$type <- factor(res$type, levels = c("Up", "Down", "Not DEG"))
res = cbind(cnt_norm[rownames(res), ], res)
# 补充注释信息
if(length(anno_columns) > 0) res = cbind(res, data_input[res$ID, anno_columns, drop=FALSE] )

# (3) 结果输出
diff_file <- paste0(output_dir, "/diff.txt")
write.table(res, diff_file, sep="\t", quote=F, col.names = T, row.names = F)


# 绘图 correlation plot
message("correlation plot")
correlation_file <- paste0(output_dir, "/correlation.pdf")
pdf(correlation_file)
ggplot(res, aes(x=log10(baseMeanA+ 0.00000000001),y= log10(baseMeanB+0.00000000001))) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color="red") + 
  scale_x_continuous(limits=c(-2.5, 5)) + 
  scale_y_continuous(limits=c(-2.5, 5)) + 
  labs(x=case_sample, y=control_sample)
dev.off()


# 绘图 MA  
message("MA plot")
ma_file <- paste0(output_dir, "/ma.pdf")
pdf(ma_file)
ggplot(res, aes(x = baseMean, y = log2FoldChange , colour = type)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 2e+04)) + 
  scale_y_continuous(limits = c(-20, 20)) + 
  theme(legend.title = element_blank()) + 
  labs(x = "Mean Expression", y="log2(FC)", tilte="MA plot")
dev.off()


# 绘图 valcano  
message("valcano plot")
valcano_file <- paste0(output_dir, "/valcano.pdf")
pdf(valcano_file)
ggplot(res, aes(x = log2FoldChange, y = -log10(res$pvalue), color =type)) + 
  geom_point() +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = bquote(paste(log[2],"(fold change)",sep="")), y = bquote(paste(-log[10],"(p value)",sep="")), title = paste("valcano(",case_sample,"/",control_sample,")",sep="")) +
  scale_x_continuous(limits=c(-10,10)) +
  annotate("rect", xmin = log2fc_cutoff, xmax = Inf, ymin = -log10(0.05), ymax = Inf, alpha=0, colour="black") +
  annotate("text", x = (log2fc_cutoff+10)/2, y = Inf, label = sum(res$type=="Up"), color = "#f8766d", vjust = 2) +
  annotate("text", x = (log2fc_cutoff-10)/2, y = Inf, label = sum(res$type=="Down"), color = "#00ba38", vjust = 2)
dev.off()
