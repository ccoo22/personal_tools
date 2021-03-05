#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: deseq2.r  -i <file> -o <dir> --case_group_name <string> --control_group_name <string> --case_sample_list <string> --control_sample_list <string> [--cor_file <file> --anno_col <string> --pvalue_cutoff <numeric> --log2fc_cutoff <numeric> --Rlib <dir>]
Options:
    -i, --input <file>              输入文件，样本reads count原始表达矩阵。 第一列是基因、转录本的ID，后面的所有列可以是样本，也可以是注释。 样本的表达量不能有缺失，缺失应该设为0。矩阵分隔符是'\\t'
    -o, --output <dir>              结果输出目录
    --case_group_name <string>      case分组名称
    --control_group_name <string>   control分组名称
    --case_sample_list <string>     case组样本列表， 样本之间用逗号分隔, 例如 a,b,c。 务必保证这些样本在input文件中存在。 注意：case + control的样本数量 要>=3。如果只有2个样本做对比，只能用老版本的deseq2或者edger的pair模式edger_two_sample_compare.r。 
    --control_sample_list <string>  control组样本列表， 样本之间用逗号分隔
    --pvalue_cutoff <numeric>       标记显著性符号时，pvalue阈值。 [default: 0.05]
    --log2fc_cutoff <numeric>       标记显著性符号时，log2 foldchange阈值, 绝对值 [default: 1]
    --anno_col <string>             input列名。从input文件中挑选指定列，放到差异结果尾部。多个列用逗号分隔。
    --cor_file <file>               协变量矫正。至少包含两列，有表头，不能有缺失值。第一列样本名，第二列及之后的列：要矫正的变量，可以是字符型、离散型数值、连续型数值。样本顺序没有限制。
    --Rlib <dir>                    R包路径  [default: /home/genesky/software/r/3.5.1/lib64/R/library] " -> doc

opts                     <- docopt(doc, version = '对两组样本的表达量做差异分析 \n          甘斌 129\n')
input                    <- opts$input
output_dir               <- opts$output
case_group_name          <- opts$case_group_name
control_group_name       <- opts$control_group_name
case_sample_list         <- opts$case_sample_list
control_sample_list      <- opts$control_sample_list
pvalue_cutoff            <- as.numeric(opts$pvalue_cutoff)
log2fc_cutoff            <- as.numeric(opts$log2fc_cutoff)
anno_col                 <- opts$anno_col
cor_file                 <- opts$cor_file
Rlib                     <- opts$Rlib
.libPaths(Rlib)

# 部分信息提前准备
case_samples = unlist(strsplit(case_sample_list, split=','))
control_samples = unlist(strsplit(control_sample_list, split=','))
all_samples = c(case_samples, control_samples)  # 顺序不要变
anno_columns = c()
if(! is.null(anno_col) ) anno_columns = unlist(strsplit(anno_col, split=','))


# input='gene_count_matrix.txt'
# output_dir='./'
# case_group_name = 'R'
# control_group_name       <- 'C'
# case_sample_list         <- 'R11H,R14H,R15H,R1H,R3H'
# control_sample_list      <- 'C1H,C2H,C3H,C4H,C5H'
# pvalue_cutoff            <- 0.05
# log2fc_cutoff            <- 1

# 加载R包
library(gplots)
library(ggplot2)
library(DESeq2)

message("load reads count matrix")
data_input <- read.table(input, header=T, sep="\t", row.names=1, comment.char="", check.names=F)

# 矫正数据
data_cor_input = data.frame()
if(! is.null(cor_file)) data_cor_input <- read.table(cor_file, header=T, sep="\t", row.names=1, comment.char="", check.names=F)

# 检查是否有丢失样本
if(sum(!all_samples %in% colnames(data_input)) > 0 )
{   
    lost_samples = all_samples[!all_samples %in% colnames(data_input)]
    message("[Error] input文件中的样本名与 case_sample_list control_sample_list 中的样本名不一致！ 丢失的样本为：", paste(lost_samples, collapse =','))
    q()
}
if(! is.null(cor_file) & sum(!all_samples %in% rownames(data_cor_input)) > 0 )
{   
    lost_samples = all_samples[!all_samples %in% rownames(data_cor_input)]
    message("[Error] cor_file文件中的样本名与 case_sample_list control_sample_list 中的样本名不一致！ 丢失的样本为：", paste(lost_samples, collapse =','))
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
data_cor_clean = data.frame()
if(! is.null(cor_file)) 
{   
    data_cor_clean = data_cor_input[all_samples, , drop=FALSE]
    if(sum(!complete.cases(data_cor_clean)))
    {
        message("[Error] cor_file文件中含有缺失值")
        q()
    }
}

# 构建差异分析模型
message("deseq2 analysis")

group = factor(c(rep(case_group_name, length(case_samples)), rep(control_group_name, length(control_samples))  ), levels = c(control_group_name, case_group_name))  # levles， control要放在前面，保证后续差异分析时，以control为对照
names(group) = all_samples
colData = data.frame(condition = group)
deseq_formula = as.formula('~ condition')

if(! is.null(cor_file))
{   
    colData = data.frame(data_cor_clean, condition = group)  # 注意，一定要把condition放在最后面
    plus_mode = paste( colnames(colData), collapse = ' + ')
    deseq_formula = as.formula( paste('~', plus_mode, sep = ' ')) 
    message("batch correction mode: ", deseq_formula)
}

# 开始deseq2处理
dds <- DESeqDataSetFromMatrix(countData = data_clean, colData = colData, design = deseq_formula )  # 初始化
dds <- DESeq(dds, betaPrior=FALSE)  # 差异分析


message("deseq2 diff output")

# 提取矫正后的表达量
cnt_norm <- as.data.frame(counts(dds, normalized=TRUE))
cnt_norm = data.frame(ID=rownames(cnt_norm), cnt_norm, stringsAsFactors=F, check.names = F)

# norm 结果输出
norm_file <- paste0(output_dir, "/diff_norm.txt")
write.table(cnt_norm, norm_file, sep="\t", quote=F, col.names = T, row.names = F)

# 差异分析结果提取
# (1) 补充case/control表达量均值
cnt_norm$baseMeanA <- apply( cnt_norm[, case_samples, drop=F], 1, mean )
cnt_norm$baseMeanB <- apply( cnt_norm[, control_samples, drop=F], 1, mean )

# (2) 差异结果提取、标记Up/Down/Not DEG
res <- as.data.frame(results(dds)) # baseMean log2FoldChange  lfcSE   stat    pvalue  padj 
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
  labs(x=case_group_name, y=control_group_name)
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
  labs(x = bquote(paste(log[2],"(fold change)",sep="")), y = bquote(paste(-log[10],"(p value)",sep="")), title = paste("valcano(",case_group_name,"/",control_group_name,")",sep="")) +
  scale_x_continuous(limits=c(-10,10)) +
  annotate("rect", xmin = log2fc_cutoff, xmax = Inf, ymin = -log10(0.05), ymax = Inf, alpha=0, colour="black") +
  annotate("text", x = (log2fc_cutoff+10)/2, y = Inf, label = sum(res$type=="Up"), color = "#f8766d", vjust = 2) +
  annotate("text", x = (log2fc_cutoff-10)/2, y = Inf, label = sum(res$type=="Down"), color = "#00ba38", vjust = 2)
dev.off()

# pvalue 排序
res      <- res[order(res$pvalue), ]
exp_diff <- data.matrix(res[res$type != "Not DEG", all_samples])
exp_all  <- data.matrix(res[, all_samples])

# 绘图 heatmap  
if(nrow(exp_diff) > 1)
{
    # 热图样本颜色
    heatmap_sample_color <- as.character(group)
    heatmap_sample_color[heatmap_sample_color == case_group_name]    <- 'red'
    heatmap_sample_color[heatmap_sample_color == control_group_name] <- 'blue'

    message("heatmap plot")
    myheatcol   = colorpanel(75, 'green','black','red')
    heatmap_file <- paste0(output_dir, "/heatmap.pdf")
    pdf(heatmap_file , width=12, height = 12 )
    heatmap.2(
      exp_diff,
      dendrogram = 'row', 
      Colv         = FALSE,
      col          = myheatcol, 
      ColSideColors = heatmap_sample_color,
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
    dev.off()

    # 绘图 Top50 heatmap 
    message("heatmap top 50 plot")
    top = 50
    if(nrow(exp_diff) < 50) top = nrow(exp_diff)
    exp_diff_top = exp_diff[1:top, ]
    
    myheatcol   = colorpanel(75, 'green','black','red')
    heatmap_file <- paste0(output_dir, "/heatmap_top50.pdf")
    pdf(heatmap_file , 12, height = 12 )
    heatmap.2(
      exp_diff_top,
      dendrogram = 'row',
      Colv         = FALSE,
      col          = myheatcol,
      ColSideColors = heatmap_sample_color,
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
    dev.off()
}else
{
    message("heatmap skipped. 差异表达量不足1个")
}

# PCA 分析（基于差异的表达量）
message("PCA plot based on diff rna")
apca <- prcomp(t(exp_diff), scale = TRUE)  # 对数据先做scale标准化处理，然后执行SVD
pc1 <- apca$x[ ,1]
pc2 <- apca$x[ ,2]
pc3 <- apca$x[ ,3]
sites <- data.frame(sample=names(pc1), PC1=pc1, PC2=pc2, PC3=pc3, Group=group[names(pc1)])
write.table(sites, paste0(output_dir, "/pc_sites.diff_rna.txt"), sep="\t", quote=F, col.names = T, row.names = F)

# 获取两个主成分贡献值
pca_importance = summary(apca)$importance
pc1_proportion <- pca_importance[2, 1]*100
pc2_proportion <- pca_importance[2, 2]*100
pca_importance = cbind(type=rownames(pca_importance),pca_importance)
write.table(pca_importance, paste0(output_dir, "/pc_cont.diff_rna.txt"), sep="\t", quote=F, col.names = T, row.names = F)


# 颜色模版、形状模版
mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
mycol <- colors()[rep(mycol, 50)]
myshape <- rep(c(15,16,17,18,19,20,21,22,23,24,25,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14),4)

# 绘图
pdf(paste(output_dir, "/pca.diff_rna.pdf", sep=""),width = 10, height = 10)
p = ggplot(data=sites, aes(x = PC1, y = PC2, colour = Group, shape = Group)) +
        geom_hline(yintercept = 0, colour = "gray65") +
        geom_vline(xintercept = 0, colour = "gray65") +
        scale_shape_manual(values = myshape ) +
        scale_color_manual(values = mycol ) +
        geom_point(size = 3, alpha = 1) +
        ggtitle("PCA plot of diff RNA") + xlab(paste("PC1 ", "(", pc1_proportion, "%)", sep="")) + ylab(paste("PC2 ", "(", pc2_proportion, "%)", sep="")) + theme(plot.title = element_text(hjust = 0.5))
p
dev.off()


# PCA 分析（基于所有的表达量）
message("PCA plot based on all rna")
apca <- prcomp(t(exp_all), scale = TRUE)
pc1 <- apca$x[ ,1]
pc2 <- apca$x[ ,2]
pc3 <- apca$x[ ,3]
sites <- data.frame(sample=names(pc1), PC1=pc1, PC2=pc2, PC3=pc3, Group=group[names(pc1)])
write.table(sites, paste0(output_dir, "/pc_sites.all_rna.txt"), sep="\t", quote=F, col.names = T, row.names = F)

# 获取两个主成分贡献值
pca_importance = summary(apca)$importance
pc1_proportion <- pca_importance[2, 1]*100
pc2_proportion <- pca_importance[2, 2]*100
pca_importance = cbind(type=rownames(pca_importance),pca_importance)
write.table(pca_importance, paste0(output_dir, "/pc_cont.all_rna.txt"), sep="\t", quote=F, col.names = T, row.names = F)


# 颜色模版、形状模版
mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
mycol <- colors()[rep(mycol, 50)]
myshape <- rep(c(15,16,17,18,19,20,21,22,23,24,25,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14),4)

# 绘图
pdf(paste(output_dir, "/pca.all_rna.pdf", sep=""),width = 10, height = 10)
p = ggplot(data=sites, aes(x = PC1, y = PC2, colour = Group, shape = Group)) +
        geom_hline(yintercept = 0, colour = "gray65") +
        geom_vline(xintercept = 0, colour = "gray65") +
        scale_shape_manual(values = myshape ) +
        scale_color_manual(values = mycol ) +
        geom_point(size = 3, alpha = 1) +
        ggtitle("PCA plot of all RNA") + xlab(paste("PC1 ", "(", pc1_proportion, "%)", sep="")) + ylab(paste("PC2 ", "(", pc2_proportion, "%)", sep="")) + theme(plot.title = element_text(hjust = 0.5))
p
dev.off()
