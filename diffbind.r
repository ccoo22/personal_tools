#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: diffbind.r  -i <file> -o <dir> --case_group_name <string> --control_group_name <string> [--cor_file <file> --pvalue_cutoff <numeric> --log2fc_cutoff <numeric> --Rlib <dir>]
Options:
    -i, --input <file>              diffbind规定的sheet.csv文件。包含以下11个表头：SampleID	Tissue	Factor	Condition	Treatment	Replicate	bamReads	ControlID	bamControl	Peaks	PeakCaller
                                    SampleID: 样本名称。case分组样本统一放在前面
                                    Tissue: 填NA即可
                                    Factor: 填写样本分组名称。
                                    Condition:NA
                                    Treatment:NA
                                    Replicate:样本编号，case组样本按顺序填写1、2、3即可。control组样本也是1、2、3
                                    bamReads：样本的原始bam文件（不要经过任何处理，bowtie2比对后、排序的bam）
                                    ControlID：input样本ID。如果没有，填写NA
                                    bamControl:input样本bam。如果没有，填写NA
                                    Peaks: macs2 分析得到的narrowPeak文件
                                    PeakCaller: 统一填写narrow
    -o, --output <dir>              结果输出目录
    --case_group_name <string>      case分组名称
    --control_group_name <string>   control分组名称
    --pvalue_cutoff <numeric>       标记显著性符号时，pvalue阈值。 [default: 0.05]
    --log2fc_cutoff <numeric>       标记显著性符号时，log2 foldchange阈值, 绝对值 [default: 1]
    --cor_file <file>               协变量矫正。至少包含两列，有表头，不能有缺失值。第一列样本名，第二列及之后的列：要矫正的变量，可以是字符型、离散型数值、连续型数值。样本顺序没有限制。
    --Rlib <dir>                    R包路径  [default: /home/genesky/software/r/3.5.1/lib64/R/library] " -> doc

opts                     <- docopt(doc, version = '基于diffbind + deseq2，做ATAC的两组之间的差异分析 \n          甘斌 129\n')
input                    <- opts$input
output_dir               <- opts$output
case_group_name          <- opts$case_group_name
control_group_name       <- opts$control_group_name
pvalue_cutoff            <- as.numeric(opts$pvalue_cutoff)
log2fc_cutoff            <- as.numeric(opts$log2fc_cutoff)
cor_file                 <- opts$cor_file
Rlib                     <- opts$Rlib
.libPaths(Rlib)



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
library(DiffBind)

# （1）读入sheet，并对peak区域做交集、合并处理
message("read samplesheet")
db_obj <- dba(sampleSheet=input) # db_obj$binding 这个数据就是合并后的区域

case_samples = as.character(db_obj$samples$SampleID[db_obj$samples$Factor == case_group_name])  # case分组样本名称
control_samples = as.character(db_obj$samples$SampleID[db_obj$samples$Factor == control_group_name]) # control分组样本名称
all_samples = c(case_samples, control_samples)
if(length(case_samples) == 0 || length(control_samples) == 0)
{
    message("[Error] 在sheet文件中，case_group_name 或 ：control_group_name 没有匹配到, 请仔细检查:", paste(c(case_group_name, control_group_name), collapse =','))
    q()
}


# 矫正数据
data_cor_input = data.frame()
if(! is.null(cor_file)) data_cor_input <- read.table(cor_file, header=T, sep="\t", row.names=1, comment.char="", check.names=F)
if(! is.null(cor_file) & sum(!all_samples %in% rownames(data_cor_input)) > 0 )
{   
    lost_samples = all_samples[!all_samples %in% rownames(data_cor_input)]
    message("[Error] cor_file文件中的样本名与 case_sample_list control_sample_list 中的样本名不一致！ 丢失的样本为：", paste(lost_samples, collapse =','))
    q()
}

# （2）计算每个peaks/regions的count信息
message("read bam to get reads count")
db_obj <- dba.count(db_obj, bUseSummarizeOverlaps=TRUE)

# 输出reads count数据 
message("output raw reads count")
reads_count = sapply(db_obj$peaks, function(x){x$Reads})
reads_count[reads_count < 1] = 1  # 把0赋值为1  # diffbind规则
reads_count = data.frame(id = paste('peak_', 1:db_obj$totalMerged, sep=''), db_obj$peaks[[1]][, c('Chr', 'Start', 'End')], reads_count, stringsAsFactors=F)
colnames(reads_count) = c('id', 'Chr', 'Start', 'End', db_obj$samples$SampleID)
rownames(reads_count) = reads_count$id
reads_count_file <- paste0(output_dir, "/reads_count.txt")
write.table(reads_count, reads_count_file, sep="\t", quote=F, col.names = T, row.names = F)

# 输出bam reads count数据
bam_reads = as.numeric(db_obj$class['Reads', ])
names(bam_reads) = db_obj$samples$SampleID
reads_count_bam_file <- paste0(output_dir, "/reads_count_bam.txt")
write.table(bam_reads, reads_count_bam_file, sep="\t", quote=F, col.names = F, row.names = T)


# 由于diffbind不支持批次效应矫正，故放弃他们的差异分析过程，直接用他们提取的reads count数据 + deseq2计算过程，我们自己算差异
# # （3） diffbind 差异分析
# db_obj <- dba.contrast(db_obj, categories=DBA_FACTOR,minMembers = 2)  # 基于sheet文件的Factor列建立分组对比模型
# db_obj <- dba.analyze(db_obj, method=DBA_ALL_METHODS, bSubControl=F)  # 开始差异分析。 method：差异分析的算法选择，有DBA_EDGER、DBA_DESEQ2、DBA_DESEQ三种。DBA_ALL_METHODS表示同时做 DBA_EDGER和DBA_DESEQ2. 另外，因为没有input样本，所以吧bSubControl设为F

# # (4) diffbind 差异分析 结果提取
# # edger
# comp1.edgeR <- dba.report(db_obj, method=DBA_EDGER, contrast = 1, th=1)
# comp1.edgeR_count <- dba.report(db_obj, method=DBA_EDGER, contrast = 1, th=1, bCounts = TRUE)  
# # comp1.edgeR_count <- dba.report(db_obj, method=DBA_EDGER, contrast = 1, th=1, bCounts = TRUE, bNormalized=FALSE)  # 返回原始count
# meta_colnames = colnames(comp1.edgeR_count@elementMetadata)  # 提取原始meta信息的列名称,包括样本名，一些统计名称等。
# out_count = as.data.frame(comp1.edgeR_count, stringsAsFactors=F)
# colnames(out_count)[6:ncol(out_count)] = meta_colnames

# # deseq2
# # comp1.deseq2 <- dba.report(db_obj, method=DBA_DESEQ2, contrast = 1, th=1)
# # comp1.deseq2_count <- dba.report(db_obj, method=DBA_DESEQ2, contrast = 1, th=1, bCounts = TRUE)  



# 甘斌：解析原始代码，自己计算差异分析，因为需要导入环境变量做矫正

# 提取需要的样本、去掉表达量为0的基因
data_raw = reads_count[, all_samples]
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

# # 开始deseq2处理
libsize = as.numeric(bam_reads[all_samples])
dds <- DESeq2::DESeqDataSetFromMatrix(countData = data_clean, colData = colData, design = deseq_formula )  # 初始化
DESeq2::sizeFactors(dds) = libsize/min(libsize)
dds <- DESeq2::estimateDispersions(dds,fitType='local')
dds <- DESeq2::nbinomWaldTest(dds)


message("deseq2 diff output")
# 提取矫正后的表达量
cnt_norm <- as.data.frame(DESeq2::counts(dds, normalized=TRUE))
cnt_norm = data.frame(ID=rownames(cnt_norm), cnt_norm, stringsAsFactors=F, check.names=F)

# norm 结果输出
norm_file <- paste0(output_dir, "/diff_norm.txt")
write.table(cnt_norm, norm_file, sep="\t", quote=F, col.names = T, row.names = F)

# 差异分析结果提取
# (1) 补充case/control表达量均值
cnt_norm$baseMeanA <- apply( cnt_norm[, case_samples, drop=F], 1, mean )
cnt_norm$baseMeanB <- apply( cnt_norm[, control_samples, drop=F], 1, mean )

# (2) 差异结果提取、标记Up/Down/Not DEG
res <- as.data.frame(DESeq2::results(dds)) # baseMean log2FoldChange  lfcSE   stat    pvalue  padj 
res$type <- "Not DEG"
res$type[res$pvalue < pvalue_cutoff & res$log2FoldChange >= log2fc_cutoff ]   <- "Up"
res$type[res$pvalue < pvalue_cutoff & res$log2FoldChange <= -(log2fc_cutoff)] <- "Down"
res$type <- factor(res$type, levels = c("Up", "Down", "Not DEG"))
res = cbind(cnt_norm[rownames(res), ], res)
 
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

# 热图绘制数据，最多20000条，否则占用内存太大，导致R死掉
exp_heatmap <- exp_diff
if(nrow(exp_heatmap) > 20000) exp_heatmap <- exp_diff[1:20000,]

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
      exp_heatmap,
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
    if(nrow(exp_heatmap) < 50) top = nrow(exp_heatmap)
    exp_diff_top = exp_heatmap[1:top, ]
    
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