#!/usr/bin/env Rscript

library(docopt)
"Usage: Deseq2_diff_kallisto.R [options] INPUT CASE_GROUP CONTROL_GROUP CASE CONTROL OUTPUT
Options:
   -w --width=width    the width of viewport  [default: 12]
   -h --height=width   the height of viewport [default: 12]
Arguments:
  INPUT          the input file name
  CASE_GROUP     the case group name
  CONTROL_GROUP  the control group name
  CASE           the case group sample name
  CONTROL        the control group sample name
  OUTPUT         the output directory name" -> doc


opts          <- docopt(doc)
input         <- opts$INPUT
case_group    <- opts$CASE_GROUP
control_group <- opts$CONTROL_GROUP
case          <- opts$CASE
control       <- opts$CONTROL
output        <- opts$OUTPUT


library(tximport)
library(DESeq2)
library("EnsDb.Hsapiens.v86")
library("ensembldb")

# (1) 读入kallisto数据
file_table = read.table(input, head = F, stringsAsFactors = F) # 两列，第一列：样本名，第二列：文件路径，没有表头
file_list  = file_table[[2]]
names(file_list) = file_table[[1]] # 样本名

# 数据库信息获取
db <- EnsDb.Hsapiens.v86
tx <- transcripts(db, return.type = "DataFrame")
tx <- tx [, c("tx_id", "gene_id")]

gene      <- genes(db, return.type = "DataFrame")
gene_type <- gene[, c("symbol", "gene_biotype")]
gene_type <- gene_type[!duplicated(gene_type$symbol) ,] # 安基因名，去重复
gene      <- gene[, c("gene_id", "symbol")]

tx2genes <- merge(tx, gene, by = "gene_id")
tx2genes <- tx2genes[, 2:3]

# 读入表达量,并在基因水平汇总
txi <- tximport(files = file_list, type = "kallisto", tx2gene = tx2genes, ignoreTxVersion = T) # 必须加上ignoreTxVersion参数，防止转录本版本号导致无法匹配

# 基因biotype注释
gene_anno = data.frame(symbol = row.names(txi$counts))
gene_anno = merge(gene_anno, gene_type, by = 'symbol', all.x = T, sort = F) # 获取基因的类型
colnames(gene_anno) = c('Gene', 'gene_biotype')

# 输出 基因biotype
write.table(gene_anno, file = paste0(output, '/kallisto.gene.anno_gene_biotype.txt'), sep = "\t", quote = F, row.names = FALSE)

# (2) 差异分析
case_sample = unlist(strsplit(case, split=","))
control_sample = unlist(strsplit(control, split=","))
sampletable   <- data.frame(row.names = c(control_sample, case_sample), 
                      condition = rep(
                                  c(control_group, case_group),
                                  times  = c(length(control_sample), length(case_sample)) ,
                                  levels = c(control_group, case_group)
                                 )          
                      )
sampletable = data.frame(row.names = colnames(txi$counts), condition = sampletable[colnames(txi$counts), ])
sampletable$condition <- relevel(sampletable$condition, ref = control_group)

dds  <- DESeqDataSetFromTximport(txi, sampletable, ~condition)
# PCA
rld <- rlogTransformation(dds)   #DEseq2自己的方法标准化数据`
exprSet_new <- assay(rld)   #提取DEseq2标准化后的数据
normalization_file <- paste(output, "/DESeq2.reads_count_normalization.pca.txt", sep = "")
write.table(exprSet_new, file = normalization_file, sep = "\t",quote = F)

# PCA图
pca_pdf <- paste(output, "/pca.pdf", sep = "")
pdf(pca_pdf)
plotPCA(rld, intgroup = c("condition"), returnData = FALSE)  # DEseq2自带函数
dev.off()
# PCA 坐标
pcaData <- plotPCA(rld, intgroup = c("condition"), returnData = TRUE)  # DEseq2自带函数
pcaData = data.frame(Sample = rownames(pcaData), pcaData)
pca_file <- paste(output, "/pca.xls", sep = "")
write.table(pcaData, file = pca_file, sep = "\t", quote = F, row.names = FALSE)


# 差异分析
dds <- DESeq(dds,betaPrior=FALSE)
cnt <- as.data.frame(counts(dds, normalized=TRUE))
sample_num <- ncol(cnt)
cnt$baseMeanA <- apply( sapply(unlist(strsplit(control, split=",")), function(t){cnt[[t]]}), 1, mean )
cnt$baseMeanB <- apply( sapply(unlist(strsplit(case, split=",")), function(t){cnt[[t]]}), 1, mean )
res <- as.data.frame(results(dds))
res <- cbind(cnt, res)
res$type <- "Not DEG"
res$type[res$pvalue < 0.05 & res$log2FoldChange >= 1 ] <- "Up"
res$type[res$pvalue < 0.05 & res$log2FoldChange <= -1] <- "Down"
res$type <- factor(res$type, levels = c("Up", "Down", "Not DEG"))
res = data.frame(Feature = rownames(res), res, check.names = FALSE); # 第一列添加行名。
diff_file <- paste(output, "/diff.xls", sep="")
write.table(res, diff_file, sep="\t", quote=F, row.names = FALSE)


# correlation plot
library(ggplot2)

correlation_file <- paste(output, "/correlation.pdf", sep="")
pdf(correlation_file)
ggplot(res, aes(x=log10(baseMeanA+ 0.00000000001),y= log10(baseMeanB+0.00000000001))) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color="red") + 
  scale_x_continuous(limits=c(-2.5, 5)) + 
  scale_y_continuous(limits=c(-2.5, 5)) + 
  labs(x="log10(baseMeanA)", y="log10(baseMeanB)")
dev.off()

# MA plot
ma_file <- paste(output, "/MA.pdf", sep="")
pdf(ma_file)
ggplot(res, aes(x = baseMean, y = log2FoldChange , colour = type)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 2e+04)) + 
  scale_y_continuous(limits = c(-20, 20)) + 
  theme(legend.title = element_blank()) + 
  labs(x = "baseMean", y="log2(FC)", tilte="MA plot")
dev.off()

# valcano plot
valcano_file <- paste(output, "/valcano.pdf", sep="")
pdf(valcano_file)
ggplot(res, aes(x = log2FoldChange, y = -log10(res$pvalue), color =type)) + 
  geom_point() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = bquote(paste(log[2],"(fold change)",sep="")), y = bquote(paste(-log[10],"(p value)",sep="")), title = "genes: CASE/CONTROL") +
  scale_x_continuous(limits=c(-10,10)) 
dev.off()


# Top50 heatmap plot
library(gplots)
res <- res[res$type != 'Not DEG', ]
if( nrow(res) < 50 ){
               data <- data.matrix(res[ , 2:(sample_num+1)])          # 小于50,全画
	}else{
               res <- res[order(res$pvalue), ]
               data <- data.matrix(res[1:50, 2:(sample_num+1)])           #  相对丰度最高的50个
	}
heatmap_file <- paste(output, "/heatmap.top50.pdf", sep="")
pdf(heatmap_file , width=as.numeric(opts$w), height = as.numeric(opts$h) )

# 样本颜色
my_plot_color <- as.character(sampletable[colnames(data),])
my_plot_color[my_plot_color == case_group]    <- 'red'
my_plot_color[my_plot_color == control_group] <- 'blue'
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
library(gplots)
data <- data.matrix(res[res$type != "Not DEG", 2:(sample_num+1)])
heatmap_file <- paste(output, "/heatmap.pdf", sep="")
pdf(heatmap_file , width=as.numeric(opts$w), height = as.numeric(opts$h) )

# 样本颜色
my_plot_color <- as.character(sampletable[colnames(data),])
my_plot_color[my_plot_color == case_group]    <- 'red'
my_plot_color[my_plot_color == control_group] <- 'blue'
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