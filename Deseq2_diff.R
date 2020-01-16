#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)
"Usage: Deseq2_diff.R -i <file> -o <dir> --case_group_name <string> --control_group_name <string> --case_sample <string> --control_sample <string> [--Rlib <dir> --heatmap_width <int>  --heatmap_height <int>]
Options:
  -i, --input <file>              表达量矩阵，每一行是基因，每一列是样本，第一列是基因名
  -o, --output_dir <dir>          结果输出目录
  --case_group_name <string>      case组名称 
  --control_group_name <string>   control组名称 
  --case_sample <string>          case组样本列表，用“逗号分隔”
  --control_sample <string>       control组样本列表，用“逗号分隔” 
  --Rlib <dir>                    R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]
  --heatmap_width <int>           热图宽度 [default: 12]
  --heatmap_height <int>          热图高度 [default: 12]" -> doc

opts           <- docopt(doc)
input          <- opts$input
case_group     <- opts$case_group_name
control_group  <- opts$control_group_name
case           <- opts$case_sample
control        <- opts$control_sample
output         <- opts$output_dir
Rlib           <- opts$Rlib
heatmap_width  <- opts$heatmap_width
heatmap_height <- opts$heatmap_height
 
.libPaths(Rlib)
library(DESeq2)

x <- read.table(input, header=T, sep="\t", row.names=1, comment.char="", check.names=F)

control_count           <- sapply(unlist(strsplit(control, split=",")), function(t){x[[t]] })
colnames(control_count) <- unlist(strsplit(control, split=","))

case_count              <- sapply(unlist(strsplit(case, split=",")), function(t){ x[[t]] })
colnames(case_count)    <- unlist(strsplit(case, split=","))

count                   <- cbind(data.frame(control_count, check.names=F), data.frame(case_count, check.names=F))
rownames(count)         <- rownames(x)


countData <- count[apply(count, 1, sum) > 0 , ]



colData   <- data.frame(row.names = colnames(countData), 
	                    condition = rep(
	                    	         	c(control_group, case_group),
	                    	         	times  = c(ncol(control_count), ncol(case_count)) ,
	                    	         	levels = c(control_group, case_group)
	                    	         ) 	        
	                    )
colData$condition <- relevel(colData$condition, ref = control_group)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

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
pdf(heatmap_file , width=as.numeric(heatmap_width), height = as.numeric(heatmap_height) )

# 样本颜色
my_plot_color <- as.character(colData[colnames(data),])
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
pdf(heatmap_file , width=as.numeric(heatmap_width), height = as.numeric(heatmap_height) )

# 样本颜色
my_plot_color <- as.character(colData[colnames(data),])
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
