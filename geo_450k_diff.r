#!/home/genesky/software/r/3.5.1/bin/Rscript

.libPaths("/home/genesky/software/r/3.5.1/lib64/R/library/")
library(docopt)
"Usage: geo_450k_diff.r --series_number <string> --case_sample <string> --control_sample <string> --output_prefix <file>  [--pvalue_cutoff <numeric> --case_group_name <string>  --control_group_name <string> ] 
Options:
   --series_number <string>        GEO 编号列表，包含了case_sample/control_sample 样本的芯片数据。多个GEO编号用逗号分隔。支持不同GEO中的样本进行对比，前提是芯片一致。例如： 'GSE82107'
   --case_sample <string>          case组样本编号，样本之间用逗号分隔 例如： 'GSM2183539,GSM2183540,GSM2183541,GSM2183542,GSM2183543,GSM2183544,GSM2183545,GSM2183546,GSM2183547,GSM2183548'
   --control_sample <string>       control组样本编号，样本之间用逗号分隔 例如： 'GSM2183532,GSM2183533,GSM2183534,GSM2183535,GSM2183536,GSM2183537,GSM2183538'
   --case_group_name <string>      case组名 [default: case]
   --control_group_name <string>   control组名 [default: control]
   --output_prefix <file>          输出文件的前缀，例如： ./result  。 最终生成的结果为 result.diff.xls。
   --pvalue_cutoff <numeric>       pvalue 阈值，仅输出小于该阈值的结果 [default: 1]" -> doc

opts                   <- docopt(doc, version = 'Program : geo_mrna_chip_diff v1.0 \n          甘斌 129\n')
series_number          <- opts$series_number
case_sample_list       <- opts$case_sample
control_sample_list    <- opts$control_sample
output_prefix          <- opts$output_prefix
pvalue_cutoff          <- as.numeric(opts$pvalue_cutoff)
case_group_name        <- opts$case_group_name
control_group_name     <- opts$control_group_name

case_samples    = unlist(strsplit(case_sample_list, ','))
control_samples = unlist(strsplit(control_sample_list, ','))
sample_group  <- c(rep(case_group_name, length(case_samples)), rep(control_group_name, length(control_samples)))


# https://www.ncbi.nlm.nih.gov/geo/info/geo2r.html
##################
# （1） 加载R包
##################
library(Biobase)
library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)

##################
# （2） 加载每个series的数据
##################
data_matrix <- ''  # 汇总的表达量矩阵，matrix
chip_version <- ''  # 芯片版本
probe_annotation <- ''  # 探针注释
series_count <- 0
for (series in unlist(strsplit(series_number, ',')) )
{   
    series_count <- series_count + 1
    message("loading series : ", series)
    gset_tmp <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE)
    gset_tmp <- gset_tmp[[1]]

    # 提取芯片版本、表达矩阵、注释数据库
    chip_version_tmp <- gset_tmp@annotation
    data_matrix_tmp <- exprs(gset_tmp)

    # 芯片版本核对
    if(chip_version != '' && chip_version != chip_version_tmp)
    {
        message("输入的series的芯片版本不一致，不能放在一起分析")
        q()
    }else {
       chip_version <- chip_version_tmp
    }

    # 保存芯片探针注释信息
    if(series_count == 1)
    {
        probe_annotation <- gset_tmp@featureData@data
    }
    
    # 保存表达矩阵
    data_matrix_tmp <- exprs(gset_tmp)
    # log2 转换
    qx <- as.numeric(quantile(data_matrix_tmp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0) ||
          (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) 
    { 
        data_matrix_tmp[which(data_matrix_tmp <= 0)] <- NaN
        data_matrix_tmp <- log2(data_matrix_tmp) 
    }
    # 甲基化450K数据处理，log2是不会执行的

    if(series_count == 1)
    {
        data_matrix <- data_matrix_tmp
    }else {
        # 两个批次数据行数不对
        if(nrow(data_matrix) != nrow(data_matrix_tmp)) 
        { 
            message("输入的series的探针数量不一致，不能放在一起分析"); 
            q(); 
        }
        # 两个批次数据的探针编号不对应
        if(identical(sort(rownames(data_matrix)), sort(rownames(data_matrix_tmp))) == FALSE) 
        { 
            message("输入的series的探针名称不一致，不能放在一起分析"); 
            q(); 
        }

        # 可以合并
        match_pos <- match(rownames(data_matrix), rownames(data_matrix_tmp))  # 根据探针名确定行对应关系
        data_matrix <- cbind(data_matrix, data_matrix_tmp[match_pos,])
    }
}


##################
# （3） 检查样本编号是否有问题 并 制作新的BioBase数据库
##################
message("build BioBase dataset")
if( sum( c(case_samples, control_samples) %in% colnames(data_matrix) ) != (length(case_samples) + length(control_samples)) ) { message("部分输入的case/control样本没有在series中找到，请仔细核对是否写错"); q(); }

data_matrix <- data_matrix[, c(case_samples, control_samples)]  # 仅保留需要的样本
gset <- ExpressionSet(assayData = data_matrix,
                       annotation = chip_version
                       ) 
##################
# (4) 差异分析
##################
message("diff analysis")
sml <- paste("G", c(rep(1, length(case_samples)), rep(0, length(control_samples))), sep="")

fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", n=Inf, p.value=pvalue_cutoff)  # 最终差异分析结果

##################
# (5) 差异分析结果整理输出
##################
message("write diff result")
caseExpr = apply(data_matrix[rownames(tT), case_samples], 1, mean)
controlExpr = apply(data_matrix[rownames(tT), control_samples], 1, mean)
result <- data.frame(probe = rownames(tT),
                 data_matrix[rownames(tT), ],
                 caseExpr = caseExpr,
                 controlExpr = controlExpr,
                 averageExpr = tT$AveExpr,
                 methy_diff = tT$logFC,
                 log2FC = log2(caseExpr/controlExpr),
                 pvalue = tT$P.Value,
                 fdr = tT$adj.P.Val,
                 CHR = probe_annotation[rownames(tT), 'CHR'],
                 MAPINFO = probe_annotation[rownames(tT), 'MAPINFO'],
                 UCSC_RefGene_Name = probe_annotation[rownames(tT), 'UCSC_RefGene_Name'],
                 UCSC_RefGene_Group = probe_annotation[rownames(tT), 'UCSC_RefGene_Group'],
                 Relation_to_UCSC_CpG_Island = probe_annotation[rownames(tT), 'Relation_to_UCSC_CpG_Island']
                 )
result$type <- "Not DEG"
result$type[result$pvalue < 0.05 & result$methy_diff >= 0.1 ] <- "Up"
result$type[result$pvalue < 0.05 & result$log2FC <= -0.1] <- "Down"
result$type <- factor(result$type, levels = c("Up", "Down", "Not DEG"))

file_diff = paste0(output_prefix, ".diff.xls")
message("output diff file : ", file_diff)
write.table(result, file=file_diff, row.names=F, sep="\t", quote = FALSE)

# 差异基因列表
file_diff_gene = paste0(output_prefix, ".diff.gene.xls")
diff_gene <- sort(unique(unlist(strsplit(as.character(result[result$type != 'Not DEG', 'UCSC_RefGene_Name']), ";"))))
diff_gene = data.frame(gene=diff_gene)
write.table(diff_gene, file=file_diff_gene, row.names=F, col.names=F, sep="\t", quote = FALSE)

##################
# （6）PCA  暂时不绘制，数据做过log2处理，存在极限值
##################
# mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
# mycol <- colors()[rep(mycol, 50)]
# myshape <- rep(c(15,16,17,18,19,20,21,22,23,24,25,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14),4)

# pca_coordate <- paste(output_prefix, ".pca.site.txt", sep="")
# pca_pdf <- paste(output_prefix, ".pca.pdf", sep="")
# message("plot pca : ", pca_pdf)

# pca <- prcomp(t(data_matrix))
# summaryInfo <- summary(pca)
# pc1Pro <- 100 * summaryInfo$importance['Proportion of Variance', 'PC1'] # PC1主元解释百分比
# pc2Pro <- 100 * summaryInfo$importance['Proportion of Variance', 'PC2'] # PC2主元解释百分比
# scores <- as.data.frame(pca$x)
# scores$Group = factor(sample_group, levels = unique(sample_group)) # 设定成factor数据，且定义levels与输入顺序一致，这样就可以固定pca图上的分组顺序了
# write.csv(scores, file=pca_coordate)

# pdf(file=pca_pdf, width=13, height=13)
# p = ggplot(data=scores, aes(x = PC1, y = PC2, colour = Group, shape = Group)) +
#         geom_hline(yintercept = 0, colour = "gray65") +
#         geom_vline(xintercept = 0, colour = "gray65") +
#         scale_shape_manual(values = myshape ) +
#         scale_color_manual(values = mycol ) +
#         geom_point(size = 3, alpha = 1) +
#         ggtitle('PCA') + xlab(paste("PC1 ", "(", pc1Pro, "%)", plot_x_lab, sep="")) + ylab(paste("PC2 ", "(", pc2Pro, "%)", sep="")) + theme(plot.title = element_text(hjust = 0.5))

# p
# dev.off()


##################
# （7）MA plot
##################
ma_file <- paste(output_prefix, ".MA.png", sep="")
message("plot MA : ", ma_file)
png(ma_file)
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
volcano_file <- paste(output_prefix, ".volcano.png", sep="")
message("plot volcano : ", volcano_file)
png(volcano_file)
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
if(nrow(result_diff) < 2)
{
    message("差异结果数量小于2，不绘制热图")
    q()
}

# Top50 heatmap plot
if( nrow(result_diff) < 50 ){
               data <- data.matrix(result_diff[ , c(case_samples, control_samples)])          # 小于50,全画
	}else{
               result_diff <- result_diff[order(result_diff$pvalue), ]
               data <- data.matrix(result_diff[1:50, c(case_samples, control_samples)])           #  相对丰度最高的50个
	}
heatmap_file <- paste(output_prefix, ".heatmap.top50.png", sep="")
message("plot heatmap top50 : ", heatmap_file)

png(heatmap_file  )
# 样本注释
annotation_group <- data.frame(Group = factor(sml, levels = unique(sml)))
rownames(annotation_group) <- c(case_samples, control_samples)
 
myheatcol = colorRampPalette(c('green','black','red'))(100)
pheatmap(data,
             scale = 'row',
             cluster_rows = T,
             cluster_cols = T,
             color = myheatcol,
             show_rownames = T,
             show_colnames = T,
             annotation_col = annotation_group,
    )
dev.off()



# heatmap plot
data <- data.matrix(result_diff[, c(case_samples, control_samples)])
heatmap_file <- paste(output_prefix, ".heatmap.png", sep="")
message("plot heatmap all : ", heatmap_file)

png(heatmap_file )

pheatmap(data,
             scale = 'row',
             cluster_rows = T,
             cluster_cols = T,
             color = myheatcol,
             show_rownames = F,
             show_colnames = T,
             annotation_col = annotation_group,
    )
dev.off()
