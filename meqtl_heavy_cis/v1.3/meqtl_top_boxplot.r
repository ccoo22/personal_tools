#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: meqtl_top_boxplot.r --output <file> --top_info <file> --exp_file <file> --snv_file <file> [ --y_lab <string> --rlib <dir> ]
Options:
    --output <file>              pdf输出文件，例如 result.pdf
    --top_info <file>            要绘制boxplot的cis分析结果，重点要有表头： snps gene pvalue 三列信息
    --exp_file <file>            表达量数据矩阵，第一列是id，后面n列都是样本，数值型数据，如果缺失，空着。含有表头。务必保证样本名和顺序与snp_file中的一致。
    --snv_file <file>            突变量化矩阵，第一列是id, 后面n列都是样本，数值型数据，如果缺失，空着。含有表头。通常是 基因型
    --y_lab <string>             更改y轴标签 [default: EXP]
    --rlib <dir>                 R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts            <- docopt(doc, version = 'Program : cis boxplot  \n        甘斌 129\n')
exp_file        <- opts$exp_file
snv_file        <- opts$snv_file
output          <- opts$output
top_info        <- opts$top_info
y_lab           <- opts$y_lab
rlib            <- opts$rlib

if(!is.null(rlib)) .libPaths(rlib)

# 测试参数
# exp_file        <- "exp.raw.txt"
# exp_loc_file    <- "exp.pos.txt"
# snv_file        <- "snp.code.txt"
# snv_loc_file    <- "snp.pos.txt"
# output          <- "b.txt"
# cov_file        <- NULL
# cis_pval        <- 0.05
# cisDist_val     <- 1e6
 
plot_type = 'dotplot'
#########
library(ggpubr)

message("load data")
gene_exp    = read.table(exp_file, header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t", check.names = F)
snp_val     = read.table(snv_file, header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t", check.names = F)
top_val     = read.table(top_info, header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = F)

message("start plot")
pdf(output, width = 14, height = 10)
plot_list = list()  # 图像保存在列表里
pic_count = 0
for(row in 1:nrow(top_val))
{
    pic_count = pic_count + 1
    snp_name  = as.character(top_val[row, 'snps'])
    gene_name = as.character(top_val[row, 'gene'])
    pvalue    = format(top_val[row, 'pvalue'], scientific=TRUE, digit=4) # 科学计数，更美观，字符型数据

    data_plot = data.frame(SNP = unlist(snp_val[snp_name, ]), EXP = unlist(gene_exp[gene_name, ]) )
    data_plot = data_plot[complete.cases(data_plot), , drop = FALSE] # 去掉NA
    data_plot = data_plot[apply(data_plot == "", 1, sum) == 0, ]  # 去掉空值
    if(is.null(nrow(data_plot))) next #只有一个样本或一个位点，无法进行分析

    # 柱状图
    if(plot_type == 'barplot')
    {
        p <- ggbarplot(data_plot, x="SNP", y="EXP", fill = "SNP",
          palette = "npg", #杂志nature的配色
          caption = paste0('pvalue = ', pvalue),
          # palette = "aaas", #杂志Science的配色
          # palette = "jco", #按jco杂志配色方案
          add = "mean_se", # 取均值，并添加标准差曲线
          error.plot = "upper_errorbar", # 只显示上曲线
          title = paste0(snp_name, "\n", gene_name)
          ) + ylab(y_lab) + guides(fill = FALSE) + theme(plot.title = element_text(hjust = 0.5)) # 去掉legend,标题居中
    }

    # 点状图
    if(plot_type == 'dotplot')
    {
        bin_width = (max(data_plot$EXP) - min(data_plot$EXP)) / 40 # 设定点的大小
        if(bin_width == 0) bin_width = 0.00001

        p <- ggdotplot(data_plot, x="SNP", y="EXP", fill = "SNP",
          binwidth = bin_width,
          caption = paste0('pvalue = ', pvalue),
          palette = "npg", #杂志nature的配色
          # palette = "aaas", #杂志Science的配色
          # palette = "jco", #按jco杂志配色方案
          add = c("boxplot"),  # 增加boxplot
          add.params = list(color = 'SNP'),
          title = paste0(snp_name, "\n", gene_name)
          )  + ylab(y_lab) + guides(fill = FALSE, color = FALSE) + theme(plot.title = element_text(hjust = 0.5))  # 去掉legend,标题居中

    }


    # 箱线图状图（与点状图的优化结果相似）
    if(plot_type == 'boxplot')
    {
        p <- ggboxplot(data_plot, x="SNP", y="EXP", color = "SNP",
          palette = "npg", #杂志nature的配色
          caption = paste0('pvalue = ', pvalue),
          # palette = "aaas", #杂志Science的配色
          # palette = "jco", #按jco杂志配色方案
          add = 'jitter', # 取均值，并添加标准差曲线
          add.params = list(color = 'SNP'),
          outlier.size = -1, # 隐藏离群点，否则离群点会在图上出现两次（boxplot绘制离群点，jitter再绘制所有点，从而使离群点绘制两次）
          title = paste0(snp_name, "\n", gene_name)
          ) + ylab(y_lab) + guides(fill = FALSE, color = FALSE) + theme(plot.title = element_text(hjust = 0.5))  # 去掉legend,标题居中
    }
     plot_list[[pic_count]] <- ggpar(p, font.title = c(10, "bold.italic"))
}

if(length(plot_list) > 0){
    cat('start output pic\n')
    #  ggarrange(plotlist = plot_list, nrow =arrange_pic[1], ncol = arrange_pic[2])  # 设定画布
    tmp <- capture.output( ggarrange(plotlist = plot_list, nrow = 3, ncol = 3) )
    # ggarrange 输出图片，刷屏，故抑制掉
}
dev.off()

