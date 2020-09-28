#!/home/genesky/software/r/3.5.1/bin/Rscript
 
library(docopt)
"Usage: barplot.r  -i <file> --pdf <file> --sample_group_file <file> [ --plot_type <string> --use_title_column <string>  --charactor_combine_title <string> --arrange_pic <string> --pdf_width <int> --pdf_height <int> --auto_width --pvalue_mark_file <file> --pvalue_mark_title_column <string> --pvalue_mark_anno <string> --y_lab <string> ]
Options:
   -i, --input <file>                   输入文件，行是基因，列是样本，文件必须有表头。
   --pdf <file>                         the pdf output file 
   --sample_group_file <file>           样本分组文件，第一列样本名，第二列样本分组，包含表头。仅对该文件中包含的样本绘图
   --plot_type <string>                 绘图类型，支持: barplot、dotplot、boxplot [default: barplot]
   --use_title_column <string>          使用哪几列作为行名，输入列编号 example:1,2 [default: 1]
   --charactor_combine_title <string>   上面输入多列时，用什么分隔符分隔 [default: |] 
   --arrange_pic <string>               pdf中一页中显示的方式, example: 2,3.表示一页显示6张图，2行3列 [default: 1,1]  
   --pdf_width <int>                    设定 pdf 宽度 [default: 7]  
   --pdf_height <int>                   设定 pdf 高度 [default: 7] 
   --auto_width                         自动设置pdf宽度
   --pvalue_mark_file <file>            用于注释P值的输入文件，文件必须有表头 注意：pvalue_mark的三个参数同时存在时，才能进行注释
   --pvalue_mark_title_column <string>  注释文件中,用于识别绘图名称的注释列，组合结果要与--use_title_column的结果保持一致，以此区分不同图片，否则无法在对应的图中添加注释
   --pvalue_mark_anno <string>          注释方案，3列构成的固定格式：组1名称,组2名称,注释列编号.注意：组名称要与group_name_list一致，否则也不会加注释
   --y_lab <string>                     设定y轴标签 [default: Methylation]" -> doc

opts                     <- docopt(doc, version = 'Program : barplot plot  v2.0 \n          甘斌 129\n')
input                    <- opts$input
output_pdf               <- opts$pdf
sample_group_file        <- opts$sample_group_file
plot_type                <- opts$plot_type
use_title_column         <- as.numeric(unlist(strsplit(opts$use_title_column, split=",")))
charactor_combine_title  <- opts$charactor_combine_title
arrange_pic              <- as.numeric(unlist(strsplit(opts$arrange_pic, split=",")))
pdf_width                <- as.numeric(opts$pdf_width)
pdf_height               <- as.numeric(opts$pdf_height)
auto_width               <- opts$auto_width
y_lab                    <- opts$y_lab

# P值注释参数
pvalue_mark_file         <- opts$pvalue_mark_file
pvalue_mark_title_column <- opts$pvalue_mark_title_column
pvalue_mark_anno         <- opts$pvalue_mark_anno

# # 检测用
# cat(input, "\n")
# cat(output_pdf, "\n")
# cat(sample_name_list, "\n")
# cat(group_name_list, "\n")
# cat(plot_type, "\n")
# cat(use_title_column, "\n")
# cat(charactor_combine_title, "\n")
# cat(target, "\n")
# cat(snp_stat, "\n")
# cat(arrange_pic, "\n")
# cat(auto_width, "\n")
# cat(pdf_width, "\n")
# cat(pdf_height, "\n")
# cat(pvalue_mark_file, "\n")
# cat(pvalue_mark_title_column, "\n")
# cat(pvalue_mark_anno, class(pvalue_mark_anno), "\n")
 
# q()

# # # 测试用参数
# input                   <- '/home/ganb/work/research/methyltarget_diff_program/19B0314A/report/methylation_result/2.methyl.site.txt'
# output_pdf              <- './Site_Barplot.Case.VS.Control.pdf'
# sample_name_list        <- 'SZ01,SZ02,SZ04,57_0h,L35_0h,T35-1,016_0h,30_0h,41_0h,47_0h,61_0h,46_0h,'
# group_name_list         <- 'Case,Case,Case,Case,Case,Case,Control,Control,Control,Control,Control,Control,'
# use_title_column        <- '1,2'
# charactor_combine_title <- '|'
# target                  <- 'NA'
# snp_stat                <- 'keep'
# arrange_pic             <- '1,6'
# auto_width              <- TRUE
# pdf_width               <- 7
# pdf_height              <- 7
# pvalue_mark_file          <- '/home/ganb/work/research/methyltarget_diff_program/19B0314A/report/diff_analysis/Site.diff_normal.Case.VS.Control.txt'
# pvalue_mark_title_column  <- '1,2'
# pvalue_mark_anno          <- 'Case,Control,4'
options(warn=-1)
library(ggpubr, quietly = TRUE)
suppressMessages(library(tidyr, quietly = TRUE))
 
data_input <- read.table(input, sep='\t', header = TRUE, check.names=FALSE, stringsAsFactors=F) # 读取数据，行为样本

# 获取样本、分组
data_group <- read.table(sample_group_file, sep='\t', header = TRUE, row.names = 1, check.names=FALSE, colClasses = 'character')  # 读入分组
sample_names = rownames(data_group)
sample_groups = data_group[,1]
 
# 获取指定分析的样本
data_need    = data_input[, sample_names, drop = FALSE]
data_need    = as.matrix(data_need)

# 对行进行命名
if(length(use_title_column) > 1)
{
    rownames(data_need) <- unite(data[, use_title_column], combined, sep = charactor_combine_title)$combined # 标题,有多列的话，使用unite功能进行合并
} else {
    rownames(data_need) <- data_input[, use_title_column]
}


######################## P值注释准备 ######################
# 三个参数同时存在时，才能执行
data_pvalue_mark   <- data.frame() # 注释用数据
group1_name        <- NULL # 分组1名称
group2_name        <- NULL # 分组2名称
pvalue_mark_column <- NULL # 注释信息
if(!is.null(pvalue_mark_file) && !is.null(pvalue_mark_title_column) && !is.null(pvalue_mark_anno))
{
    data_pvalue_mark <- read.table(pvalue_mark_file, sep='\t', header = TRUE, check.names=FALSE, stringsAsFactors=F) # 读取数据，一行对应一个结果
    pvalue_mark_anno <- unlist(strsplit(pvalue_mark_anno, split=",")) # 注释信息
    group1_name        <- pvalue_mark_anno[1] # 组1 名称
    group2_name        <- pvalue_mark_anno[2] # 组2 名称
    pvalue_mark_column <-  as.numeric(pvalue_mark_anno[3]) # 注释列

    pvalue_mark_title_column    <- as.numeric(unlist(strsplit(pvalue_mark_title_column, split=","))) # 注释标题的列
    if(length(pvalue_mark_title_column) > 1)
    { 
        rownames(data_pvalue_mark) <- unite(data_pvalue_mark[, pvalue_mark_title_column], combined, sep = charactor_combine_title)$combined # 标题,有多列的话，使用unite功能进行合并
    } else {
        rownames(data_pvalue_mark) <- data_pvalue_mark[, pvalue_mark_title_column]
    }
}

################### 循环绘图 ######################################
# 设定绘图规划 （一个组的柱状图大约占用2.5个宽度）
if(auto_width == TRUE)
{
    pic_count_each_row = arrange_pic[2] # 一行显示的图片数量
    pic_count_each_col = arrange_pic[1] # 一列显示的图片数量
    group_count        = length(unique(sample_groups)) # 分组数量
    pdf_width  = 2.5 * group_count * pic_count_each_row
    pdf_height = 4 * pic_count_each_col + 1
}

pdf(output_pdf, width = pdf_width, height = pdf_height)
plot_list = list()  # 图像保存在列表里
pic_count = 0
for(row_name in rownames(data_need))
{   
    pic_count = pic_count + 1
    data_plot = data.frame(Group = factor(sample_groups, levels = unique(sample_groups)), Methylation = data_need[row_name, ] )
    data_plot = data_plot[complete.cases(data_plot), , drop = FALSE] # 去掉NA
    if(is.null(nrow(data_plot))) next #只有一个样本或一个位点，无法进行分析
 
    # 柱状图
    if(plot_type == 'barplot')
    {
        p <- ggbarplot(data_plot, x="Group", y="Methylation", fill = "Group", 
        #   palette = "npg", #杂志nature的配色
          palette = "aaas", #杂志Science的配色
          # palette = "jco", #按jco杂志配色方案
          add = "mean_se", # 取均值，并添加标准差曲线
          error.plot = "upper_errorbar", # 只显示上曲线
          title = row_name
          ) + ylab(y_lab) + guides(fill = FALSE) + theme(plot.title = element_text(hjust = 0.5)) # 去掉legend,标题居中      
    }

    # 点状图
    if(plot_type == 'dotplot')
    {   
        bin_width = (max(data_plot$Methylation) - min(data_plot$Methylation)) / 40 # 设定点的大小
        if(bin_width == 0) bin_width = 0.00001
        
        p <- ggdotplot(data_plot, x="Group", y="Methylation", fill = "Group", 
          binwidth = bin_width,
        #   palette = "npg", #杂志nature的配色
          palette = "aaas", #杂志Science的配色
          # palette = "jco", #按jco杂志配色方案
          add = c("boxplot"),  # 增加boxplot
          add.params = list(color = 'Group'),
          title = row_name
          )  + ylab(y_lab) + guides(fill = FALSE, color = FALSE) + theme(plot.title = element_text(hjust = 0.5))  # 去掉legend,标题居中   

    }

    # 箱线图状图（与点状图的优化结果相似）
    if(plot_type == 'boxplot')
    {   
        p <- ggboxplot(data_plot, x="Group", y="Methylation", color = "Group", 
        #   palette = "npg", #杂志nature的配色
          palette = "aaas", #杂志Science的配色
          # palette = "jco", #按jco杂志配色方案
          add = 'jitter', # 取均值，并添加标准差曲线
          add.params = list(color = 'Group'),
          outlier.size = -1, # 隐藏离群点，否则离群点会在图上出现两次（boxplot绘制离群点，jitter再绘制所有点，从而使离群点绘制两次）
          title = row_name
          ) + ylab(y_lab) + guides(fill = FALSE, color = FALSE) + theme(plot.title = element_text(hjust = 0.5))  # 去掉legend,标题居中   
    }
 
    # 添加P值注释
    if(!is.null(group1_name)) # 赋值了，说明可以3个数据都输入了
    {
        if(group1_name %in% data_plot$Group && group2_name %in% data_plot$Group && row_name %in% rownames(data_pvalue_mark)) # 两个分组名称都在，且图片编号能找到
        {   
        	# 位置等定义
            max_value  = max(data_plot$Methylation)
            if(plot_type == 'barplot') # 如果是barplot的话，图形中最大值是mean_se值
            {
            	max_value = max(unlist(lapply(unique(data_plot$Group), function(group_name){
            		    x = data_plot$Methylation[data_plot$Group ==group_name]
            	      mean_se = mean(x) + sd(x) / sqrt(length(x))  
            	      mean_se
            	    } )))
            }
            tip_length  = max_value * 0.02      # 注释位置向下延长的线条长度
            y_position  = max_value * (1 + 0.04) # P值注释位置绘制位置
            mark       =  data_pvalue_mark[row_name, pvalue_mark_column] 
            # pvalue_sci = format(data_pvalue_mark[row_name, pvalue_mark_column], scientific=TRUE, digit=4) # 科学计数，更美观，字符型数据
            # mark       = '' # 显著性标记符号
            # if(grepl(pattern = "\\d", x = pvalue) == TRUE && pvalue <= 0.05)   mark = '*'
            # if(grepl(pattern = "\\d", x = pvalue) == TRUE && pvalue <= 0.01)   mark = '**'
            # if(grepl(pattern = "\\d", x = pvalue) == TRUE && pvalue <= 0.001)  mark = '***'
            # if(grepl(pattern = "\\d", x = pvalue) == TRUE && pvalue <= 0.0001) mark = '****'

            # 准备数据
            data_tmp = data.frame(.y. = c('Methylation'), 
                   group1 = c(group1_name), 
                   group2 = c(group2_name), 
                   P.value = c(mark), 
                   Mark   = c(mark), 
                   y.position = c(y_position)
                )
            suppressWarnings( p <- p + stat_pvalue_manual(data_tmp, xmin = 'group1', xmax = 'group2', y.position = 'y.position', label = "P.value", tip.length = 0.01 ) ) 
            # + stat_pvalue_manual(data_tmp, xmin = 'group1', xmax = NULL, y.position = 'y.position', label = "Mark") 
            # 此处会给出警告 Warning: Ignoring unknown aesthetics: xmin, xmax, annotations, y_position
            # 故，用suppressWarnings进行抑制，否则屏幕会刷满           
        }
    }

    plot_list[[pic_count]] <- p
}
 

if(length(plot_list) > 0){
    cat('start output pic\n') 
    #  ggarrange(plotlist = plot_list, nrow =arrange_pic[1], ncol = arrange_pic[2])  # 设定画布
    tmp <- capture.output( ggarrange(plotlist = plot_list, nrow =arrange_pic[1], ncol = arrange_pic[2]) )  
    # ggarrange 输出图片，刷屏，故抑制掉
}
dev.off()

