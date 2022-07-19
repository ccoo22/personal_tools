#!/home/genesky/software/r/4.1.2/bin/Rscript
library(docopt)
"Usage: boxplot.r  -i <file> -g <file> --pdf <file> [ --y_position <numeric>  --input_format <string> --add_diff <string>  --union --if_add <string> --outlier_size <int> --palette <string> --color_need <string> --pdf_width <int> --pdf_height <int> --auto_width --x_angle <int> --x_lab <string> --y_lab <string> ]
Options:
   -i, --input <file>                   the input file, each row is a point, each column is sample,行是位点，列是样本，文件必须有表头
   -g, --group <file>                   the group file,包含表头 sample group,绘图分组顺序和group列中分组出现顺序一致
   --input_format <string>              input文件格式 col/row。col表示exp文件每一列对应一个样本；row表示exp文件每一行对应一个样本 [default: col]
   --pdf <file>                         the pdf output file 
   --union                              gene结果汇总绘图
   --if_add <string>                    add character vector for adding another plot element. [default: jitter]
                                        Allowed : none, dotplot, jitter, boxplot, point, mean, mean_se, mean_sd, 
                                        mean_ci, mean_range, median, median_iqr, median_hilow, median_q1q3, median_mad, median_range
   --outlier_size <int>                 离群点大小,-1隐藏离群点 [default: -1]
   --palette <string>                   绘图颜色,允许指定调色板,杂志配色 [default: npg]
   --color_need <string>                指定颜色,需和分组数量一致,示例'#FF0000,#00FF00' [default: none]
   --pdf_width <int>                    设定 pdf 宽度 [default: 7]  
   --pdf_height <int>                   设定 pdf 高度 [default: 7] 
   --auto_width                         自动设置pdf宽度
   --add_diff <string>                  做差异分析, t.test / wilcoxon
   --y_position <numeric>               如果有 add_diff 参数， 需要指定 y_position pvalue显示的位置,该参数只适合于 --union [default: 1000]
   --x_angle <int>                      x轴上标签倾斜角度 [default: 45]
   --x_lab <string>                     设定x轴标签 [default: Gene]
   --y_lab <string>                     设定y轴标签 [default: Value]" -> doc


opts                     <- docopt(doc, version = 'Program : boxplot plot \n')
input                    <- opts$input
group                    <- opts$group
input_format             <- opts$input_format
output_pdf               <- opts$pdf
union                    <- opts$union
if_add                   <- opts$if_add
outlier_size             <- opts$outlier_size
color_used               <- opts$palette
color_need               <- opts$color_need
pdf_width                <- opts$pdf_width
pdf_height               <- opts$pdf_height
auto_width               <- opts$auto_width
x_angle                  <- opts$x_angle
x_lab                    <- opts$x_lab
y_lab                    <- opts$y_lab
add_diff                 <- opts$add_diff
y_position               <-  as.numeric(opts$y_position)
outlier_size = as.numeric(outlier_size)
x_angle     = as.numeric(x_angle)
pdf_width   = as.numeric(pdf_width)
pdf_height  = as.numeric(pdf_height)

# # # 测试用参数
# input <- "/home/lingq/script/tools/plot_ggpubr/test/gene2.txt"
# group <- "/home/lingq/script/tools/plot_ggpubr/test/group.txt"
# if_add = "jitter"
# union = "TRUE"

options(warn=-1)
library(ggpubr, quietly = TRUE)
suppressMessages(library(tidyr, quietly = TRUE))

data <- read.table(input, sep='\t', header = TRUE, row.names = 1, check.names=FALSE, stringsAsFactors = F) # 读取数据，行为样本
group <- read.table(group, sep='\t', header = TRUE, check.names=FALSE, stringsAsFactors = F, colClasses = 'character')

if(input_format == 'row')
{
    data = t(data)
}

# 数据检测
colnames(group) = c("sample","group")
sample_list <- as.character(group$sample)
if( "TRUE" %in% duplicated(sample_list) ){
    cat('[Error] Duplicate sample in group :', sample_list[duplicated(sample_list)], '\n')
    q()
}
if(sum(sample_list %in% colnames(data)) != length(sample_list)){
    cat('[Error] The sample to be analyzed in the group file does not exist in the data :', sample_list[(!sample_list %in% colnames(data))], '\n')
    q()
}
if(color_need != "none"){
    color_used <- as.character(unlist(strsplit(color_need, split=",")))
    if(length(color_used) != length(unique(group$group))){
      cat('[Error] 所给颜色数目和分组数目不一致 :', color_need, '\n')
      q()
    }
}

# 获取指定分析的样本
rownames(group) = group$sample
data_need <- data[, sample_list]
if( "character" %in% sapply(data_need, class) ){
    cat('[Error] 数据中存在非数值型变量\n')
    q()
}
data_need <-as.matrix(data_need)

################### 循环绘图 ######################################
# 设定绘图规划 （一个组的柱状图大约占用2.5个宽度）
if(auto_width == TRUE)
{
    group_count  = length(unique(group$group)) # 分组数量
    pdf_width = 2.5 * group_count
    if(union==TRUE){
      group_count = length(rownames(data_need))
      pdf_width = group_count * length(unique(group$group))
    }
}

pdf(output_pdf, width = pdf_width, height = pdf_height)
if(union==TRUE){
    data_all <- data.frame(t(data_need),check.names=FALSE)
    groups_unique  = unique(group[rownames(data_all),'group'])
    data_all$Group <- group[rownames(data_all),'group']
    data_all$Group <- factor(data_all$Group,levels=unique(data_all$Group))
    all <- gather(data_all,Gene,Value,-Group)
    all <- all[complete.cases(all), ,drop = FALSE] # 去掉NA

    cat('start output pic\n')
    p <- ggboxplot(all, x="Gene", y="Value", color = "Group",
      palette = color_used,
      #palette = c("#FF0000","#00FF00"),
      add = if_add,
      add.params = list(color = 'Group'),
      outlier.size = outlier_size, # 隐藏离群点，否则离群点会在图上出现两次（boxplot绘制离群点，jitter再绘制所有点，从而使离群点绘制两次）
      ) + labs(x=x_lab, y=y_lab, title = "ALL") + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = x_angle, vjust = 0.5))
    if(!is.null(add_diff)){
        if(add_diff == 'wilcoxon')
        {
            add_diff = NULL
        }
        stat.test <- compare_means(
           Value ~ Group, data = all, group.by = "Gene",
          method = add_diff
        )
        p = p + stat_pvalue_manual(
                stat.test, x = "Gene" , y.position = y_position,
                label = "p.signif",
                position = position_dodge(0.8))
    }
    print(p)
    dev.off()

}else{
    plot_list = list()  # 图像保存在列表里
    pic_count = 0
    for(row_name in rownames(data_need))
    {   
        pic_count = pic_count + 1
        data_plot = data.frame(Group = factor(group$group, levels = unique(group$group)), Value =data_need[row_name,])
        data_plot = data_plot[complete.cases(data_plot), , drop = FALSE] # 去掉NA
        if(is.null(nrow(data_plot))) next #只有一个样本或一个位点，无法进行分析  

        p <- ggboxplot(data_plot, x="Group", y="Value", color = "Group", 
          palette = color_used, #npg杂志nature的配色
          # palette = "aaas", #杂志Science的配色
          # palette = "jco", #按jco杂志配色方案
          add = if_add,
          add.params = list(color = 'Group'),
          #outlier.shape = 19,
          outlier.size = outlier_size, # 隐藏离群点，否则离群点会在图上出现两次（boxplot绘制离群点，jitter再绘制所有点，从而使离群点绘制两次）
          ) + labs(x=x_lab, y=y_lab, title = row_name) + guides(fill = FALSE, color = FALSE) + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = x_angle, vjust = 0.5))  # 去掉legend,标题居中   
        
        if(!is.null(add_diff)){
            if(add_diff == 'wilcoxon')
            {
                add_diff = NULL
            }
            p = p + stat_compare_means(method = add_diff)
        }
        plot_list[[pic_count]] <- p
    }
     

    if(length(plot_list) > 0){
        cat('start output pic\n') 
        #  ggarrange(plotlist = plot_list, nrow =arrange_pic[1], ncol = arrange_pic[2])  # 设定画布
        tmp <- capture.output( ggarrange(plotlist = plot_list, nrow =1, ncol = 1) )  
        # ggarrange 输出图片，刷屏，故抑制掉
    }
    dev.off()
}
