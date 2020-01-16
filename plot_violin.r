#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: plot_violin.r  -i <file> -o <pdf file> [ --arrange_pic <int,int> --x_lab <string> --y_lab <string> --title <string> --rlib <dir> ]

Options:
    -i, --input <file>        输入文件，第一列为样本名，第二列为频数或频率
    -o, --output <pdf file>   输出pdf文件路径,示例：./a.pdf
    --arrange_pic <int,int>   图像展示方式,一页显示n*m个样本 [default: 1,5]
    --x_lab <string>          x轴名称 [default: sample]
    --y_lab <string>          y轴名称 [default: freq]
    --title <string>          标题    [default: Sample SNV frequency violin plot]
    --rlib <dir>              R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，小提琴图\n')
input             <- opts$input
output            <- opts$output
arrange_pic       <- opts$arrange_pic
x_lab             <- opts$x_lab
y_lab             <- opts$y_lab
title             <- opts$title

arrange_pic <- as.integer(unlist(strsplit(arrange_pic, ',')))

# 加载R包
message('加载ggpubr')
library(ggpubr, quietly = TRUE)

# 读入数据
message('读入数据')
data <- read.table(input, sep='\t', header = TRUE, check.names=FALSE) # 读取数据，行为样本
colnames(data) <- c('sample', 'freq')

# 包含的样本
sample_name_list <- sort(as.character(unique(data$sample)))
sample_count     <- length(sample_name_list)

# 颜色模版
mycol <- c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)# 颜色模版
mycol <- colors()[rep(mycol, 50)] 

# 开始绘图
message("开始绘图")
pdf_width  = 2 * arrange_pic[2] # 一个样本2个宽度单位
pdf_height = 7 * arrange_pic[1] # 一行7个单位高度
pdf(output, width = pdf_width, height = pdf_height)

plot_list = list()  # 图像保存在列表里
for(pic_count in 1:ceiling(sample_count / arrange_pic[2]))
{	
    # 当前行要绘制的样本
    start_index = (pic_count - 1) * arrange_pic[2] + 1
    end_index   = pic_count * arrange_pic[2]
    if(end_index > sample_count) end_index = sample_count

    plot_sample = sample_name_list[start_index:end_index]
    show_sample = plot_sample # 绘图显示样本名称，主要用于处理empty样本
    color_sample = mycol[start_index:end_index]

    message("绘制样本： ", paste0(plot_sample, collaple=','))
    # 绘图数据提取
	data_plot = data[data$sample %in% plot_sample, ]

    # 最后一幅图可能需要填充空值，从而保证图形一致性
    if(length(plot_sample) != arrange_pic[2])
    {   
        empty_sample = paste0('unknown', 1:(arrange_pic[2]-length(plot_sample)))
        empty_data = data.frame(sample = empty_sample, freq = 0)
        data_plot = rbind(data_plot, empty_data)
       
        plot_sample <- c(plot_sample, empty_sample)
        show_sample = c(show_sample, rep('', length(empty_sample)))
 
    }
   
    p <- ggplot(data = data_plot, aes(x = sample, y = freq)) +
        geom_violin(aes(fill = sample) ) +
        scale_fill_manual(breaks = plot_sample, values = color_sample) +
        guides(fill = FALSE) + 
        geom_hline(aes(yintercept = 0.75), colour="#990000", linetype="dashed") + 
        geom_hline(aes(yintercept = 0.5), colour="#990000", linetype="dashed") + 
        geom_hline(aes(yintercept = 0.25), colour="#990000", linetype="dashed") + 
        theme(text = element_text(size = 20)) + labs(x = x_lab, y = y_lab) + 
        ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) + scale_x_discrete(breaks = plot_sample, labels = show_sample)

	# p <- ggviolin(data_plot, x="sample", y="freq",
	# 	color = "black", fill = 'sample', palette = color_sample
	# 	) + guides(fill = FALSE, color = FALSE)

    plot_list[[pic_count]] <- p

}
if(length(plot_list) > 0){
    cat('start output pic\n') 
    #  ggarrange(plotlist = plot_list, nrow =arrange_pic[1], ncol = arrange_pic[2])  # 设定画布
    tmp <- capture.output( ggarrange(plotlist = plot_list, nrow =arrange_pic[1], ncol = 1) )  
    # ggarrange 输出图片，刷屏，故抑制掉
}
dev.off()







    # $R->send(qq` mycol <- c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)\n`); # 颜色模版
    # $R->send(qq` mycol <- colors()[rep(mycol, 50)] \n`);

    # $R->send(qq` plot_data <- read.table("$freq_vioplot", head = T, colClasses = "character") \n`); # 样本名存在0开头的纯数字格式，R会自动删除，所以按照字符读入，后面再转成数字
    # $R->send(qq` plot_data\$Frequency <- as.numeric(plot_data\$Frequency) \n`);
    # $R->send(qq` plot_data\$Sample <- factor(plot_data\$Sample) \n`);
    # $R->send(qq` sample_name <- unique(plot_data\$Sample) \n`);
    # $R->send(qq` sample_color <- mycol[1 : length(sample_name)] \n`);

    # $R->send(qq` pdf("$pdf_vioplot", width = $pdf_width, height = 10)\n`);
    # $R->send(qq` ggplot(data = plot_data, aes(x = Sample, y = Frequency)) +
    #     geom_violin(aes(fill = Sample) ) +
    #     scale_fill_manual(breaks = sample_name, values = sample_color) +
    #     guides(fill = FALSE) + 
    #     geom_hline(aes(yintercept = 0.75), colour="#990000", linetype="dashed") + 
    #     geom_hline(aes(yintercept = 0.5), colour="#990000", linetype="dashed") + 
    #     geom_hline(aes(yintercept = 0.25), colour="#990000", linetype="dashed") + 
    #     theme(text = element_text(size = 20)) +
    #     ggtitle("Sample SNV frequency violin plot") + theme(plot.title = element_text(hjust = 0.5)) \n`);
    # $R->send(qq`dev.off()\n`); 
