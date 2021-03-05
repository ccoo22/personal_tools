#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)
"Usage: barplot_PICRUSt2.r  --inputfile <file> --groupfile <file> --outputfile <file>  [--markfile <file> --categorize1 <string> --categorize2 <string>  --ylab <string>  --rlib <string> ]
Options:
    --inputfile <file>           前两列数据：第一层分类, 第二层分类, 其余列是每个样本的数值，含表头。按照前两列顺序绘图。示例: /home/dongxj/work/develop/16s/v3.1/utils/fun-analysis/Ttest/v3.0/result/input1.txt
    --groupfile <file>           样本分组文件，含表头。两列，第一列样本名称，第二列分组名称。脚本根据分组文件中的样本进行绘图。 示例: /home/dongxj/work/develop/16s/v3.1/utils/fun-analysis/Ttest/v3.0/result/sample_group2.txt
    --markfile <file>            标记文件，两列内容，第一列要添加的第二层分类名称，第二列要添加的标记 例如 * **, 含表头。  示例：/home/dongxj/work/develop/16s/v3.1/utils/fun-analysis/Ttest/v3.0/result/mark1.txt  [default: NA]      
    --categorize1 <string>       PDF显示时，指定第一层显示的名称 [default: group]
    --categorize2 <string>       PDF显示时，指定第二层显示的名称 [default: OTU]
    --ylab <string>              barplot图y轴标签 [default: OTU Values]
    --outputfile <file>          输出pdf文件                                                      
    --rlib <string>              R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library/]"-> doc

opts         <- docopt(doc, version = 'PICRUSt2差异柱状图\n')
inputfile    <- opts$inputfile
groupfile    <- opts$groupfile
markfile     <- opts$markfile 
categorize1  <- opts$categorize1
categorize2  <- opts$categorize2
ylab         <- opts$ylab
outputfile   <- opts$outputfile
rlib         <- opts$rlib
.libPaths(rlib)

library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)

##测试数据
# inputfile = "/home/dongxj/work/develop/16s/v3.1/utils/fun-analysis/Ttest/v3.0/result/input1.txt"
# groupfile = "/home/dongxj/work/develop/16s/v3.1/utils/fun-analysis/Ttest/v3.0/result/sample_group1.txt"
# markfile  = "/home/dongxj/work/develop/16s/v3.1/utils/fun-analysis/Ttest/v3.0/result/mark1.txt"
# outputfile = "/home/dongxj/work/develop/16s/v3.1/utils/fun-analysis/Ttest/v3.0/test1.pdf"
# categorize1 = "group"
# categorize2 = "OTU"
# ylab = "test"

# inputfile = "/home/dongxj/work/develop/16s/v3.1/utils/fun-analysis/Ttest/v3.0/result/input2.txt"
# groupfile = "/home/dongxj/work/develop/16s/v3.1/utils/fun-analysis/Ttest/v3.0/result/sample_group2.txt"
# markfile  = "/home/dongxj/work/develop/16s/v3.1/utils/fun-analysis/Ttest/v3.0/result/mark2.txt"
# outputfile = "/home/dongxj/work/develop/16s/v3.1/utils/fun-analysis/Ttest/v3.0/test2.pdf"
# categorize1 = "group"
# categorize2 = "OTU"
# ylab = "test"

##读取数据
input = read.table(inputfile, header = TRUE, check.name = F, quote = "", comment.char = "", sep = "\t", fill = T,  stringsAsFactors = FALSE)
colnames(input)[1:2] = c("group_level1", "group_level2") 

sampleInfo = read.table(groupfile, header = T, sep = "\t", quote="", colClasses = "character")
colnames(sampleInfo) = c("Sample", "Group")
rownames(sampleInfo) = sampleInfo[,1]
group = as.character(sampleInfo[,2])
samplegroup = unique(sampleInfo$Group)

##检查样本
lost_samples = sampleInfo$Sample[!sampleInfo$Sample %in% colnames(input)]
if (length(lost_samples) > 0){
    message("[Error] 部分分组文件中的样本在input文件中没有找到 : ", lost_samples)
    q()
}

##检查第一层同一组别的数据是否连续排列
level1 = unique(input$group_level1)
data <- list()
for (i in level1){
    groupnum = max(which(input$group_level1 == i)) - min(which(input$group_level1 == i)) +1 
    if (length(which(input$group_level1 == i)) != groupnum) {
       message("[WARNING] 第一层分类，相同组别没有连续排列，排序后进行后续绘图 : ", i )
    }
    data = rbind(data, input[which(input$group_level1 == i),])
}

##准备绘图数据
group_data = data[,c("group_level1", "group_level2", rownames(sampleInfo)), drop = F]
group_data$rank = 1:nrow(group_data)   ## 根据group_level2编号，第一、二层分类 标签的位置由此确定
group_data = gather(group_data, group, value, -group_level1, -group_level2, -rank) ##五列数据：group_level1 group_level2 rank group value 
for (i in 1:nrow(sampleInfo)){
    group_data$group = gsub(sampleInfo$Sample[i], sampleInfo$Group[i], group_data$group) ##分组名称替换样本名称
}

group_data$group_level1 = factor(group_data$group_level1, levels = unique(group_data$group_level1))
group_data$group_level2 = factor(group_data$group_level2 , levels = unique(group_data$group_level2))

##PDF设置
pdf(outputfile, width = 10, height = 10)

##第一张图，第一层分类绘制（分组名字，竖线）
text_x = sapply(level1, function(g) {median(group_data[which(as.character(group_data$group_level1) == g), "rank"])}) ## 根据rank列，得到每组group的位置（取中位数）
text_y = rep(0.14, length(level1 ))
seg_x1 = sapply(level1, function(g) {min(group_data[which(as.character(group_data$group_level1) == g), "rank"])-0.4}) ##这里的最小、最大值对应下面第三张图y轴bar的坐标，线段需要在这个位置上下扩展0.5bar长度（这里定的0.4，所以第三张图柱子的宽度设为0.8），否则一个样本一个分组的情况只会展示一个点
seg_x2 = sapply(level1, function(g) {max(group_data[which(as.character(group_data$group_level1) == g), "rank"])+0.4})
seg_y1 = rep(0.17, length(level1))
seg_y2 = rep(0.17, length(level1))
text_xlim   = c(1,length(unique(group_data$group_level2))) ##获取text_xlim，即rank范围。绘图时固定x轴取值范围

p1 = ggplot() + 
     annotate("text", x = text_x, y = text_y, label = level1) + ##annotate函数，通过设置X和Y的值来控制添加文本的具体位置; "text" 表示添加文本
     annotate("segment", x = seg_x1, xend = seg_x2, y = seg_y1, yend = seg_y2, size = 1) +  ##添加线段，一共4个值，分别对应 x起始， x终止，y起始，y终止; "segment" 表示添加线段
     theme_bw() + ##ggplot2 自带主题，theme_grey()为默认主题，theme_bw()为白色背景主题，theme_classic()为经典主题
     theme(axis.text = element_text(colour = "white"), axis.ticks = element_line(colour = "white"), panel.border = element_rect(colour = "white"), panel.grid = element_blank(), axis.title.y = element_blank()) +
     scale_x_continuous(expand = c(0.05, 0.05)) + ##expand = c(0.05, 0.05)，柱状图与x轴起始、终止位置均保留0.05的距离；另外两个图也需同样的限制，确保三幅图的位置不会错位；如果x是数字，则添加scale_x_continuous(); 如果x是字符/因素，则添加scale_x_discrete()
     ylab(categorize1) +
     coord_flip(clip = "off", xlim = text_xlim) + ##coord_flip()横向转换坐标：x轴和y轴互换; clip = "off" 可以画在xlim限制的范围外,因为线段坐标均向外扩展了0.4，一定会超出text_xlim
     ylim(0.12, 0.18)

##第二张图，第二层分类绘制（名字）
text2_x = unique(group_data$rank)
text2_y = rep(0.15, length(unique(group_data$rank)))
text2_label = unique(as.character(group_data$group_level2))

p2 = ggplot() + 
     annotate("text", x = text2_x, y = text2_y, label = text2_label) +
     coord_flip() +
     theme_bw() + 
     theme(axis.text = element_text(colour = "white"), axis.ticks = element_line(colour = "white"), panel.border = element_rect(colour = "white"), panel.grid = element_blank(), axis.title.y = element_blank()) +
     scale_x_continuous(expand = c(0.05, 0.05), limits = text_xlim) +
     ylab(categorize2) +
     ylim(0.12,0.18)

## 第三张图，ggbarplot
group_data$group = as.factor(group_data$group)
y.position  = group_data %>% group_by(group_level2, group) %>% summarise( value_75 = quantile(value, .75)[[1]]) %>% group_by(group_level2) %>% summarise(max_value_75=max(value_75)) ## 确定显著性标记添加的位置
y.position$max_value_75 =y.position$max_value_75 + max(y.position$max_value_75) * 0.04
p3 = ggbarplot(group_data, x = "group_level2", y = "value",
                fill = "group", palette = "jco", ##palette填充的颜色
                add = "mean_sd", add.params = list(group = "group"),
                position = position_dodge(0.7), orientation = "horiz") +
                theme_bw()+ 
                theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
                scale_x_discrete(expand = c(0.05, 0.05)) + 
                scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02*max(y.position$max_value_75))) +   ##expand = c(0, 0)，柱状图与y轴起始位置重叠，不留间隙
                ylab(ylab)

## 添加显著性标记
if(markfile != 'NA')
{ 
    mark = read.table(markfile, header = TRUE, check.name = F, quote = "", comment.char = "", sep = "\t", fill = T,  stringsAsFactors = FALSE)
    colnames(mark) = c("group_level2","p.signif") 
    y.position = data.frame(y.position) ##转换格式
    rownames(y.position) = y.position$group_level2
    add_info = data.frame(mark, .y.= rep("value",nrow(mark)), group1=rep(unique(group_data$group)[1],nrow(mark)), group2=rep(unique(group_data$group)[length(unique(group_data$group))],nrow(mark)), y.position=y.position[mark[, 1], 'max_value_75'], check.names=F)
    p3 = p3 + stat_pvalue_manual(add_info,
                            x = "group_level2",  xmin = 'group1', xmax = 'group2', y.position = 'y.position',
                            label = "p.signif")                     
}

## 组合    
grid.arrange(p1, p2, p3, ncol = 3, nrow = 1, widths = c(0.35*max(nchar(as.character(group_data$group_level1))), 0.35*max(nchar(as.character(group_data$group_level2))), 6))  ##一个字符占0.35width
dev.off()


