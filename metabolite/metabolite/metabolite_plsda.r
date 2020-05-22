#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)
"Usage: MetaboAnalyst_plsda -i <file> --plsda <file> --vip <file> --plsda_coordinate <file> --plsda_vip <file> --sample_group_file <file> [--color <string> --feat_num <int> --show_confidence]
Options:
	-i, --input <file>                 the input file, each row is a feature, each column is sample
	--plsda <file>                     the output plsda pdf file
	--vip <file>                       the output plsda vip pdf file
	--plsda_coordinate <file>          plsda coordinate xls file
	--plsda_vip <file>                 plsda vip xls file
	--sample_group_file <file>         the input group file, column 1 is sampleinfo, column 2 is groupinfo
	--color <string>                   format: \'group1,red;group2,blue;group3,green\' [default: NA]
	--show_confidence                  if show confidence in picture
	--feat_num <int>                   the number of feature in vip for ploting [default: 15]" -> doc

opts              <- docopt(doc, version = 'Program : MetaboAnalyst PLSDA plot v1.2 \n      lhj 272\n')
inputfile         <- opts$input
plsdaplot         <- opts$plsda
vipplot           <- opts$vip
plsda_coordinate  <- opts$plsda_coordinate
plsda_vip         <- opts$plsda_vip
groupfile         <- opts$sample_group_file
color             <- opts$color
confidence        <- opts$show_confidence
feat_num          <- as.numeric(opts$feat_num)

# 测试用参数
# inputfile         <- '/home/lhj/research/MetaboAnalyst/PLSDA/test/input/normalized_data.txt'
# plsdaplot         <- '/home/lhj/research/MetaboAnalyst/PLSDA/test/output/PLSDA.pdf'
# vipplot           <- '/home/lhj/research/MetaboAnalyst/PLSDA/test/output/PLSDAvip.pdf'
# plsda_coordinate  <- '/home/lhj/research/MetaboAnalyst/PLSDA/test/output/PLSDA.sites.xls'
# plsda_vip         <- '/home/lhj/research/MetaboAnalyst/PLSDA/test/output/PLSDA.vip.xls'
# # groupfile         <- '/home/lhj/research/MetaboAnalyst/PLSDA/test/input/sample_group.xls'
# # color             <- '1,Red;3,Black;2,Purple;4,Gray;5,skyblue;6,orangered;7,White;8,Yellow;9,Blue'
# # groupfile         <- '/home/lhj/research/MetaboAnalyst/PLSDA/test/input/sample_group_9.xls'
# # color             <- 'B1,Red;3,Black;2,Purple;A4,Gray;5,skyblue;6,orangered;7,White;8,Yellow;9,Blue'
# # groupfile         <- '/home/lhj/research/MetaboAnalyst/PLSDA/test/input/sample_group_2.xls'
# # color             <- 'B1,Red;9,Blue'
# groupfile         <- '/home/lhj/research/MetaboAnalyst/PLSDA/test/input/sample_group_4.xls'
# color             <- 'B1,Red;9,Blue;6,orangered;A4,Gray'
# confidence        <- TRUE
# feat_num          <- 20

# options(warn=-1) #忽视任何警告
message("\nstart MetaboAnalyst_plsda.r\n")

library(mixOmics, quietly = TRUE)

# 数据前处理
data <- read.table(inputfile, sep = '\t', header = TRUE, check.names = FALSE) # 读取定量数据，列为样本
sample.groups <- read.table(groupfile, sep = "\t", header = FALSE, quote = "", stringsAsFactors = FALSE) # 读取样本以及分组数据
# PLSDA分析要求
# （1）输入数据，行为特征值，列为样本
# （2）每一列不能存在缺失数据

# 获取指定分析的样本
# 判断指定的样本是否存在于原始数据
sample_need  = sample.groups[,1]
if (length(sample_need[!(sample_need %in% colnames(data))]) != 0)
{
	cat('[Error] please check the sample_name_list !\n')
	q()
}
data_need    = data[, sample_need, drop = FALSE]
sample_group = as.character(sample.groups[,2])

# 判断分组信息是否有缺失值
if (NA %in% sample_group)
{
	cat('[Error] please check the group_name_list !\n')
	q()
}

data_new <- t(as.matrix(data_need))               # 预处理完用于后续分析的数据
colnames(data_new) <- data[,1]                    # 原始数据id号（第一列）作为列名,确保无重复

# PLS-DA analysis
plsda <- plsda(data_new, sample_group, ncomp = 2) #  选择展示前 2 个排序轴
sites <- as.data.frame(plsda$variates$X)          #  提取样本点坐标
sites$Group = factor(sample_group, levels = unique(sample_group)) # factor数据，levels与输入顺序一致，以固定plsda图上的分组顺序
colnames(sites) = c("plsda1", "plsda2", "Group")
result <- cbind(rownames(sites), sites)
colnames(result) = c("Sample", "plsda1", "plsda2", "Group")
write.table(result, file = plsda_coordinate, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

message("绘制 PLS_DA")
comp1 <- round(plsda$explained_variance$X[1]*100, 2) # comp 1主元解释百分比
comp2 <- round(plsda$explained_variance$X[2]*100, 2) # comp 2主元解释百分比

if (color != 'NA')
{
	message("使用指定颜色")
	Groups = unique(sample.groups[,2])
	group_col = unlist(strsplit(unlist(strsplit(color, ';')), ','))
	if (length(group_col) != 2 * length(Groups))
	{
		message("指定颜色数量与分组数不一致")
		q()
	}
	color_define = t(matrix(group_col, nrow = 2))
	mycol        = color_define[, 2]
	names(mycol) = color_define[, 1]  # 组名与颜色一一对应

}else{
	message("使用默认颜色设置")
	mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)  
	mycol <- colors()[rep(mycol, 50)]
}

# 形状模版
myshape <- rep(c(15,16,17,18,19,20,21,22,23,24,25,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14),4)

pdf(plsdaplot, width = 13, height = 13)

if (confidence == TRUE)
{
	message("展示95%置信区间")
	p = ggplot(data=sites, aes(x = plsda1, y = plsda2, colour = Group, shape = Group)) +  # 映射
        geom_hline(yintercept = 0, colour = "gray65") +                                   # 添加水平线
        geom_vline(xintercept = 0, colour = "gray65") +                                   # 添加垂直线
        scale_shape_manual(values = myshape) +                                            # 指定形状
        scale_color_manual(values = mycol) +                                              # 指定颜色
        geom_point(size = 3, alpha = 1) +                                                 # 设置点的大小、透明度
        stat_ellipse(level = 0.95, show.legend = F)  +                                    # 添加置信区间
        ggtitle("PLSDA plot of MetaboAnalyst") + xlab(paste("X-variate 1 ", "(", comp1, "%)", sep="")) + ylab(paste("X-variate 2 ", "(", comp2, "%)", sep="")) + theme(plot.title = element_text(hjust = 0.5))   # 定义标题、横纵坐标
}else{
	message("不展示置信区间")
	p = ggplot(data=sites, aes(x = plsda1, y = plsda2, colour = Group, shape = Group)) +  # 映射
        geom_hline(yintercept = 0, colour = "gray65") +                                   # 添加水平线
        geom_vline(xintercept = 0, colour = "gray65") +                                   # 添加垂直线
        scale_shape_manual(values = myshape) +                                            # 指定形状
        scale_color_manual(values = mycol) +                                              # 指定颜色
        geom_point(size = 3, alpha = 1) +                                                 # 设置点的大小、透明度
        ggtitle("PLSDA plot of MetaboAnalyst") + xlab(paste("X-variate 1 ", "(", comp1, "%)", sep="")) + ylab(paste("X-variate 2 ", "(", comp2, "%)", sep="")) + theme(plot.title = element_text(hjust = 0.5))   # 定义标题、横纵坐标
}

p
dev.off()


# 计算vip值
p_vip <- vip(plsda)                       # Variable Importance in the Projection (VIP)
comp1vip <- as.matrix(p_vip[,c("comp 1")])
res <- cbind(data[,1:2], comp1vip)        # 输出原始数据NO. Peak 以及对应的vip值
colnames(res) <- c("NO.", "Peak", "vip scores")
write.table(res, file = plsda_vip, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

message("绘制 VIP")

# 特征数目不足，则全部绘制
if (feat_num > nrow(res))
{
	feat_num <- nrow(res)
}

sample_group_uniq   <- unique(sample_group)

# 取前 n 个特征值展示
top.value   <- res[order(res[,"vip scores"], decreasing = T)[1:feat_num],]
imp.top.vec <- top.value[,"vip scores"] # vip
names(imp.top.vec) <- top.value[,"NO."] # 编号 NO.
feat_name   <- rev(top.value[,"Peak"])  # 特征名 Peak
imp.top.vec <- rev(imp.top.vec)         # 升序展示
# 对某特征值，数据按分组取均值 并转置
mns <- by(data_new[, names(imp.top.vec)], sample_group, function(x){ apply(x, 2, mean, trim = 0.1) })
mns <- sapply(mns, function(x){x})

pdf(vipplot, width = 13, height = 13)

layout(matrix(c(1,2,3), nrow = 1, byrow = T), widths = c(60 - (3 + 2 * length(sample_group_uniq)), 2*length(sample_group_uniq), 3))  # 画布布局
par(mar = c(5,18,2,0.1))                  # bottom, left, top, and right

vip.nms <- substr(feat_name, 1, 24)       # 仅展示特征名前 24 个字符
names(imp.top.vec) <- NULL

# 绘制vip值点图
dotcolor <- "blue"
xlbl     <- "VIP scores"
x <- imp.top.vec
y <- 1:feat_num  
plot(1:3, xlim = c(x[1], rev(x)[1]), ylim = c(0, feat_num + 1), type = "n", yaxs = "i", axes = T, yaxt="n", xlab = xlbl, ylab = "", cex.lab = 1.5)  # 坐标图
abline(h = y, col = "gray", lwd = 2, lty = 2)                         # 坐标图里加虚线
points(x, y, bty = "n", pch = 21, bg = dotcolor, cex = 2)             # 坐标图里加点
mtext(vip.nms, side = 2, at = 1:feat_num, las = 2, line = 1, cex = 1) # 左侧加特征名标签. side 1=下面，2=左边，3=上边，4=右边; line 文字距边轴的距离;las 标签是否平行于（0）或垂直于（2）坐标轴

y_min = par("usr")[3]
y_max = par("usr")[4]

# 绘制分组点图
# 对每个特征值，各组按序分配颜色，值最大的展示红色，最小为绿色
par(mar = c(5,0.3,2,0))
plot(1:3, ylim = c(y_min, y_max), type = "n", xaxs = "i", yaxs = "i", axes = F, xlab = "", ylab = "") # 空图
nc  <- ncol(mns)
col <- colorRampPalette(c("red", "green"))(nc)

bg  <- matrix("", nrow(mns), nc)
colnames(bg) <- colnames(mns)
for (m in 1:nrow(mns))
{
	colnames(mns)
	bg[m,] <- (col[nc:1])[rank(mns[m,])]
}
# colnames(bg) <- names(rank(mns[1,])) # 列名为rank后的names

shiftx <- (2.8 - 1.2)/length(sample_group_uniq) # 点间距
shifty <- (y_max - y_min)/(feat_num + 1)

# 设定点坐标
x <- rep(1.2, feat_num)
y <- 1:feat_num

for (n in 1:ncol(mns)){
    points(x, y, bty = "n", pch = 22, bg = bg[,sample_group_uniq[n]], cex = 3)            # 根据坐标、颜色 绘制点
    text(x[1], (feat_num+0.4)*shifty, sample_group_uniq[n], srt = 45, adj = c(0.2,0.5))   # 加分组名
    x <- x + shiftx                                                                       # 各组点图间隔
}

# 绘制红绿色图例
# 确定坐标
par(mar = c(5,0,2,0))
col <- colorRampPalette(c("red", "green"))(50)
nc  <- length(col)
x   <- rep(1.5, nc)  # 定横坐标
# 定纵坐标
starty <- 3 * (y_max - y_min)/6
endy   <- 4 * (y_max - y_min)/6
y <- seq(from = starty, to = endy, length = nc)

plot(1:2, ylim = c(y_min, y_max), type = "n", xaxs = "i", yaxs = "i", axes = F, xlab = "", ylab = "")
points(x, y, bty = "n", pch = 15, col = rev(col), cex = 3)

text(x[1],  endy+shifty/4, "High")                                         # 加标签
text(x[1], starty-shifty/4, "Low")

dev.off()

message("\nfinish MetaboAnalyst_plsda.r\n")