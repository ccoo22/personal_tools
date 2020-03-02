# 检测 -> 脚本输入
ARGS <- commandArgs(T)
if(length(ARGS) != 3){
	cat("
Program: example pca
Version: v1.0
Contact: 297 赵若玮

Usage:   Rscript", SCRIPT, "ARGS1 ARGS2 ARGS3

Options:

         INPUTFILE     输入文件
		 GROUP         分组信息
         OUTDIR        结果文件输出目录
	\n");
	q()
}
INPUTFILE <- normalizePath(ARGS[1])
GROUP <- normalizePath(ARGS[2])
OUTDIR <- normalizePath(ARGS[3])

# 导入 -> package
library("factoextra")
library("ggpubr")

###################################################################### 主程序

# 读取分组信息
sampleinfo <- read.table(GROUP, header = FALSE, sep = "\t") # 读取分组信息
colnames(sampleinfo) <- c("sample", "group") # 分组信息列表列名赋值
rownames(sampleinfo) <- sampleinfo[,1] # 分组信息列表行名赋值

# 读取原始文件
inputdata <- read.table(INPUTFILE, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, comment.char = "", quote = "")

# 数据过滤
# 获取用于分析的分组样本数据,去除多余数据（如：OTU表后的注释，其他组数据等）
newdata <- inputdata[,rownames(sampleinfo)]
# 数据过滤，去除所有样本中相加和为0的数据
newdata <- newdata[rowSums(newdata) != 0, ]

# 数据专制成行为样本，列为变量
pcadata <- t(newdata)

#PCA 分析
# 对转置后的数据行进行求和，即求每个样本中的所有变量之和
num <- apply(pcadata, 1, sum)

# 判断是否将数据进行标准化	
if(length(unique(num)) == 1){
	apca <- prcomp(pcadata)
}else{
	# 去除常量列
	good <- sapply(1:ncol(pcadata), function(x){
		if(length(unique(pcadata[,x])) == 1){
			return(FALSE)
		}
		return(TRUE)
	})
	pcadata <- pcadata[,good]
	apca <- prcomp(pcadata, scale = TRUE)
}


# 输出结果
# 获取两个主成分得分
pc1 <- apca$x[ ,1]
pc2 <- apca$x[ ,2]
sites <- cbind(pc1, pc2)
colnames(sites) <- c("PC1", "PC2")
write.csv(sites, paste(OUTDIR, "/pc_sites.csv", sep=""), quote=FALSE)
	
# 获取两个主成分贡献值
pc1_proportion <- summary(apca)$importance[2,1]*100
pc2_proportion <- summary(apca)$importance[2,2]*100
write.csv(summary(apca)$importance, paste(OUTDIR, "/pc_cont.csv", sep=""), quote=FALSE)

# 画图
pdf(paste(OUTDIR, "/pca.pdf", sep=""),width = 10, height = 10)

# factoextra包绘图
p <- fviz_pca_ind(apca,
				 geom.ind = c("text","point"), # 图中的样本点以何种形式展示，包含三种参数"text"以文本形式展示
											   # "point"以点形式展示
											   # "arrow"以箭头形式展示
											   # 三种形式可以单独展示，也可以任意两两组合进行展示，三种同时使用也是可以的
											   
				 col.ind = sampleinfo$group,   # 样本颜色根据分组信息来区分
				 
				 palette = "category20",       # 包含20中颜色的调色版，基于ggsci包，同样的支持ggsci中的其它调色板
											   # 比如"aaas","jco","ucscgb"等
											   
				 addEllipses = TRUE,           # 在图中添加分组椭圆
				 ellipse.type = "t",           # 设置分组椭圆的类型，可选择的类型有:"convex","confidence","t","norm","euclid"
				 ellipse.level = 0.66          # 设置分组椭圆置信区间

				 )

# 调用ggpubr包调整图片标题
ggpar(p,
      title = "Principal Component Analysis",
      subtitle = "Iris data set",
      caption = "Source: factoextra",
      xlab = paste("PC1 [",pc1_proportion,"%]", sep = ""),
	  ylab = paste("PC2 [",pc2_proportion,"%]", sep = ""),
      legend.title = "Groups", legend.position = "top",
      ggtheme = theme_bw()
	  )
	  
dev.off()