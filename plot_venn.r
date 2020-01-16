#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: veen.r  -i <file> -o <string> [--Rlib <dir>]
Options:
   -i, --input_file <file>                 输入文件，一列是一组，包含表头，允许不同的列长度不同
   -o, --output_prefix <string>            结果输出前缀
                                           注意：根据项目不同，最好单独优化

   --Rlib <dir>                            R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts                   <- docopt(doc, version = 'veen图绘制 \n')
input_file             <- opts$input_file
output_prefix          <- opts$output_prefix
Rlib                   <- opts$Rlib
.libPaths(Rlib)

# 加载
library(VennDiagram)


# 读取数据
message("读入数据")
data <- read.table(input_file, header = T, sep="\t", stringsAsFactors=F, check.names = F)

# 转换成列表
message("数据转换成列表，为绘图做准备")
data_plot <- lapply(colnames(data), function(x){ 
  # 去掉NA
  data_tmp = data[, as.character(x)] # as.character防止列名是纯数字
  data_tmp = data_tmp[complete.cases(x)]
  # 去掉空白 
  data_tmp = data_tmp[data_tmp != ""]
  data_tmp
  }  )
names(data_plot) = colnames(data) # 列表命名
class_count = length(data_plot) # 数量

# 颜色模版
mycol <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62) 
mycol <- colors()[mycol] 

message("绘制")
# VennDiagram 会产生log文件，下述代码能够抑制
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
# 甘斌：根据不同维度，分别进行了优化
# 绘制一维韦恩图
if(class_count == 1)
{ 
  # 读入数据为list类型，默认图片保存为tiff格式
  venn.plot <- venn.diagram(
    x = data_plot,
    filename = paste0(output_prefix, '.png'), imagetype = "png",
    col = "black", #边框颜色
    lwd = 5, #边框线宽度
    # fontface = "bold", #标签字体
    fill = "grey", #填充色
    alpha = 0.7, #透明度
    cex = 4, #标签字体大小
    cat.cex = 3, #类名字体大小
    # cat.fontface = "bold", #类名字体
    margin = 0.04 #边际距离
  )
} else if(class_count == 2)
{ # 绘制二维韦恩图
  venn.plot <- venn.diagram(
    x = data_plot,
    filename = paste0(output_prefix, '.png'), imagetype = "png",
    lwd = 3,
    fill = c("cornflowerblue", "darkorchid1"),
    alpha = 0.6,
    label.col = "white",
    cex = 1.5,
    fontfamily = "serif",
    # fontface = "bold",
    cat.col = c("cornflowerblue", "darkorchid1"),
    cat.cex = 2,
    cat.fontfamily = "serif",
    # cat.fontface = "bold",
    margin = 0.05,
    cat.dist = c(0.03, 0.03), # 调整文字与圆圈的距离
    cat.pos = c(-20, 20) # 控制文字在x轴的相对位置，0是在圆的正上方
  )
} else if(class_count == 3)
{ # 绘制三维韦恩图
  venn.plot <- venn.diagram(
    x = data_plot,
    filename = paste0(output_prefix, '.png'), imagetype = "png",
    col = "transparent",
    fill = c("red", "blue", "green"),
    alpha = 0.5,
    # label.col = c("darkred", "white", "darkblue", "white",
   #              "white", "white", "darkgreen"),
    cex = 2.5,
    fontfamily = "serif",
    # fontface = "bold",
    # cat.default.pos = "text",
    cat.col = c("darkred", "darkblue", "darkgreen"),
    cat.cex = 2.7,
    cat.fontfamily = "serif",
    # cat.fontface = "bold", # 类名字体
    margin = 0.07, # 边际距离
    # cat.dist = c(0.06, 0.06, 0.03),
  )
} else if(class_count == 4)
{ # 绘制四维韦恩图
  venn.plot <- venn.diagram(
    x = data_plot,
    filename = paste0(output_prefix, '.png'), imagetype = "png",
    col = "black",
    lty = "dotted", #边框线型改为"dotted"虚线
    lwd = 3, # 边框线的宽度
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
    alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4", "white", "white", "white",
                "white", "white", "darkblue", "white",
                "white", "white", "white", "darkgreen", "white"),
    cex = 2.0,
    fontfamily = "serif",
    # fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
    cat.cex = 1.8,
    # cat.fontface = "bold",
    cat.fontfamily = "serif"
  )
} else if(class_count == 5)
{ # 绘制五维韦恩图
  venn.plot <- venn.diagram(
    x = data_plot,
    filename = paste0(output_prefix, '.png'), imagetype = "png",
    lty = "dotted",
    lwd = 2,
    col = "black",  #"transparent",
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    alpha = 0.60,
    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.cex = 0.8,
    # cat.fontface = "bold",
    margin = 0.07,
    cex = 0.8
  )
} else
{ # 绘制n维韦恩图
  venn.plot <- venn.diagram(
    x = data_plot,
    filename = paste0(output_prefix, '.png'), imagetype = "png",
    # col = "black", #边框颜色
    lwd = 3, #边框线宽度
    fill = mycol[1:length(data_plot)], #填充色
    alpha = 0.6, #透明度
    cex = 1.5, # 标签字体大小
    fontfamily = "serif",
    # fontface = "bold", #标签字体
    cat.col = mycol[1:length(data_plot)],# 类名颜色
    cat.cex = 1.7, # 类名字体大小
    cat.fontfamily = "serif",
    # cat.fontface = "bold", # 类名字体
    margin = 0.05, # 边际距离
  ) 
}

