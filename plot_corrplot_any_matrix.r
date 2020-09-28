#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: plot_corrplot_any_matrix.r  -i <file> -o <pdf file> [--show_color_legend --is_corr --width <int> --shape <string>  --show_type <string> --show_coef_value <color> --item_color <string>   --height <int> --title <string> --rlib <dir>]

Options:
    -i, --input <file>              输入文件，一个完整的 m * n 矩阵，第一行、第一列都是名称。务必保证矩阵的完整性。
    -o, --output <pdf file>         输出pdf文件路径,示例：./a.pdf
    --shape <string>                每一个cell的形状。circle, square, ellipse, number, shade, color, pie [default: circle] 
    --show_type <string>            cell显示方式，全部显示，还是只显示右上角，还是左下角。 full/lower/upper [default: full]
    --show_coef_value <color>       在矩阵中显示相关性数值，同时设置颜色 [default: black]
    --item_color <string>           行、列名称颜色 [default: black]
    --is_corr                       输入的文件是不是相关性值矩阵？ 如果是，热图会按照[-1, 1]控制颜色梯度，如果不是，则以实际矩阵的数值设定梯度。 
    --show_color_legend             显示颜色示例
    --width <int>                   pdf宽度 [default: 7]
    --height <int>                  pdf高度 [default: 7]
    --title <string>                标题    [default: NA]
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，相关性图\n')
input              <- opts$input
output             <- opts$output
width              <- as.integer(opts$width)
height             <- as.integer(opts$height)
title              <- opts$title
shape              <- opts$shape
show_type          <- opts$show_type
show_coef_value    <- opts$show_coef_value
item_color         <- opts$item_color
is_corr            <- opts$is_corr
show_color_legend  <- opts$show_color_legend
rlib               <- opts$rlib
 
# 加载R包
.libPaths(rlib)
library(corrplot)

# 读入数据
message('读入数据')
data       <- read.table(input, sep='\t', header = TRUE, row.names = 1, check.names=FALSE) # 读取数据，行为样本

# 颜色
    ##  different color series
     col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                                "cyan", "#007FFF", "blue","#00007F"))
     col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                "#4393C3", "#2166AC", "#053061"))
     col3 <- colorRampPalette(c("blue", "white", "red"))
     col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                                "cyan", "#007FFF", "blue", "#00007F"))


# corrplot 绘图
if(title == 'NA') title = " "

pdf(file=output, width=width, height=height)

corrplot(as.matrix(data), 
        method = shape, 
        type = show_type, 
        title = title, 
        addCoef.col = show_coef_value, 
        tl.col = item_color, 
        is.corr = is_corr,
        cl.pos  = show_color_legend,
        col = col3(200)
        )
dev.off()
# method: 图形 "circle", "square", "ellipse", "number", "shade", "color", "pie"
# type： 显示方式，全显示还是半边 "full", "lower", "upper"
# title: 最上方的标题
# outline：外边框颜色 black/red
# addgrid.col: 网格的颜色
# addCoef.col: 显示相关性数据值，并设置颜色
# order: 对特征进行排序，选择方法： "original" "AOE" "FPC" "hclust" "alphabet"
# tl.col: 字体颜色
# p.mat: pvalue矩阵
# sig.level: pvalue阈值，大于该阈值，图形上会打一个 错误符号
# pch.col： 错误符号的颜色
# col: 相关性数值颜色梯度
# is.corr ： 输入矩阵是不是相关性值矩阵？ TRUE 按照[-1, 1] 设置颜色梯度。FALSE 则按照实际数据变动设置颜色梯度。





# ggcorrplot 绘图，个人感觉不太美观
# ggcorrplot(corr, method = 'square') # 方形,用颜色区分;circle 圆形；除了颜色、还有大小区分
# 参数
#     corr: 相关性矩阵 the correlation matrix to visualize 

#   method: 圆形、方形 character, the visualization method of correlation matrix to
#           be used. Allowed values are "square" (default), "circle".

#     type: 显示方式 全部、右下部分、左上部分 character, "full" (default), "lower" or "upper" display.

#  ggtheme: ggplot主题 ，例如： ggtheme = ggplot2::theme_gray； ggplot2 function or theme object. Default value is
#           `theme_minimal`. Allowed values are the official ggplot2
#           themes including theme_gray, theme_bw, theme_minimal,
#           theme_classic, theme_void, .... Theme objects are also
#           allowed (e.g., `theme_classic()`).

#    title: character, title of the graph.

# show.legend: logical, if TRUE the legend is displayed.

# legend.title: a character string for the legend title. lower
#           triangular, upper triangular or full matrix.

# show.diag: logical, whether display the correlation coefficients on the
#           principal diagonal.

#   colors: 颜色梯度，长度为3的颜色向量。a vector of 3 colors for low, mid and high correlation
#           values. 例如： colors = c("#6D9EC1", "white", "#E46726")

# outline.color: 外边框颜色 the outline color of square or circle. Default value is
#           "gray".

# hc.order: 是否使用聚类方法进行排序 logical value. If TRUE, correlation matrix will be hc.ordered
#           using hclust function.

# hc.method: 聚类方法 the agglomeration method to be used in hclust (see ?hclust).
# lab: 是否把相关性值显示出来 logical value. If TRUE, add correlation coefficient on the
#           plot.

# lab_col, lab_size: size and color to be used for the correlation
#           coefficient labels. used when lab = TRUE.
#    p.mat: pvalue矩阵输入。matrix of p-value. If NULL, arguments sig.level, insig, pch,
#           pch.col, pch.cex is invalid.
# sig.level: pvalue阈值，大于改阈值的认为是不显著，标记特殊符号。significant level, if the p-value in p-mat is bigger than
#           sig.level, then the corresponding correlation coefficient is
#           regarded as insignificant.

#    insig: 不显著的点标记特殊图形来源。 character, specialized insignificant correlation
#           coefficients, "pch" (default), "blank". If "blank", wipe away
#           the corresponding glyphs; if "pch", add characters (see pch
#           for details) on corresponding glyphs.
#           blank值：挖掉对应的模块

#      pch: 不显著的点标记特殊图形 add character on the glyphs of insignificant correlation
#           coefficients (only valid when insig is "pch"). Default value
#           is 4.
