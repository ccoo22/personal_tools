#!/home/genesky/software/r/3.5.1/bin/Rscript

.libPaths("/home/genesky/software/r/3.5.1/lib64/R/library/")
library(docopt)
"Usage: plot_plink_pedigrees.r --ped_file <file> --output_prefix <pdf>  [--pdf_width <float> --pdf_height <float> --symbolsize <float> --cex <float>] 
Options:
   --ped_file <file>           plink ped格式家系信息文件，至少包含6列信息，不要加表头
   --pdf_width <float>         PDF图像宽度 [default: 7]
   --pdf_height <float>        PDF图像高度 [default: 7]
   --symbolsize <float>        图像中，样本形状相对于cex的大小 [default: 1]
   --cex <float>               图像中，所有形状、文字大小 （--symbolsize 1.5 --cex 0.5 效果 = 缩小文字）[default: 1]
   --output_prefix <pdf>       pdf输出文件前缀， 示例： ./pedigree " -> doc

opts                     <- docopt(doc, version = 'Program : plot pca v1.0 \n          甘斌 129\n')
ped_file                 <- opts$ped_file
output_prefix            <- opts$output_prefix
pdf_width                <- as.numeric(opts$pdf_width)
pdf_height               <- as.numeric(opts$pdf_height)
symbolsize               <- as.numeric(opts$symbolsize)
cex                      <- as.numeric(opts$cex)
 
library(FamAgg)
# ped_file="./test.ped"
ped_data = read.table(ped_file, head = F, check.names = F, stringsAsFactors=F, sep = "", fill=T)
ped_data = ped_data[, 1:6]
colnames(ped_data) = c('family', 'id', 'father', 'mother', 'sex', 'affected')


# ped数据初始化，符合kinship2格式
format_data <- ped_data
format_data$father[ped_data$father == '0' | ped_data$father == '-9'] = ''
format_data$mother[ped_data$mother == '0' | ped_data$mother == '-9'] = ''
format_data$sex[ped_data$sex == '0' | ped_data$sex == '-9'] = 3
format_data$affected[ped_data$affected == '0' | ped_data$affected == '-9'] = NA
format_data$affected[ped_data$affected == '1'] = 0
format_data$affected[ped_data$affected == '2'] = 1

affected <- format_data$affected
names(affected) <- format_data$id

# 创建数据
fad <- FAData(pedigree = format_data, trait = affected) 


# 绘图
switchPlotfun("ks2paint")
for(family in unique(format_data$family))
{   
    file = paste0(output_prefix, ".", family, ".pdf")
    message('plot : ', file)
    pdf(file, width = pdf_width, height = pdf_height)

    # 左右两侧留出空白（原图优化的不好，当样本ID字符较长时，被遮挡一部分）
    layout(matrix(c(1:3), 1, 3), widths=c(1, 10, 1))  # 按照 1 ： 10 ： 1 的比例关系分配两侧空白区域

    # 空白
    par(mar = c(5,0, 1, 0))
    plot(1:10,xlim = c(1,10), ylim=c(0,10), type = "n", xaxs = "i", yaxs = "i", axes = F, xlab="", ylab="") # 创造一个空图

    # 家系图
    par(mar = c(5,0, 1, 0))
    plotPed(fad, family = family, device = "plot", symbolsize = symbolsize, cex = cex) 

    # 空白
    par(mar = c(5,0, 1, 0))
    plot(1:10,xlim = c(1,10), ylim=c(0,10), type = "n", xaxs = "i", yaxs = "i", axes = F, xlab="", ylab="") # 创造一个空图
    dev.off()

}



 
