#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: show_circos_color.r  -i <file> -o <pdf>
Options:
   -i, --input_file <file>                    circos 颜色模版输入文件，例如：/home/ganb/soft/circos-0.69-5/etc/colors.brewer.conf
   -o, --output_pdf <file>                    结果输出pdf文件" -> doc

opts                   <- docopt(doc, version = '基于methylKit的RRBS数据差异分析,位点水平 \n')
input_file             <- opts$input_file
output_pdf             <- opts$output_pdf
 
 
# 测试
# input_file = '/home/ganb/soft/circos-0.69-5/etc/colors.brewer.conf'
# output_pdf = "colors.brewer.pdf"

# 读入数据
data = read.table(input_file, head = F, stringsAsFactors = F)
colnames(data) = c('name', 'split', 'RGB')

# 提取含有RGB信息的行
data_choose = data[grep(pattern="[0-9]{1,},[0-9]{1,},[0-9]{1,}", x=data[,3], value=FALSE), ]
data_choose$RGB16 = sapply(data_choose$RGB, function(x){ 
   color_255 <- as.numeric(unlist(strsplit(x, ','))) 
   color_16  <- rgb(color_255[1], color_255[2], color_255[3], maxColorValue=255)
   color_16
   } )

# 绘制
pdf(output_pdf, width = 30, height = 15)
# 颜色总数
color_count         = nrow(data_choose)
# 每页显示数量
color_set_each_page = 100

# 坐标模版：颜色x轴，文字x轴，y轴
x_color_tmp <- rep(seq(from = 2, by = 2, length.out = 10), each = 10)
x_txt_tmp   <- rep(seq(from = 1, by = 2, length.out = 10), each = 10)
y_tmp <- rep(c(10:1), 10)
for(page in 1:ceiling(color_count / color_set_each_page))
{  
   start = (page - 1) * color_set_each_page + 1
   end   = page * color_set_each_page
   if(end > color_count) end   = color_count 

   # 要展示的颜色
   color_show = data_choose[, 'RGB16'][start:end]
   # 要展示的文字
   text_chow  = data_choose[, 'name'][start:end]
   # 当前要展示的数量
   count_show = 1:length(color_show)

   # 10 * 10 展示
   plot(x = c(1,10), y = c(1,10), xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n', type = 'n', ann = F, xlim = c(0,21), ylim=c(0,11))
   rect(xleft = x_color_tmp[count_show] - 0.5, xright = x_color_tmp[count_show] + 0.5, ybottom = y_tmp[count_show] - 0.5 , ytop = y_tmp[count_show] + 0.5, col = color_show, border = 1)
   text(x_txt_tmp[count_show], y_tmp[count_show], text_chow)

    
}
dev.off()
