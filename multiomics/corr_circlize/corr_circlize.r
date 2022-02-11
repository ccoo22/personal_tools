#!/home/genesky/software/r/4.0.3/bin/Rscript

library(docopt)

"Usage: corr_circlize.r  --table_file <file> --sample_group <file> --case <string> --control <string> --output_dir <pdf file> [--method <string> --width <int>  --height <int> --filter_pvalue <numeric> --filter_cor <numeric> --rho_color_value <numeric> --rlib <dir>]

Options:
    --table_file <file>            绘制多组学相关性热图的必要文件。含表头，每一行表示一个组学相关信息，第一列组学名称，第二列表达量文件。
                                   表达量文件：含表头，每一行是一个特征，每一列是一个样本，第一行是样本名，第一列是特征名。
    --sample_group <file>          样本分组文件。含表头，第一列样本名称，第二列分组名称，流程对--case --control指定组的样本进行后续相关性计算及绘图
    --case <string>                指定case组名称，即--sample_group第二列的某一分组名称
    --control <string>             指定control组名称，即--sample_group第二列的某一分组名称
    --output_dir <file>            输出路径
    --method <string>              相关性计算方法，可选 'pearson', 'spearman', 'kendall' [default: spearman]
    --filter_pvalue <numeric>      pvalue值过滤 [default: 0.05]
    --filter_cor <numeric>         相关性值（+表示正相关，-表示负相关）过滤条件，绘图时link线条颜色，也由相关性值决定，大于0线条为蓝色，表示正相关，小于0线条为红色，表示负相关 [default: 0.8]
    --width <int>                  PDF宽度 [default: 10]
    --height <int>                 PDF高度 [default: 10]
    --rlib <dir>                   R包路径 [default: /home/genesky/software/r/4.0.3/lib64/R/library]" -> doc

opts           <-  docopt(doc, version='circlize包 -- 绘制相关性 弦图\n')
table_file     <-  opts$table_file
sample_group   <-  opts$sample_group
case           <-  opts$case
control        <-  opts$control
output_dir     <-  opts$output_dir
method         <-  opts$method
filter_pvalue  <- as.numeric(opts$filter_pvalue)
filter_cor     <- as.numeric(opts$filter_cor)
width          <-  as.integer(opts$width)
height         <-  as.integer(opts$height)  
rlib           <-  opts$rlib

# 加载R包
.libPaths(rlib)
library(psych) ##corr.test
library(circlize) ##circos.par()
library(dplyr)  ##distinct
library(tidyr) ##max min
set.seed(123)

# 测试数据
# table_file = "/home/dongxj/work/research/plot/circlize/data/test/table_file.txt"
# sample_group = "/home/dongxj/work/research/plot/circlize/data/test/sample_group.txt"
# case    = "Post"
# control = "Baseline"
# method  = 'spearman'
# output_dir = "/home/dongxj/work/research/plot/circlize/data/test/result"

# 读入数据
table_data  <- read.table(table_file,  sep='\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
sample_data <- read.table(sample_group, sep='\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
colnames(table_data)[1]  <-  'exp_txt'
colnames(sample_data) <-  'group'

# 读取指定分组的样本
cases <- rownames(sample_data[sample_data[,'group']==case,,drop=F])
controls <- rownames(sample_data[sample_data[,'group']==control,,drop=F])
samples <- c(cases, controls)

# 读取表达量数据
data_case_exp    <- list()
data_control_exp <- list()
data_all_exp         <- list()
for (i in 1:nrow(table_data)) {
    data_name <- rownames(table_data)[i]
    data_exp  <- read.table(table_data[,'exp_txt'][i],  sep='\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
    data_samples <- colnames(data_exp)

    lost_samples = samples[!samples %in% data_samples]
    if (length(lost_samples) > 0)
    {
        message("[Error] case control分组文件中的样本在", data_name, "表达量文件中没有找到 : ", lost_samples)
        q()
    }
    ## 表达量数据，对--case --control指定组的样本及所有样本分别进行后续相关性计算
    data_case_exp[[i]]    <- data_exp[,cases]
    data_control_exp[[i]] <- data_exp[,controls]
    data_all_exp[[i]]     <- data_exp[,samples]
}


# case control分别 两两计算相关性值、pvalue值, 并绘图
plot_corr_circlize <- function(group_name, sample_name, exp_data){
    circlize_all_links <- data.frame()
    for (i in 1:(nrow(table_data)-1)) {
        for (j in (i+1):nrow(table_data) ){
            data1_name <- rownames(table_data)[i]
            data2_name <- rownames(table_data)[j]
            data1_exp    <- exp_data[[i]]
            data2_exp    <- exp_data[[j]]

            # 计算相关性值、pvalue值，分别输出
            result = matrix(nrow = nrow(data1_exp) * nrow(data2_exp), ncol = 7);       
            result_estimate = matrix(nrow = nrow(data1_exp), ncol = nrow(data2_exp))
            result_pvalue   = matrix(nrow = nrow(data1_exp), ncol = nrow(data2_exp))
            
            rownames(result_estimate) = rownames(data1_exp)
            colnames(result_estimate) = rownames(data2_exp)
            
            rownames(result_pvalue) = rownames(data1_exp)
            colnames(result_pvalue) = rownames(data2_exp)

            row = 0
            for(data1 in rownames(data1_exp))
            {
                for(data2 in rownames(data2_exp))
                {
                    row  <- row + 1
                    data <- cbind(name1 = as.numeric(data1_exp[data1, sample_name]), name2 = as.numeric(data2_exp[data2, sample_name]))
                    data <- data.frame(data[complete.cases(data),]) #去掉缺失值
                    data <- data[apply(data == "", 1, sum) == 0, ]  #去掉空值               
                    nmiss <- nrow(data) #nmiss必须>1
                    pmt <- NA
                    cmt <- NA
                    if(nmiss > 2){
                        cor <- corr.test(data[,1], data[,2], method = method, adjust = "none") ### corr.test()进行特征间相关性分析（按列数据进行相关性分析，文件内容需要转置）;如果不矫正，即adjust ="none"，其相关系数、P值和cor.test()结果一样。      
                        cmt <- cor$r
                        pmt <- cor$p
                    }
                    result[row, ] <- c(data1_name, data1, data2_name, data2, nmiss, pmt, cmt)
                    result_estimate[data1, data2] <- cmt
                    result_pvalue[data1, data2] <- pmt
                }
            }
            # 输出 相关性值、pvalue值
            result_estimate <- data.frame(id = rownames(result_estimate), result_estimate, check.names=F)
            result_pvalue   <- data.frame(id = rownames(result_pvalue), result_pvalue, check.names=F)
            corr_file       <- paste(output_dir, '/', group_name, '_', data1_name, '_', data2_name, '_corr.txt', sep = '')
            write.table(result_estimate, corr_file,  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
            pvalue_file     <- paste(output_dir, '/', group_name, '_', data1_name, '_', data2_name, '_pvalue.txt', sep = '')
            write.table(result_pvalue, pvalue_file,  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

            ##所有的link
            circlize_all_links <- rbind(result, circlize_all_links)
        }
    }
    circlize_all_links = circlize_all_links[complete.cases(circlize_all_links), ]  # 一定要先清理掉含有NA的记录
    colnames(circlize_all_links) <- c('source_class', 'source', 'target_class', 'target',  'NMISS', 'Pvalue', 'Estimate')
    circlize_all_links$Pvalue   <- as.numeric(circlize_all_links$Pvalue)
    circlize_all_links$Estimate <- as.numeric(circlize_all_links$Estimate)
    circlize_filter_links <-  circlize_all_links[circlize_all_links$Pvalue < filter_pvalue & abs(circlize_all_links$Estimate) > filter_cor, ]
    circlize_filter_links <-  circlize_filter_links[order(abs(circlize_filter_links$Estimate),decreasing=T),] #相关性值降序
    circlize_filter_links$direction[circlize_filter_links$Estimate > 0] <- 'plus correlation'
    circlize_filter_links$direction[circlize_filter_links$Estimate < 0] <- 'minus correlation'
    circlize_filter_links <- circlize_filter_links[,c('source_class', 'source', 'target_class', 'target', 'Pvalue', 'Estimate','direction')]
    ##输出绘图数据
    circlize_filter_links_file  <- paste(output_dir, '/', group_name, '_circlize.links', sep = '')
    write.table(circlize_filter_links, circlize_filter_links_file,  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    ##排序（去重复），绘图指定绘图顺序
    orderdata1 <- circlize_filter_links[,c('source_class','source')]
    colnames(orderdata1) <- c('group','description')
    orderdata2 <- circlize_filter_links[,c('target_class','target')]
    colnames(orderdata2) <- c('group','description')
    orderdata <- distinct(rbind(orderdata1,orderdata2))  ##dplyr包distinct函数 去重复，最终两列数据 group description
 
    ##同一组别的数据连续排列
    group <- unique(orderdata$group)
    order_data <- list()
    for (i in group){
        order_data <- rbind(order_data, orderdata[which(orderdata$group == i),])
    }
 
    ##统计group列每个类别的数量，从1开始计数，用于绘图时的位置坐标
    ##添加绘图时类别名称label字体颜色（每个类别的颜色相同）
    group_count = table(order_data[,'group']) ## table()函数统计每个类别的数量
    color = c('#5a5ba6', '#ea4a46', '#a8883c', '#c53f2d', '#26479b') ##每个类别标签的颜色（暂定5个类别的颜色）
    for(i in 1:length(group_count))
    {   
        group = names(group_count)[i]
        order_data[order_data$group == group, 'pos'] = 1:group_count[i]
        order_data[order_data$group == group, 'color'] = rep(color[i], group_count[i])
    }

    ## 自定义染色体（每个类别视为一个染色体）
    df = data.frame()
    for(i in 1:length(group_count))
    {   
        df = rbind(df, c(names(group_count)[i], 0, group_count[i]))
    }
    colnames(df) = c('chr', 'start', 'end') ##三列数据

    ##开始绘图
    outputfile  <- paste(output_dir, '/', group_name,'_corr_circlize.pdf', sep = '')
    cairo_pdf(outputfile, width = width, height = height, family = 'GB1')

    ## 初始化
    circos.par(track.height = 0.1, gap.degree = 0.7, start.degree = 90, cell.padding=c(0.02, 0, 0.02, 0), track.margin = c(-0.2, 0.2))  #gap.degree = 0.7不同类别染色体之间的间隙, start.degree = 90开始绘图位置（角度）, circos.par() 设置全局变量，所有轨道的宽度都是0.1半径（完整半径等于1），cell.padding=c(0.02, 0, 0.02, 0)绘图区域四周边距，track.margin = c(-0.2, 0.2)轨道边距参数
    circos.genomicInitialize(df, plotType = NULL) # 初始化，声明扇形（自定义的染色体（每个类别视为一个染色体））

    ##绘制轨道及添加文字
    circos.track(factors = order_data$group, x = order_data$pos, y = order_data$pos, bg.border = F, ## bg.border = F 不显示扇形边框
        panel.fun = function(x, y) {
            labels = order_data[order_data$group == CELL_META$sector.index, 'description']
            color = order_data[order_data$group == CELL_META$sector.index, 'color']
            circos.text(x, rep(0, length(x)), labels, cex = 0.5, col = color, facing = "reverse.clockwise", adj=1, niceFacing = TRUE) ##文字右对齐adj=1, 左对齐0，居中显示0.5, niceFacing = TRUE 最合适的角度展示文字（方便查看）
    }) 

    # 添加连线
    bed1 = data.frame()
    bed2 = data.frame()
    for (i in 1:nrow(circlize_filter_links))
    {
        position1 = order_data[order_data[,'description'] == circlize_filter_links[i, 'source'], 'pos']
        link_color = ifelse(circlize_filter_links[i, 'direction'] == 'plus correlation', '#4c91c1', '#ff0c0c') 
        bed1 = rbind(bed1, c(circlize_filter_links[i, 'source_class'], position1, position1, link_color))

        position2 = order_data[order_data[,'description'] == circlize_filter_links[i, 'target'], 'pos']
        bed2 = rbind(bed2, c(circlize_filter_links[i, 'target_class'], position2, position2, link_color))
    }
    colnames(bed1) = c('chr', 'start', 'end', 'color')
    colnames(bed2) = c('chr', 'start', 'end', 'color')
    bed1$start = as.integer(bed1$start)
    bed1$end = as.integer(bed1$end)
    bed2$start = as.integer(bed2$start)
    bed2$end = as.integer(bed2$end)

    circos.genomicLink(bed1, bed2, col = bed1$color, border = NA)
    circos.clear()
    dev.off()
}

##case control 所有样本分别执行
plot_corr_circlize(case, cases, data_case_exp)
plot_corr_circlize(control, controls, data_control_exp)
plot_corr_circlize('All', samples, data_all_exp)
