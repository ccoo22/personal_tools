#!/home/genesky/software/r/4.0.3/bin/Rscript

library(docopt)

"Usage: corr_heatmap.r  --table_file <file> --sample_group <file> --case <string> --control <string> --output_dir <pdf file> [--stop_cluster_rows --stop_cluster_columns --method <string> --width <int>  --height <int> --heatmap_proportion <float> --rlib <dir>]

Options:
    --table_file <file>            绘制多组学相关性热图的必要文件。含表头，每一行表示一个组学相关信息，第一列组学名称，第二列表达量文件，第三列标签文件。
                                   表达量文件：含表头，每一行是一个特征，每一列是一个样本，第一行是样本名，第一列是特征名。
                                   标签文件：含表头，第一列特征名称，其余列是标签列legend1、legend2、legend3..，热图图例根据标签文件列名添加，不同组学的标签文件若存在相同列名只能显示其中一个图例；若标签文件不存在，table_file相应位置不填写即可
    --sample_group <file>          样本分组文件。含表头，第一列样本名称，第二列分组名称，流程对--case --control指定组的样本进行后续相关性计算及绘图
    --case <string>                指定case组名称，即--sample_group第二列的某一分组名称
    --control <string>             指定control组名称，即--sample_group第二列的某一分组名称
    --output_dir <file>            输出路径
    --method <string>              相关性计算方法，可选 'pearson', 'spearman', 'kendall' [default: spearman]
    --width <int>                  PDF宽度， 默认自适应调整（如果出图不好看，可以直接指定手工调整）
    --height <int>                 PDF高度， 默认自适应调整
    --heatmap_proportion <float>   热图的主体部分占据整张图的比例 [default: 0.9]
                                   因为特征的名称可能非常长，会被pdf边缘覆盖，此时可以适当降低这个数值，以保证边缘留出更大的空白
    --stop_cluster_rows            不做行聚类
    --stop_cluster_columns         不做列聚类
    --rlib <dir>                   R包路径 [default: /home/genesky/software/r/4.0.3/lib64/R/library]" -> doc

opts           <- docopt(doc, version='ComplexHeatmap包 -- 绘制相关性热图\n')
table_file     <-  opts$table_file
sample_group   <-  opts$sample_group
case           <-  opts$case
control        <-  opts$control
output_dir     <-  opts$output_dir
method         <-  opts$method
stop_cluster_rows    <-  opts$stop_cluster_rows
stop_cluster_columns <-  opts$stop_cluster_columns
width          <-  as.integer(opts$width)
height         <-  as.integer(opts$height)
heatmap_proportion <-  as.numeric(opts$heatmap_proportion)
rlib           <-  opts$rlib

# 加载R包
.libPaths(rlib)
library(psych) ##corr.test
library(ggpubr)
library(ComplexHeatmap)
set.seed(123)

# 测试数据
# table_file = "/home/dongxj/work/research/plot/heatmap/test/table_file2.txt"
# sample_group = "/home/dongxj/work/research/plot/heatmap/test/sample_group.txt"
# case    = "Post"
# control = "Baseline"
# method   = 'spearman'


# 读入数据
table_data  <- read.table(table_file,  sep='\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
sample_data <- read.table(sample_group, sep='\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
colnames(table_data)  <-  c('exp_txt', 'legend_txt')
colnames(sample_data) <-  'group'

# 检查样本, 读取数据
cases <- rownames(sample_data[sample_data[,'group']==case,,drop=F])
controls <- rownames(sample_data[sample_data[,'group']==control,,drop=F])
samples <- c(cases, controls)

data_name_exp <- list()
data_name_legend <- list()
for (i in 1:nrow(table_data)) {
    data_name <- rownames(table_data)[i]
    data_exp  <- read.table(table_data[,'exp_txt'][i],  sep='\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
    data_samples = colnames(data_exp)

    lost_samples = samples[!samples %in% data_samples]
    if (length(lost_samples) > 0)
    {
        message("[Error] case control分组文件中的样本在", data_name, "表达量文件中没有找到 : ", lost_samples)
        q()
    }

    ## 表达量数据，对--case --control指定组的样本进行后续相关性计算
    data          <- data_exp[,samples]
    data_name_exp[[i]] <- data

    ## legend数据，计算log2foldchange，判断每个特征在哪一组富集
    if (table_data[,'legend_txt'][i] != ''){##判断legend文件是否存在
        data_legend    <- read.table(table_data[,'legend_txt'][i],  sep='\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
        legend_feature <- rownames(data_legend)
        exp_feature    <- rownames(data)
        lost_features  <- exp_feature[!exp_feature %in% legend_feature]
        if (length(lost_features) > 0){ ##检查标签文件的特征是否有缺失
            message("[Error]", data_name,"标签文件中的部分特征在其表达量文件中没有找到 :", lost_features)
            q()   
        }else{
            data_log2FoldChange <- as.data.frame(log(apply(data[,cases], 1, mean)/apply(data[,controls], 1, mean),2))
            colnames(data_log2FoldChange) <- 'log2FoldChange'
            data_legend <- data_legend[rownames(data),,drop=F] ##legend文件与exp文件特征名称顺序保持一致
            data_legend$enrich_group[data_log2FoldChange$log2FoldChange > 0] <- paste('enriched in', case, 'samples')
            data_legend$enrich_group[data_log2FoldChange$log2FoldChange < 0] <- paste('enriched in', control, 'samples')    
            data_name_legend[[i]] <- data_legend  
        }
    }else{
        message("[Warning]",rownames(table_data)[i],"标签文件为空") 
        data_log2FoldChange <- as.data.frame(log(apply(data[,cases], 1, mean)/apply(data[,controls], 1, mean),2))
        colnames(data_log2FoldChange) <- 'log2FoldChange'
        data_legend <- data_log2FoldChange[rownames(data),,drop=F] ##legend文件与exp文件特征名称顺序保持一致，如果legend文件缺失，则legend只有一列log2Foldchange信息，后续根据该列得到enrich_group信息
        data_legend$enrich_group[data_log2FoldChange$log2FoldChange > 0] <- paste('enriched in', case, 'samples')
        data_legend$enrich_group[data_log2FoldChange$log2FoldChange < 0] <- paste('enriched in', control, 'samples')    
        data_name_legend[[i]] <- data_legend[,'enrich_group',drop=F]  
    }
}

# 两两计算相关性值，并绘制热图
for (i in 1:(nrow(table_data)-1)) {
    for (j in (i+1):nrow(table_data) ){
        data1_name <- rownames(table_data)[i]
        data2_name <- rownames(table_data)[j]
        data1_exp    <- data_name_exp[[i]]
        data2_exp    <- data_name_exp[[j]]
        data1_legend <- data_name_legend[[i]]
        data2_legend <- data_name_legend[[j]]

        # 计算相关性值、pvalue值
        cor <- corr.test(t(data1_exp), t(data2_exp), method = method, adjust = "none") ### corr.test()进行特征间相关性分析（按列数据进行相关性分析，文件内容需要转置）;如果不矫正，即adjust ="none"，其相关系数、P值和cor.test()结果一样。
        cmt <- cor$r  
        pmt <- cor$p
        
        # 输出 相关性值、pvalue值
        corr_file     <- paste(output_dir, '/',data1_name, '_', data2_name, '_corr.txt', sep = '')
        output_corr   <- data.frame(id = rownames(cmt), cmt, check.names = F)
        write.table(output_corr, corr_file,  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        pvalue_file   <- paste(output_dir, '/',data1_name, '_', data2_name, '_pvalue.txt', sep = '')
        output_pvalue <- data.frame(id = rownames(pmt), pmt, check.names = F)
        write.table(output_pvalue, pvalue_file,  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")          

        # 显著性P值表示方式 '*' '**''
        p_mark = matrix('', nrow(pmt), ncol(pmt))
        p_mark[pmt<0.01] <- '**'
        p_mark[pmt>0.01 & pmt <0.05] <- '*'
        
        #### 绘图注释颜色设定
        # 定义相关性值颜色
        cor_color <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))

        # enrich组富集颜色设定
        enrich_group_name  <- unique(data1_legend$enrich_group)
        enrich_group_color <- c('#cf5641','#90900a')
        names(enrich_group_color) <- enrich_group_name

        # 注释列颜色
        legend_color <- c('ucscgb', 'aaas', 'lancet', 'jco', 'npg', 'uchicago', 'simpsons', 'rickandmorty')

        # data1 legend
        data1_color <- list()
        if(ncol(data1_legend) < 2){
            data1_color$enrich_group <- enrich_group_color #legend只有一列enrich_group信息
        }else{
            for (a in 1:(ncol(data1_legend)-1)){
                legend <- unique(data1_legend[,a])
                colors = get_palette(palette = legend_color[a], k = length(legend) )
                names(colors) <- legend
                data1_color[[a]] <- colors
            }
            names(data1_color) <- colnames(data1_legend)[1:(ncol(data1_legend)-1)]
            # data1 所有legend
            data1_color$enrich_group <- enrich_group_color
        }

        # data2 legend
        data2_color <- list()
        if(ncol(data2_legend) < 2){
            data2_color$enrich_group <- enrich_group_color #legend只有一列enrich_group信息
        }else{
            for (b in 1:(ncol(data2_legend)-1)){
                legend <- unique(data2_legend[,b])
                colors = get_palette(palette = legend_color[b+4], k = length(legend) )
                names(colors) <- legend
                data2_color[[b]] <- colors
            }
            names(data2_color) <- colnames(data2_legend)[1:(ncol(data2_legend)-1)]
            # data2 所有legend
            data2_color$enrich_group <- enrich_group_color
        }

        #### 行列注释内容
        ha_row     <- HeatmapAnnotation(df = data1_legend, which = 'row', col = data1_color, show_annotation_name = FALSE, show_legend = TRUE, annotation_name_gp = gpar(fontsize = 7), gap = unit(0.8, "points"))
        ha_column  <- HeatmapAnnotation(df = data2_legend, which = 'col', col = data2_color, show_annotation_name = FALSE, show_legend = TRUE, annotation_name_gp = gpar(fontsize = 7), gap = unit(0.8, "points"))

        ####显著性差异标记legend
        lgd_sig = Legend(pch = c("*", "**"), type = "points", labels = c("P < 0.05", "P < 0.01"))

        ####pdf宽度自适应模块
        if(length(width) == 0)
        {
            width = ncol(cmt) / 2  # 1 inches -> 2个特征
            height = nrow(cmt) / 4 # 纵向没有legend展示，要求不高
            # 如果太小，则优先使用12 inches
            width = ifelse(width < 12, 12, width)
            height = ifelse(height < 12, 12, height)
            # 如果太大，则设置为最大值
        }
        # 绘图
        output_pdf <- paste(output_dir, '/', data1_name, '_', data2_name, '_corr_heatmap.pdf', sep='')
        cairo_pdf(output_pdf, width = width, height = height, family="GB1")  # inches
        par(oma=c(2,2,2,2),mar=c(2,2,2,2))

        a = Heatmap(cmt,
            name = 'Spearman\'s p', ##相关性值的legend 名称(legend_title)
            border = NA, ##cell间不加边框
            cluster_rows = !stop_cluster_rows, ##行聚类
            show_row_dend = FALSE, ##展示行聚类树
            cluster_columns = !stop_cluster_columns, #列聚类
            show_column_dend = FALSE, ##展示列聚类树
            row_names_side = 'left', ##行名在左侧
            column_names_side = 'top', ##列名在右侧
            row_names_max_width = unit(15, "cm"), ##设置行名的最大宽度，避免行名太长显示不全
            column_names_max_height = unit(15, "cm"),
            heatmap_width = unit(heatmap_proportion * width, "inches"), ##整体宽度，包括热图、行名、图例
            heatmap_height = unit(heatmap_proportion * height, "inches"),##整体高度，包括热图、列名、图例
            row_names_gp = gpar(fontsize = 9),
            column_names_gp = gpar(fontsize = 8),
            top_annotation = ha_column, ##列的注释
            left_annotation = ha_row, ##行的注释
            col = rev(cor_color(100)), ##相关性值的颜色范围
            heatmap_legend_param = list(at = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)), ##默认at = c(-1, -0.5, 0, 0.5, 1)
            layer_fun = function(j, i, x, y, width, height, fill) { ##cell填充，这里填充 * ** 
                    v = pindex(p_mark, i, j)
                    grid.text(v, x, y)})

        draw(a, merge_legend = TRUE, heatmap_legend_side = 'right', annotation_legend_list = list(lgd_sig)) ## merge_legend = TRUE, 一列展示
        dev.off()        
    }
}  

