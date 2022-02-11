#!/home/genesky/software/r/4.0.3/bin/Rscript

library(docopt)
"Usage: batch_correct_combat.r --input <file> --batch <file> --batch_variables <string>  --output_prefix <string> [ --keep_variables <string> --keep_na]
Options:
    --input, -i <file>          表达量矩阵文件，每一行是一个基因、每一列是一个样本、第一行是样本，第一列是基因名。行、列名称不能有重复。
                                可同时输入多个文件（逗号分隔），脚本自动根据行名进行合并，仅保留所有文件共有的基因，同时，样本名不要有重复。
    --batch, -b <file>          批次信息文件，用于矫正的批次。第一列是样本名，第二列及之后的列是批次信息，不能有缺失值。第一行是表头。
                                仅分析该文件包含的样本
    --batch_variables <string>  需要矫正的变量名称，多个变量用逗号分隔. 
    --keep_variables <string>   需要保留的变量名称，多个变量用逗号分隔（如果你的样本中包含癌、癌旁样本，则一定要把这个信息作为 keep_variables，否则矫正过程中会把癌、癌旁这种分组信息矫正掉） [default: NA]
                                也可以不填写这个变量。 
    --keep_na                   保留表达量存在缺失值的基因，默认去除
    --output_prefix, -o <string>    输出文件的前缀，例如： ./result  。 最终生成的结果为 result.diff.xls。
" -> doc

opts                   <- docopt(doc, version = '批次效应矫正 \n          甘斌 129\n')
input            <- opts$input
batch            <- opts$batch
batch_variables  <- opts$batch_variables
keep_variables   <- opts$keep_variables
output_prefix    <- opts$output_prefix
keep_na    <- opts$keep_na
# input = 'GSE10358_series_matrix.format.txt,GSE6891_series_matrix.format.txt'
batch_variables = unlist(strsplit(batch_variables, ','))
keep_variables  = unlist(strsplit(keep_variables, ','))
if(keep_variables[1] == 'NA') keep_variables = NULL

library(sva)
library(factoextra)
library(ggpubr)
library(tidyr)
##### (1) 数据预处理

# 读入表达量矩阵
data_list = list()
gene_accum = c()
sample_accum = c()
count = 0
for(file in unlist(strsplit(input, ',')))
{   
    count = count + 1
    message('[process] 读取表达量文件: ', file)
    data_tmp = read.table(file, header=T, row.names = 1, sep='\t', check.names=F, quote='')
    gene_count = nrow(data_tmp)
    sample_count = ncol(data_tmp)
    message("          包含 ",gene_count, " 个基因，", sample_count, " 个样本")
    data_list[[file]] = data_tmp
    gene_accum = c(gene_accum, rownames(data_tmp))
    sample_accum = c(sample_accum, colnames(data_tmp))
}

# 确认共有的基因
message('[process] 查找共有的基因')
gene_common = c()
gene_appear  = as.data.frame(table(gene_accum), stringsAsFactors=F)
gene_common = gene_appear[ gene_appear$Freq == length(data_list), 1]
message('          共有的基因有: ', length(gene_common), " 个")
if(length(gene_common) < 10)
{
    message("[Error] 共有基因数量太少，不能分析")
    q()
}

# 确认样本是否有重复
message('[process] 确认是否有重复样本')
sample_appear  = as.data.frame(table(sample_accum), stringsAsFactors=F)
if(sum(sample_appear$Freq > 1) > 0)
{
    message("[Error] 输入的矩阵中存在重复的样本 ：",  paste(sample_appear[ sample_appear$Freq > 1, 1], collapse=','))
    q()
}

# 合并数据
message('[process] 表达量矩阵合并')
data_raw = data_list[[1]][gene_common, ,drop=F]
if(length(data_list) > 1)
{
    for(i in 2:length(data_list))
    {
        data_raw = data.frame(data_raw, data_list[[i]][gene_common, ,drop=F], check.names=F)
    }
}

# 读入批次数据
message('[process] 读入批次数据')
batch_info = read.table(batch, header=T, row.names=1, sep='\t', check.names=F)

if(sum(!complete.cases(batch_info))  > 0 )
{
    message("[Error] 批次信息不能有缺失值")
    q()
}

# 检查样本是否缺失
message('[process] 检查样本名是否异常')
if(sum(!rownames(batch_info) %in% colnames(data_raw)) > 0 )
{
    lost_samples = rownames(batch_info)[!rownames(batch_info) %in% colnames(data_raw)]
    message("[Error] batch_info 中的部分样本在表达文件中缺失：", paste(lost_samples, collapse =','))
    q()
}else{
    message('          样本都存在，共纳入 ', nrow(batch_info), ' 个样本')
}

# 检查矫正的表型是否有缺失
message('[process] 检查矫正的批次名称是否异常')
variables = c(batch_variables, keep_variables)
if(sum(!variables %in% colnames(batch_info)) > 0 )
{
    lost_variables = variables[!variables %in% colnames(batch_info)]
    message("[Error] batch_info 中的缺少指定的批次信息：", paste(lost_variables, collapse =','))
    q()
}

# 表达矩阵清理，仅保留需要的样本,同时按照batch信息排序
message('[process] 输出合并后的原始表达量矩阵（仅保留batch文件包含的样本），同时去掉有缺失的基因')
data_raw = data_raw[, row.names(batch_info)]
complete_condition = complete.cases(data_raw)

if(! keep_na)
{   
    message('          去除有缺失值的基因')

    data_lost = data_raw[!complete_condition, ]
    data_raw = data_raw[complete_condition, ]
    message('          去掉基因 ', nrow(data_lost), ' 个', '，剩余 ', nrow(data_raw), '个')
    if(nrow(data_lost) > 0)
    {
        file_lost = paste0(output_prefix, '.Before_combate.lost_gene.txt')
        write.table(data.frame(id=rownames(data_lost), data_lost, check.names=F), file_lost, sep='\t', row.names=F, quote=F)
    }
}else{
    message('          按照参数要求，保留有缺失值的基因')
    message('          数据中共检测到 ', sum(!complete_condition) ,' 个有缺失值的基因')
}

file_before = paste0(output_prefix, '.Before_combate.txt')
write.table(data.frame(id=rownames(data_raw), data_raw, check.names=F), file_before, sep='\t', row.names=F, quote=F)


# 绘图函数

select_top_n_samples <- function(info, variable, topn)
{ 
    # info: data.frame格式，行名为样本名
    # variable: info的某一个列名称
    # topn: 整数
    # 根据info矩阵的variable列信息，提取样本名称，相同variable值下，最多取topn个
    sample_topn = c()
    for(variable_value in unique(info[, variable]))
    {
        sample_count = sum(info[, variable] == variable_value)
        sample_count = ifelse(sample_count > topn, topn, sample_count)
        sample_topn = c(sample_topn, rownames(info)[info[, variable] == variable_value][1:sample_count] )
    }
    return(sample_topn)
}

data_plot <- function(data_input, batch_info, variables, output_prefix, title)
{   
    surfix = gsub(' ','_', title)

    data_input = data_input[complete.cases(data_input), ]   # 去除有缺失的基因
    data_input <- data_input[rowSums(data_input) != 0, ]  # 去除表达量为0的基因
    pca_obj <- prcomp(t(data_input))
    pc1_proportion <- summary(pca_obj)$importance[2, 1]*100
    pc2_proportion <- summary(pca_obj)$importance[2, 2]*100
    for(variable in variables)
    {   
        message('          boxplot 绘图，legend 批次信息为：', variable)
        # boxplot 图
        # 每个批次仅展示最多20个样本，不然就太多了
        samples_boxplot = select_top_n_samples(batch_info, variable, 20)
        data_boxplot = gather(data_input[, samples_boxplot], sample_id, expression)  # 坍塌
        data_boxplot$batch = batch_info[data_boxplot$sample_id, variable]  # 加样本分组
        data_boxplot$batch = factor(data_boxplot$batch, levels = unique(data_boxplot$batch))  # batch 也设为factor
        data_boxplot$sample_id = factor(data_boxplot$sample_id, levels = unique(data_boxplot$sample_id))  # 固定样本顺序: 注意一定要放在最后，否则会影响上一步 加样本分组的正确执行
        
        # 绘图
        pdf_boxplot = paste0(output_prefix, '.', surfix, '.boxplot.', variable, '.pdf')
        pdf(pdf_boxplot, width=ifelse(length(samples_boxplot) * 0.14 > 7, length(samples_boxplot) * 0.14, 7), height=7)
        p1 <- ggboxplot(data_boxplot, "sample_id", "expression",
                color = "black", 
                fill = "batch",
                outlier.shape = NA,
                palette = 'ucscgb',
                error.plot = 'errorbar',
                legend = "right",
                legend.title = variable,
                title = title
                )
        p1  <- p1 + theme(axis.ticks.x=element_blank(), axis.text.x = element_blank())  # 去掉x轴刻度、文本
        print(p1)
        dev.off()
    
        # PCA 图
        # factoextra 包绘图
        message('          PCA 绘图，legend 批次信息为：', variable)
        pdf_raw_pca = paste0(output_prefix, '.', surfix, '.pca.', variable, '.pdf')
        pdf(pdf_raw_pca)
        p2 <- fviz_pca_ind(pca_obj,
                                     geom.ind = c("point"), # 图中的样本点以何种形式展示，包含三种参数"text"以文本形式展示
                                                                                               # "point"以点形式展示
                                                                                               # "arrow"以箭头形式展示
                                                                                               # 三种形式可以单独展示，也可以任意两两组合进行展示，三种同时使用也是可以的
    
                                     col.ind = batch_info[, variable],   # 样本颜色根据分组信息来区分
    
                                     palette = "ucscgb",       # 包含20中颜色的调色版，基于ggsci包，同样的支持ggsci中的其它调色板
                                                                                               # 比如"aaas","jco","ucscgb"等
    
                                     addEllipses = TRUE,           # 在图中添加分组椭圆
                                     ellipse.type = "t",           # 设置分组椭圆的类型，可选择的类型有:"convex","confidence","t","norm","euclid"
                                     ellipse.level = 0.66          # 设置分组椭圆置信区间
                                     )
        p2 <- p2 +  scale_shape_manual(values=0:18)  # 设定形状候选列表
        # 调用ggpubr包调整图片标题
        p2 <- ggpar(p2,
            title = title,
            xlab = paste("PC1 [",pc1_proportion,"%]", sep = ""),
            ylab = paste("PC2 [",pc2_proportion,"%]", sep = ""),
            legend.title = variable, 
            legend.position = "right",
            ggtheme = theme_bw()
              )
        print(p2)
        dev.off()
    }

}

##### (2) 原始数据绘图
message('[process] 原始数据绘图，查看是否存在批次效应')
# 每一个批次单独出图，方便查看是否存在批次问题
# PCA / boxplot 绘图
data_plot(data_raw, batch_info, batch_variables, output_prefix, 'Before combate')


##### (3) 批次效应矫正
message('[process] 批次效应矫正')
data_combat = data_raw
for(i in 1:length(batch_variables))
{
    message("          矫正批次 ",batch_variables[i])
    # 本次批次矫正，需要保留的变量有
    mod = NULL
    mod_variables = c(keep_variables)
    if(i+1 <= length(batch_variables))
    {
        mod_variables = c( mod_variables,  batch_variables[(i+1):length(batch_variables)])
    }
    if(length(mod_variables) > 0)
    {
        formdf = as.formula(paste(" ~ ",paste(mod_variables, collapse=" + "),sep=""))
        mod <- model.matrix(formdf, data=batch_info)
    }

    data_combat <- ComBat(dat = data_combat, batch = batch_info[, batch_variables[i]], mod = mod)
}
data_combat = as.data.frame(data_combat)

##### (4) 矫正结果输出
message('[process] 输出矫正后的表达量矩阵')
file_after = paste0(output_prefix, '.After_combate.txt')
write.table(data.frame(id=rownames(data_combat), data_combat, check.names=F), file_after, sep='\t', row.names=F, quote=F)


##### (5) 矫正后的结果绘图
message('[process] 矫正后的数据绘图，查看是否存在批次效应')
data_plot(data_combat, batch_info, batch_variables, output_prefix, 'After combate')
