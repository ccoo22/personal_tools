#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: diff_two_data_match_gene.r --class <file> --class_format <string> --exp <file> --exp_format <string> -o <file> -p <pdf> [ --class_order <string> --sample_file <file> --feature_file <file> --pdf_y_lab <string> --pdf_width <numeric> --pdf_height <numeric>  ]
Options:
    --class <file>                  样本分类数据矩阵。每一行对应一个特征，每一列对应一个样本，第一行是样本名，第一列是特征名。 注：也可以转置一下。具体的数据类型通过 class_format 声明
                                    这个数据，可以想象为一个矩阵，记录了每一个基因在每一个样本上发生拷贝数变化与否。离散型数值。
    --class_format <string>         class文件格式 col/row。col表示class文件每一列对应一个样本；row表示class文件每一行对应一个样本 [default: col]
    --exp <file>                    样本表达量数据矩阵。每一行对应一个特征，每一列对应一个样本，第一行是样本名，第一列是特征名。 注：也可以转置一下。具体的数据类型通过 exp_format 声明
                                    这个数据，可以想象为一个矩阵，记录了每一个基因在每一个样本上的表达量。连续性数值。
    --exp_format <string>           exp文件格式 col/row。col表示class文件每一列对应一个样本；row表示class文件每一行对应一个样本 [default: col]
    -o, --output <file>             pvalue结果文件
    -p, --pdf_file <pdf>            绘图文件pdf

    --class_order <string>          class文件中，数据在绘图时横坐标显示的顺序。 如果class每一个特征的分类都不一样，不建议使用。如果都一样，可以指定顺序。这里的值可以比class中的多。多个value用逗号分隔。如果提供该参数，务必保证包含所有种类的value.
    --sample_file <file>            样本文件，一列数据，含有表头，表示使用class/exp中哪些样本进行分析。
                                    默认使用class文件中的所有样本
    --feature_file <file>           特征文件，一列数据，含有表头，表示使用class/exp中哪些特征进行分析。
                                    默认使用class文件中的所有特征进行分析
    --pdf_width <numeric>           PDF宽度 [default: 7]
    --pdf_height <numeric>          PDF高度 [default: 7]
    --pdf_y_lab <string>            绘制boxplot图时，y轴名称 [default: exp] 
" -> doc

opts              <- docopt(doc, version = '两个文件中，对相同特征的两个数据做差异分析（例如：同一个基因的CNV 、表达量）， 差异分析用anova \n')
class             <- opts$class
class_format      <- opts$class_format
exp               <- opts$exp
exp_format        <- opts$exp_format
output            <- opts$output
class_order       <- opts$class_order
pdf_file          <- opts$pdf_file
pdf_width         <- as.numeric(opts$pdf_width)
pdf_height        <- as.numeric(opts$pdf_height)
pdf_y_lab         <-  opts$pdf_y_lab 
sample_file       <-  opts$sample_file 
feature_file      <-  opts$feature_file 


library(ggpubr)

# 读入数据 统一改为 每一行是样本
message("read class data")
data_class = read.table(class, header = T, row.names = 1, check.names = F, sep = '\t')
if(class_format == 'col') data_class = t(data_class)

message("read exp data")
data_exp = read.table(exp, header = T, row.names = 1, check.names = F, sep = '\t')
if(exp_format == 'col') data_exp = t(data_exp)

# 确定分析的样本与特征
samples = rownames(data_class)
features = colnames(data_class)
if(!is.null(sample_file))
{
    data_sample   = read.table(sample_file, header = T, check.names = F, sep = "\t", stringsAsFactors=F, colClasses = 'character') 
    samples       = data_sample[, 1]  
}
if(!is.null(feature_file))
{
    data_feature   = read.table(feature_file, header = T, check.names = F, sep = "\t", stringsAsFactors=F, colClasses = 'character') 
    features       = data_feature[, 1]  
}

# 检查class/exp中样本、特征是否一致
# 样本
if(sum(!samples %in% rownames(data_class)) > 0)
{
    losts = samples[!samples %in% rownames(data_class)]
    message("[Error] class 文件中样本缺失:", losts)
    message("请仔细检查class文件，以及class_format参数是否符合你的数据格式")
    q()
}
if(sum(!samples %in% rownames(data_exp)) > 0)
{
    losts = samples[!samples %in% rownames(data_exp)]
    message("[Error] exp 文件中样本缺失:", losts)
    message("请仔细检查exp文件，以及exp_format参数是否符合你的数据格式")
    q()
}
# 特征
if(sum(!features %in% colnames(data_class)) > 0)
{
    losts = features[!features %in% colnames(data_class)]
    message("[Error] class 文件中特征缺失:", losts)
    message("请仔细检查class文件，以及class_format参数是否符合你的数据格式")
    q()
}
if(sum(!features %in% colnames(data_exp)) > 0)
{
    losts = features[!features %in% colnames(data_exp)]
    message("[Error] exp 文件中特征缺失:", losts)
    message("请仔细检查exp文件，以及exp_format参数是否符合你的数据格式")
    q()
}
# class绘图显示顺序确认
if(!is.null(class_order))
{   
    message("确定分组显示顺序")
    class_order = unlist(strsplit(class_order, ','))
    for(feature in features)
    {
        class = unique(data_class[, feature])
        class = class[!is.na(class)]
        if(sum(!class %in% class_order) > 0)
        {   
            lost = class[!class %in% class_order]
            message("[Error] class_order 没有包含所有的特征：", lost)
            q()
        }
    }
}


# 开始差异分析与绘图
pdf(pdf_file, width = pdf_width, height = pdf_height)

result = matrix(NA, nrow = length(features), ncol = 3)
colnames(result) = c('feature', 'sample count', 'pvalue')
row = 0
for(feature in features)
{   
    message("[process] feature ", feature)
    row = row + 1
    x = data_class[samples, feature]
    y = data_exp[samples, feature]
    analysis_data <- data.frame(x=x, y=y)
    analysis_data <- analysis_data[complete.cases(analysis_data),]
    x <- analysis_data$x
    y <- analysis_data$y

    pvalue = NA
    sample_count = nrow(analysis_data)
    # 去掉缺失值后，只有1组数据了，不能计算
    if(length(unique(x)) > 1)
    {
        fit <- aov(y~x)
        result_all       <- summary(fit) # 提取p
        pvalue <- result_all[[1]]['x', 'Pr(>F)']
    }
    result[row, ] = c(feature, sample_count, pvalue)

    # 绘图
    # 设定显示顺序
    if(!is.null(class_order))
    {   
        classes = class_order[class_order %in% levels(analysis_data$x)]
        analysis_data$x = factor(as.character(analysis_data$x), levels = classes)
    }
    pvalue_show = NA
    if(!is.na(pvalue)) pvalue_show = format(pvalue, scientific = TRUE, digits=4)
    p <- ggboxplot(analysis_data, x="x", y="y", color = "x",
          palette = "aaas", #杂志Science的配色
          add = 'jitter', # 取均值，并添加标准差曲线
          add.params = list(color = 'x'),
          outlier.size = -1, # 隐藏离群点，否则离群点会在图上出现两次（boxplot绘制离群点，jitter再绘制所有点，从而使离群点绘制两次）
          title = feature,
          subtitle  = paste0('pvalue = ', pvalue_show)
          ) + ylab(pdf_y_lab) + guides(fill = FALSE, color = FALSE) + theme(plot.title = element_text(hjust = 0.5))  # 去掉legend,标题居中
    print(p)
}

write.table(result, output, row.names = F, quote = F, sep = '\t')
dev.off()