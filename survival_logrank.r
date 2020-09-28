#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: survival_logrank.r --surv <file> --class <file> -o <file> -p <pdf> [--class_format <string> --pdf_width <numeric> --pdf_height <numeric> --pdf_surv_median_line <string> --pdf_conf_int --pdf_pvalue --pdf_fun <string> --pdf_risk_table --pdf_add_all --pdf_palette <string> --pdf_xlab <string> --pdf_legend <string> --pdf_legend_x <numeric> --pdf_legend_y <numeric> --pdf_break_x_by <numeric> --pdf_legend_title <string>]
Options:
    -s, --surv <file>               生存信息矩阵，第一列样本名，第二列time, 第三列 status（0 alive, 1 dead）， 有表头
    --class <file>                  logrank分析，样本分类矩阵。每一行对应一个特征，每一列对应一个样本，第一行是样本名，第一列是特征名称。不允许缺失。如果有缺失，会自动排除包含缺失的样本。
                                    注意1：流程只会分析--surv文件中包含的样本
                                    注意2：另外，输入的文件也可以是每一行对应一个样本，每一列对应一个特征。数据格式请通过
    --class_format <string>         class文件格式 col/row。col表示exp文件每一列对应一个样本；row表示exp文件每一行对应一个样本 [default: col]
    -o, --output <file>             pvalue结果文件
    -p, --pdf_file <pdf>            统计结果绘图文件pdf
    --pdf_pvalue                    增加pvalue， 建议添加
    --pdf_risk_table                图像下面添加随时间变化，样本生存数量表， 建议添加

    --pdf_width <numeric>           PDF宽度 [default: 7]
    --pdf_height <numeric>          PDF高度 [default: 7]
    --pdf_surv_median_line <string> 图中添加50%生存时间线，包括水平、垂直两种线. none, hv, h, v [default: hv] 
    --pdf_conf_int                  增加置信区间
    --pdf_fun <string>              生存分析绘图形式 'cumhaz' plots the cumulative hazard function; 'pct' for survival probability in percentage。 默认为NULL, 与pct等同，只是纵坐标*100.
    --pdf_add_all                   添加一条所有样本放在一起的生存曲线（不会影响pvalue）
    --pdf_palette <string>          样本分组调色板 hue, grey, npg , aaas , lancet , jco ,  ucscgb , uchicago , simpsons 和 rickandmorty [default: npg]
    --pdf_xlab <string>             X坐标轴名称 [default: Time]
    --pdf_legend <string>           图例显示位置 top, bottom, left, right, none [default: top]
    --pdf_legend_x <numeric>        图例显示位置，x相对坐标轴, 取值范围 [0, 1]， 例如 0.8。 人工指定legend显示位置，如果声明了pdf_legend_x与pdf_legend_y两个参数，则忽略pdf_legend参数。
    --pdf_legend_y <numeric>        图例显示位置，y相对坐标轴, 取值范围 [0, 1]， 例如 0.75。 (0.8, 0.75) 坐标位于图像内部的右上角
    --pdf_break_x_by <numeric>      对x轴刻度时间点设定显示间隔， 例如 ： 100， 则以100为时间间隔显示x轴刻度
    --pdf_legend_title <string>     图例标题。默认是特征的名字。
    --pdf_title <string>            图片标题。
" -> doc

opts               <- docopt(doc, version = 'survival  log-rank生存分析软件 \n')
surv              <- opts$surv
class             <- opts$class
class_format      <- opts$class_format
output            <- opts$output
pdf_file          <- opts$pdf_file
pdf_width         <- as.numeric(opts$pdf_width)
pdf_height        <- as.numeric(opts$pdf_height)
pdf_surv_median_line <-  opts$pdf_surv_median_line 
pdf_conf_int         <-  opts$pdf_conf_int 
pdf_pvalue           <-  opts$pdf_pvalue
pdf_fun              <-  opts$pdf_fun 
pdf_risk_table       <-  opts$pdf_risk_table 
pdf_add_all          <-  opts$pdf_add_all 
pdf_palette          <-  opts$pdf_palette 
pdf_xlab             <-  opts$pdf_xlab 
pdf_legend           <-  opts$pdf_legend 
pdf_legend_x         <-  opts$pdf_legend_x 
pdf_legend_y         <-  opts$pdf_legend_y 
pdf_break_x_by       <-  opts$pdf_break_x_by
pdf_legend_title     <-  opts$pdf_legend_title 
pdf_title            <-  opts$pdf_title 


library(survival)
library(survminer)

# 读入数据
message("read surv data")
surv_data = read.table(surv, head = T, row.names = 1, check.names = F, sep = '\t')
surv_data = surv_data[, c(1,2)]
colnames(surv_data)[1:2] = c('time', 'status') 
samples = rownames(surv_data)

message("read class data")
data_class   = read.table(class, head = T, row.names = 1, check.names = F, sep = "\t")  # 分类表型信息
if(class_format == 'col') data_class   = t(data_class)


if(sum(!samples %in% rownames(data_class)) > 0)
{
    losts = samples[!samples %in% rownames(data_class)]
    message("[Error] surv 中的样本在 class 文件中缺失:", losts)
    message("请仔细检查class文件，以及class_format参数是否符合你的数据格式")
    q()
}
data_class_clean = data_class[samples, , drop=F] 

# 开始logrank分析
count = 0
pvalue_data = matrix(0,nrow = ncol(data_class_clean), ncol=3)
colnames(pvalue_data) = c('feature', 'sample_count', 'Pvalue')

pdf(pdf_file, width = pdf_width, height = pdf_height)
for(feature in colnames(data_class_clean))
{
    message("[process Log-Rank] ", feature)
    count = count + 1
    data_analysis = cbind(surv_data, data_class_clean[, feature])
    colnames(data_analysis)[3] = 'variable' # 由于列名存在‘-’字符，不支持，所以表型统一命名为variable,ggsurvplot中添加legend.title

    data_analysis <- data_analysis[complete.cases(data_analysis), ]
    data_analysis <- data_analysis[apply(data_analysis, 1, function(x){sum(x=='')}) == 0, ]
    sample_count  = nrow(data_analysis)
    # 只有一组数据，无法分析
    if(length(unique(data_analysis[, 3])) == 1 )
    {
        pvalue_data[count, ] = c(feature, sample_count, 'NA')
        next
    }
        
    fit <- survfit(Surv(time, status) ~ variable, data = data_analysis)
    surv_diff <- survdiff(Surv(time, status) ~ variable, data = data_analysis) # 差异分析
    pvalue    <- pchisq(surv_diff$chisq, length(surv_diff$n)-1, lower.tail = FALSE) # 获取P值
    surv_info <- surv_summary(fit, data = data_analysis) # 拟合详细结果 

    # 记录pvalue
    pvalue_data[count, ] = c(feature, sample_count, pvalue)

    # 部分绘图参数预处理
    if(!is.null(pdf_legend_x) & !is.null(pdf_legend_y)) pdf_legend = c(as.numeric(pdf_legend_x), as.numeric(pdf_legend_y))
    if(!is.null(pdf_break_x_by)) pdf_break_x_by = as.numeric(pdf_break_x_by)

    # 生存曲线绘图
    p <- ggsurvplot(fit, data = data_analysis,
                    fun = pdf_fun,
                    pval = pdf_pvalue, 
                    conf.int = pdf_conf_int,
                    risk.table = pdf_risk_table, # Add risk table
                    # risk.table.col = "strata", # Change risk table color by groups
                    # linetype = "strata", # Change line type by groups
                    surv.median.line = pdf_surv_median_line, # Specify median survival
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    add.all = pdf_add_all,
                    palette = pdf_palette, #杂志nature的配色
                    legend.title= ifelse(is.null(pdf_legend_title), feature, pdf_legend_title),
                    legend = pdf_legend,
                    xlab = pdf_xlab, break.x.by = pdf_break_x_by,
                    title = pdf_title,
                    

                    )
        
        # 第一幅图要设置newpage = FALSE, 否则会多一个空白页
        result = tryCatch(print(p, newpage = count !=1),
                            error=function(e){cat("ERROR :",conditionMessage(e),"\n")}
                            ) # 注部分图可能无法绘制，报错
        
}
dev.off()

write.table(pvalue_data, output, row.names = F, quote = F, sep = '\t')

