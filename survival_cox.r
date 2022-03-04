#!/home/genesky/software/r/4.0.3/bin/Rscript
library(docopt)
"Usage: survival_cox.r --surv <file> --exp <file> -o <file> -p <pdf> [--pvalue_cutoff <numeric>  --gene_file <file> --exp_format <string> -m <string> --pdf_width <numeric> --pdf_height <numeric> --sort_by_hazard_ratio]
Options:
    -s, --surv <file>               生存信息矩阵，第一列样本名，第二列time, 第三列 status（0 alive, 1 dead）， 有表头
    --exp <file>                    分析表达量文件矩阵，每一行对应一个特征，每一列对应一个样本，第一行是样本名，第一列是特征名称。不允许缺失。如果有缺失，会自动排除包含缺失的样本。
                                    注意1：流程只会分析--surv文件中包含的样本
                                    注意2：另外，输入的文件也可以是每一行对应一个样本，每一列对应一个特征。数据格式请通过exp_format参数进行声明
    --exp_format <string>           exp文件格式 col/row。col表示exp文件每一列对应一个样本；row表示exp文件每一行对应一个样本 [default: col]
    -m, --model <string>            分析模式，有两种:univariate / multivariate,，一种是所有exp特种放在一起，做多因素生存分析；另一种是每一个因素单独分析，做单因素生存分析 [default: univariate]
    -o, --output <file>             统计结果文件
    -p, --pdf_file <pdf>            统计结果绘图文件pdf
    --pvalue_cutoff <numeric>       输出pvalue 小于该值的基因文件 [default: 0.05]
    --gene_file <file>              指定分析的基因，一列数据，一行一个基因,没有表头。如果不指定，默认用 exp 的所有基因
    --pdf_width <numeric>           PDF宽度 [default: 7]
    --pdf_height <numeric>          PDF高度 [default: 7]
    --sort_by_hazard_ratio          PDF绘图时，按照hazard ratio值排序。
" -> doc

opts               <- docopt(doc, version = 'COX生存分析软件 \n')
surv              <- opts$surv
exp               <- opts$exp
exp_format        <- opts$exp_format
output            <- opts$output
pdf_file          <- opts$pdf_file
model             <- opts$model
gene_file         <- opts$gene_file
pdf_width         <- as.numeric(opts$pdf_width)
pdf_height        <- as.numeric(opts$pdf_height)
pvalue_cutoff     <- as.numeric(opts$pvalue_cutoff)
sort_by_hazard_ratio        <- opts$sort_by_hazard_ratio




library(survival)
library(forestplot)

# 读入数据
message("read surv data")
surv_data = read.table(surv, head = T, row.names = 1, check.names = F, sep = '\t', stringsAsFactors = F)
surv_data = surv_data[, c(1,2)]
colnames(surv_data)[1:2] = c('time', 'status') 
samples = rownames(surv_data)

message("read exp data")
data_exp   = read.table(exp, head = T, row.names = 1, check.names = F, sep = "\t", stringsAsFactors = F)  # 表型信息
if(exp_format == 'col') data_exp   = t(data_exp)

if(sum(!samples %in% rownames(data_exp)) > 0)
{
    losts = samples[!samples %in% rownames(data_exp)]
    message("[Error] surv 中的样本在 exp 文件中缺失:", losts)
    q()
}
data_exp_clean = data_exp[samples, , drop=F] 

# 基因选择
if(! is.null(gene_file))
{
    gene_file_data <- read.table(gene_file, sep='\t', header = FALSE, check.names=FALSE, stringsAsFactors = F) 
    genes = gene_file_data[,1]
    lost_genes = genes[!genes %in% colnames(data_exp_clean)]
    if(length(lost_genes) > 0)
    {
        message("指定的基因不存在：", paste(lost_genes, collapse=','))
        q()
    }
    data_exp_clean = data_exp_clean[, genes , drop=F] 
    message("指定分析 ", length(genes), " 个基因")
}


# 特征名称临时替换，防止存在空格、-等特殊字符，导致报错
feature_names_input         = colnames(data_exp_clean)
feature_names_new           = paste0('gene', 1:ncol(data_exp_clean))
colnames(data_exp_clean)    = feature_names_new
names(feature_names_input)  = feature_names_new

##################
# COX分析
##################
# exp(coef) 就是Hazrd Ratio. HR = 1 无影响。HR < 1, 降低风险，保护因素。HR > 1, 增加风险，危险因素。
message("start COX analysis")
result = matrix(nrow=length(feature_names_new), ncol=6)
if(model == 'univariate')
{   
    row = 0
    for(feature in feature_names_new)
    {   
        row = row + 1
        # 准备输入数据、去掉缺失
        data_analysis = cbind(surv_data, data_exp_clean[, feature])
        colnames(data_analysis)[3] = feature
        data_analysis <- data_analysis[complete.cases(data_analysis), ]
        data_analysis <- data_analysis[apply(data_analysis, 1, function(x){sum(x=='')}) == 0, ]
       
        # 拟合
        formula <- as.formula( paste0('Surv(time, status) ~ ',  feature ))
        fit <- coxph(formula, data = data_analysis)
        res_coef <- summary(fit)$coefficients[1, c(1, 2, 5)]
        res_confint <- summary(fit)$conf.int[1, c(3,4)]
        result[row, ] = c(feature_names_input[feature], res_coef, res_confint)

    }
}else
{
    # 准备输入数据、去掉缺失
    data_analysis = cbind(surv_data, data_exp_clean)
    colnames(data_analysis)[3:ncol(data_analysis)] = feature_names_new
    data_analysis <- data_analysis[complete.cases(data_analysis), ]
    data_analysis <- data_analysis[apply(data_analysis, 1, function(x){sum(x=='')}) == 0, ]
       
    # 拟合
    formula <- as.formula( paste0('Surv(time, status) ~ ', paste(feature_names_new, collapse = ' + ') ))
    fit <- coxph(formula, data = data_analysis)

    res_coef <- summary(fit)$coefficients[, c(1, 2, 5)]
    res_confint <- summary(fit)$conf.int[, c(3,4)]
    result <- cbind(feature = feature_names_input[rownames(res_coef)], res_coef, res_confint)
}

colnames(result) = c('feature', 'coef', 'exp(coef)', 'pvalue', 'lower .95', 'upper .95')
result = data.frame(result, check.names=F)
result[,'pvalue'] = as.numeric(result[,'pvalue'])
write.table(result, output, row.names = F, quote = F, sep = '\t')

sig_gene = result[,'feature'][result[, 'pvalue'] < pvalue_cutoff]
message('共 ', length(sig_gene), ' 个显著差异基因')
sig_gene = data.frame(gene=sig_gene, stringsAsFactors = F)
write.table(sig_gene, paste0(output, ".sig_gene.txt"), col.names = F, row.names=F, sep = "\t", quote=F)

# 绘图
message("start plot")
if(sort_by_hazard_ratio) result = result[order(as.numeric(result[, 'exp(coef)'])), ]


feature    = result[, 'feature']
pvalue     = format(as.numeric(result[, 'pvalue']), scientific = TRUE, digits=4)
hr         = round(as.numeric(result[, 'exp(coef)']), 3)
hr_low     = round(as.numeric(result[, 'lower .95']), 3)
hr_up      = round(as.numeric(result[, 'upper .95']), 3)

data_plot = cbind(c(NA, feature), 
                    c('pvalue', pvalue),
                    c('Hazard ratio', paste0(hr, '(', hr_low, '-', hr_up, ')'))
                 )
pdf(pdf_file, width = pdf_width, height = pdf_height)
forestplot(
    labeltext = data_plot,
    mean = c(NA, hr),
    lower = c(NA, hr_low),
    upper = c(NA, hr_up),
    zero = 1,  # 显示y=n的垂直线
    xlab = 'Hazard ratio',  # x轴的标题
    boxsize = 0.3, # 误差条中的正方形大小
    lty.ci = 7, # 误差条的线的线型
    lwd.ci = 3,# 误差条的线的宽度
    col=fpColors(line = "#000000", #误差条的线的颜色
                        box="#0000CD"), #误差条的正方形的颜色
    new_page=FALSE,
)
dev.off()
