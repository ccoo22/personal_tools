#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: lasso_survival_cox.r --surv <file> --exp <file> -o <dir> --prefix <string>  [--gene_file <file>  --lambda <string> --lambda_value <numeric> --score_data <file> --rlib <dir>]

Options:
    -s, --surv <file>               lasso生存信息矩阵，第一列样本名，第二列time, 第三列 status（0 alive, 1 dead）， 有表头
    --exp <file>                    lasso分析表达量文件矩阵，每一行对应一个特征，每一列对应一个样本，第一行是样本名，第一列是特征名称。不允许缺失。如果有缺失，会自动排除包含缺失的样本。
                                    注意：流程只会分析--surv文件中包含的样本
    --gene_file <file>              指定分析的基因，一列数据，一行一个基因,没有表头。如果不指定，默认用 exp 的所有基因
    -o, --output_dir <dir>          输出目录
    -p, --prefix <string>           输出文件前缀    
                                    输出的文件有：
                                        prefix.coef.txt 特征系数
                                        prefix.lambda.lambda.txt 当前系数对应的lambda值
                                        prefix.log_lambda_coefficients.pdf 系数随lambda变动的曲线
                                        prefix.log_lambda_partial_likelihood_deviance.pdf 异常值随lambda变动曲线
                                        prefix.sample_score_based_on_coef.txt 基于选定的lambda值得到的系数，计算每个样本的打分（>0是case组， <0是control组）
    --lambda <string>               lasso分析lambda值选择 min/1se [default: 1se]
                                    min：即选择错了率最低值对应的lambda
                                    1se: 指在min一个方差范围内得到最简单模型的那一个lambda值，1se给出的是一个具备优良性能且自变量个数最少的模型。
    --lambda_value <numeric>        指定lambda值。当指定lambda值时，--lambda参数会被忽略
    --score_data <file>             与input文件行的顺序完全一样的矩阵，对该矩阵中的样本使用lasso分析的系数进行打分
                                    该文件可以不提供                   
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，lasso分析， 基于生存信息\n')
surv                <- opts$surv
exp                 <- opts$exp
output_dir          <- opts$output_dir
prefix              <- opts$prefix
gene_file           <- opts$gene_file

lambda              <- opts$lambda
lambda_value        <- opts$lambda_value
score_data          <- opts$score_data
rlib                <- opts$rlib
.libPaths(rlib)
 
# 导入 -> package
library(glmnet)
# library(ROCR)
# library(pROC)
set.seed(91)

###################################################################### 主程序
message("read input")
data_surv  = read.table(surv, head = T, row.names = 1, check.names = F, sep = "\t")  # 生存信息
data_exp   = read.table(exp, head = T, row.names = 1, check.names = F, sep = "\t")  # 表型信息
data_exp   = t(data_exp)
samples = rownames(data_surv)
colnames(data_surv) = c('time', 'status')  # 必须这么命名，lasso要用到

if(sum(!samples %in% rownames(data_exp)) > 0)
{
    losts = samples[!samples %in% rownames(data_exp)]
    message("[Error] surv 中的样本在 exp 文件中缺失:", losts)
    q()
}

# 基因选择
if(! is.null(gene_file))
{
    gene_file_data <- read.table(gene_file, sep='\t', header = FALSE, check.names=FALSE, stringsAsFactors = F) 
    genes = gene_file_data[,1]
    lost_genes = genes[!genes %in% colnames(data_exp)]
    if(length(lost_genes) > 0)
    {
        message("指定的基因不存在：", paste(lost_genes, collapse=','))
        q()
    }
    data_exp = data_exp[, genes , drop=F] 
    message("指定分析 ", length(genes), " 个基因")
}


# lasso分析准备
message("prepare lasso input")

x_select          = data_exp[samples, ]
y_select          = data_surv

# 对特征重命名，防止有特殊符号，导致报错
feature_names_input         = colnames(x_select)
feature_names_new           = paste0('gene', 1:ncol(x_select))
colnames(x_select)          = feature_names_new
names(feature_names_input)  = feature_names_new

# 去除缺失
x_clean = data.matrix(x_select)
y_clean = data.matrix(y_select)
non_miss_sample = complete.cases(x_select)
if(sum(!non_miss_sample) > 0)
{   
    lost = sum(!non_miss_sample)
    message("存在缺失数据: ", lost, " 个");
    x_clean =  data.matrix(x_select[non_miss_sample, ])
    y_clean =  data.matrix(y_select[non_miss_sample, ] )
}


# lasso 分析: 这里只做二分类的模型。binomial
# lasso支持单个连续因变量y "gaussian"
# lasso支持多个连续因变量y mgaussian
# lasso支持生存分析cox
message("lasso fit")
cvfit  = cv.glmnet(x_clean, y_clean, family = "cox")

# lasso 系数确认/输出
lasso_lambda_value = ifelse(lambda == 'min', cvfit$lambda.min, cvfit$lambda.1se)
if(length(lambda_value) > 0) 
{
    lasso_lambda_value = as.numeric(lambda_value)  # 使用指定lambda值
    lambda = 'manual_defined'
}
 
lasso_coef         = data.matrix(coef(cvfit, s=lasso_lambda_value))  # 系数为0或者.，表示在当前lambda值下，该系数被消掉。注意：该系数只是lasso的系数，不是逻辑回归的系数。
# 特征名称恢复
lasso_coef_rename  = data.frame(feature = feature_names_input[rownames(lasso_coef)], coef = lasso_coef[,1], row.names = NULL, stringsAsFactors=F)
lasso_coef_rename$feature[which(is.na(lasso_coef_rename$feature))] = '(Intercept)'
# 输出到文件
lasso_coef_file = paste0(output_dir, "/", prefix, '.coef.txt')
lasso_lambda_file = paste0(output_dir, "/", prefix, '.lambda.', lambda, ".txt")
write.table(lasso_lambda_value, lasso_lambda_file, quote = FALSE, row.names = FALSE, col.names=F, sep = '\t')
write.table(lasso_coef_rename, lasso_coef_file, quote = FALSE, row.names = FALSE, sep = '\t')

  
 # 绘图 : 随着lambda的变化，筛选出来的变量组合造成的损失，损失越小越好（最左侧虚线），当然也经常用1se值（最右侧虚线）
pdf(paste0(output_dir, "/", prefix, '.log_lambda_partial_likelihood_deviance.pdf'))
plot(cvfit)
dev.off()

# 绘图 ：直接用原始数据拟合，查看随着lambda变化，每一个特征的系数变化
fit<-glmnet(x_clean, y_clean, family = "cox")
pdf(paste0(output_dir, "/", prefix, '.log_lambda_coefficients.pdf'))
plot(fit, xvar="lambda")
dev.off()

if(sum(lasso_coef[,1] != 0) == 0)
{
    message("[Warings]  当前lambda值下，没有筛选出任何变量, 请更换lambda值\n")
    q()
}
# 基于lambda值、x值，预测样本分类
# 注意： b = x_clean %*% lasso_coef[2:nrow(lasso_coef),] + lasso_coef[1,1] 
# pred = as.integer(predict(fit, newx = x_clean, type = "class", s = lasso_lambda_value)) # 

# 基于系数计算样本的得分
lasso_sample_score_file = paste0(output_dir, "/", prefix, '.sample_score_based_on_coef.txt')
score = data_exp %*% lasso_coef  # 根据系数、样本数据，得到样本的打分，> 0是样本设为1，小于0时样本设为0， 这样得到的向量与pred结果完全一样
score_final = data.frame(sample=rownames(score), score=score[,1])
write.table(score_final, lasso_sample_score_file, quote = FALSE, row.names = FALSE, sep = '\t')


if(!is.null(score_data))
{
    message('使用input数据得到的系数，对 score_data 中的样本打分')
    data_input_score = t(read.table(score_data, head = T, row.names = 1, check.names = F, sep = "\t"))
    if(sum(!feature_names_input %in% colnames(data_input_score)) > 0 )
    {
        message('[Error] score_data文件的特征与input不符，无法完成打分任务')
        q()
    }

    # 排序
    data_input_score = data_input_score[,feature_names_input]

    # 打分输出
    lasso_sample_score_file = paste0(output_dir, "/", prefix, '.score_data.txt')
    score = data_input_score %*% lasso_coef  #  
    score_final = data.frame(sample=rownames(score), score=score[,1])
    write.table(score_final, lasso_sample_score_file, quote = FALSE, row.names = FALSE, sep = '\t')


}
