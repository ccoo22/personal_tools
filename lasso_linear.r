#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)

"Usage: lasso.r -y <file> -x <file> -o <dir> -p <string> [--lambda <string> --lambda_value <numeric> --score_data <file> --rlib <dir>]

Options:
    -y <file>              因变量，两列数据，第一列样本名，第二列表型数据（连续型变量）；第一行为表头。
                           注：分析时，样本以该文件为准。自变量中的样本只能比该文件多。样本顺序没有要求。
    -x <file>              自变量，每一行为一个样本，每一列为一个特征，第一行为表头。第一列为样本名。
                           注：不允许缺失，如果有缺失，流程自动删除存在缺失数据的样本。
    -o, --output_dir <dir>          输出目录
    -p, --prefix <string>           输出文件前缀    
                                    输出的文件有：
                                        prefix.coef.txt 特征系数
                                        prefix.lambda.lambda.txt 当前系数对应的lambda值
                                        prefix.log_lambda_coefficients.pdf 系数随lambda变动的曲线
                                        prefix.log_lambda_mean-squared-error.pdf 分类错误率随lambda变动曲线
                                        prefix.sample_score_based_on_coef.txt 基于选定的lambda值得到的系数，计算每个样本的打分（>0是case组， <0是control组）
    --lambda <string>               lasso分析lambda值选择 min/1se [default: 1se]
                                    min：即选择错了率最低值对应的lambda
                                    1se: 指在min一个方差范围内得到最简单模型的那一个lambda值，1se给出的是一个具备优良性能且自变量个数最少的模型。
    --lambda_value <numeric>        指定lambda值。当指定lambda值时，--lambda参数会被忽略
    --score_data <file>             与input文件行的顺序完全一样的矩阵，对该矩阵中的样本使用lasso分析的系数进行打分
                                    该文件可以不提供                   
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]" -> doc

opts   <- docopt(doc, version='甘斌，lasso分析, 基于连续变量\n')
input_x             <- opts$x
input_y             <- opts$y
output_dir          <- opts$output_dir
prefix              <- opts$prefix

lambda              <- opts$lambda
lambda_value        <- opts$lambda_value
score_data          <- opts$score_data
rlib                <- opts$rlib
.libPaths(rlib)

# 导入 -> package
library(glmnet)
set.seed(91)

###################################################################### 主程序
message("read input")
tmp = read.table(input_x, head = T, row.names = 1, check.names = F, sep = "\t", colClass='character')
data_input_x = read.table(input_x, head = T, row.names = 1, check.names = F, sep = "\t", colClass=c('character', rep('numeric', ncol(tmp))))  # 第一列必须是字符型，防止R自动转换
data_input_y = read.table(input_y, head = T, row.names = 1, check.names = F, sep = "\t", colClass=c('character', 'numeric'))

# 错误检测
message("sample/group check")
if(sum(!rownames(data_input_y) %in% rownames(data_input_x)) > 0)
{   
    lost_samples = rownames(data_input_y)[!rownames(data_input_y) %in% rownames(data_input_x)]
    message("[Error]  y中部分样本在x中缺失，请务必保证x中包含y中的所有样本")
    message("[Error]  丢失的样本有： ",  paste(lost_samples, collapse=','))
    q()
}

# 数据准备
data_x = data_input_x[rownames(data_input_y), ]
non_miss_sample = complete.cases(data_x)

y_clean = as.matrix(data_input_y[non_miss_sample, , drop =F])
x_clean = as.matrix(data_x[non_miss_sample, ])


# 对特征重命名，防止有特殊符号，导致报错
feature_names_input         = colnames(x_clean)
feature_names_new           = paste0('gene', 1:ncol(x_clean))
colnames(x_clean)           = feature_names_new
names(feature_names_input)  = feature_names_new
 

# lasso 分析: 这里只做连续型因变量的模型。 gaussian
message("lasso fit")
cvfit  = cv.glmnet(x_clean, y_clean, type.measure = "mse", nfolds = 20)

# lasso 系数确认/输出
lasso_lambda_value = ifelse(lambda == 'min', cvfit$lambda.min, cvfit$lambda.1se)
lasso_lambda_value
if(!is.null(lambda_value)) 
{   
    message('in')
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

if(sum(lasso_coef[,1] != 0) == 0)
{
    message("[Warings]  当前lambda值下，没有筛选出任何变量, 请更换lambda值\n")
}
    
# 绘图 : 随着lambda的变化，筛选出来的变量组合造成的损失，损失越小越好（最左侧虚线），当然也经常用1se值（最右侧虚线）
pdf(paste0(output_dir, "/", prefix, '.log_lambda_mean-squared-error.pdf'))
plot(cvfit)
dev.off()

# 绘图 ：直接用原始数据拟合，查看随着lambda变化，每一个特征的系数变化
fit<-glmnet(x_clean, y_clean, family = "gaussian")
pdf(paste0(output_dir, "/", prefix, '.log_lambda_coefficients.pdf'))
plot(fit, xvar="lambda")
dev.off()

# 基于lambda值、x值，预测样本分类
# 注意： b = x_clean %*% lasso_coef[2:nrow(lasso_coef),] + lasso_coef[1,1] 
# pred = as.integer(predict(fit, newx = x_clean, type = "class", s = lasso_lambda_value)) # 

# 基于系数计算样本的得分
lasso_sample_score_file = paste0(output_dir, "/", prefix, '.sample_score_based_on_coef.txt')
score = x_clean %*% lasso_coef[2:nrow(lasso_coef),] + lasso_coef[1,1]  # 根据系数、样本数据，得到样本的打分，> 0是样本设为1，小于0时样本设为0， 这样得到的向量与pred结果完全一样
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
    score = data_input_score %*% lasso_coef[2:nrow(lasso_coef),] + lasso_coef[1,1]  #  
    score_final = data.frame(sample=rownames(score), score=score[,1])
    write.table(score_final, lasso_sample_score_file, quote = FALSE, row.names = FALSE, sep = '\t')


}