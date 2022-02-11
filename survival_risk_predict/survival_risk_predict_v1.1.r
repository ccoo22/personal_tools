#!/home/genesky/software/r/3.5.1/bin/Rscript

library(docopt)
"Usage: survival_analysis.r --surv <file> --exp <file> -o <dir> [--exp_format <string> --do_log2 --gene <file> --gene_anno <file> --method <string> --time_cutoff <string> --time_label <string> --rlib <dir>]
Options:
    -s, --surv <file>               生存信息矩阵，第一列样本名，第二列time, 第三列 status（0 alive, 1 dead），有表头。不允许缺失。如果有缺失，会自动排除包含缺失的样本。
                                    示例：/home/zhangshuang/work/research/survival_hazard_analysis_forestplot/script/v2.0/example/input/sample.survival.txt
    --exp <file>                    分析表达量矩阵，每一行对应一个特征，每一列对应一个样本，第一行是样本名，第一列是特征名称。不允许缺失。如果有缺失，会自动排除包含缺失的样本。
                                    注意1：流程只会分析--surv文件中包含的样本
                                    注意2：另外，输入的文件也可以是每一行对应一个样本，每一列对应一个特征。数据格式请通过exp_format参数进行声明
    -o, --output <dir>              结果输出目录
    --exp_format <string>           exp文件格式 col/row。col表示exp文件每一列对应一个样本；row表示exp文件每一行对应一个样本 [default: col]
    --do_log2                       绘制热图时，对表达量数据做log2处理
    --gene <file>                   基因列表文件，一列数据，无表头，指定生存分析的基因，不提供则默认分析表达量矩阵中的所有基因
    --gene_anno <file>              基因注释文件，两列数据，无表头，第一列为表达量矩阵中基因ID，第二列为相应基因在结果中展示的基因名称，不提供则使用表达量矩阵中基因名称
    --method <string>               ROC fitting方法: NNE 或 KM，默认 KM [default: KM]
    --time_cutoff <string>          suvival ROC 与 suvival logrank 研究、统计的截止时间，年月日等。逗号分隔:\'365,1095,1825,2920\', \'6,12\' [default: 365,1095,1825,2920]
    --time_label <string>           对应截止时间的标签名，逗号分隔: \'1 year,3 year,5 year,8 year\', \'6 months,12 months\' [default: 1 year,3 year,5 year,8 year]
    --rlib <dir>                    R包路径 [default: /home/genesky/software/r/3.5.1/lib64/R/library]
" -> doc

opts              <- docopt(doc, version = 'suvival risk predict v1.1 \n zhangshuang 298\n')
surv              <- opts$surv
exp               <- opts$exp
output            <- opts$output
exp_format        <- opts$exp_format
do_log2           <- opts$do_log2
gene              <- opts$gene
gene_anno         <- opts$gene_anno
method            <- opts$method
time_cutoff       <- opts$time_cutoff
time_label        <- opts$time_label
rlib              <- opts$rlib


# surv              <- "/home/zhangshuang/work/other_project/20B1217C_survival_analysis/data/Overall_Survival_Time.txt"
# exp               <- "/home/zhangshuang/work/other_project/20B1217C_survival_analysis/data/exp.txt"
# output            <- "/home/zhangshuang/work/other_project/20B1217C_survival_analysis/test"
# exp_format        <- "col"
# do_log2           <- TRUE
# gene              <- "/home/zhangshuang/work/other_project/20B1217C_survival_analysis/data/gene.txt"
# gene_anno         <- "/home/zhangshuang/work/other_project/20B1217C_survival_analysis/data/gene_anno.txt"
# method            <- "KM"
# time_cutoff       <- '365,1095,1825,2920'
# time_label        <- '1 year,3 year,5 year,8 year'
# rlib              <- "/home/genesky/software/r/3.5.1/lib64/R/library"

############
# 加载R包
############
.libPaths(rlib)
library(survival)
library(survminer)
library(survivalROC)
library(patchwork)
library(pheatmap)
library(ggplot2)
library(ggpubr)

############
# 读入数据并去除缺失
############
message("read surv data")
data_surv = read.table(surv, header = T, row.names = 1, stringsAsFactors = FALSE, check.names = F, sep = '\t')  # 生存信息
data_surv = data_surv[, c(1,2)]
colnames(data_surv)[1:2] = c('time', 'status')
samples = rownames(data_surv)

message("read exp data")
data_exp = read.table(exp, header = T, row.names = 1, stringsAsFactors = FALSE, check.names = F, sep = "\t")
if(exp_format == 'col') data_exp   = t(data_exp)

# 指定分析基因
if(!is.null(gene))
{
    data_gene <- read.table(gene, header = F, stringsAsFactors = FALSE, check.names = F, sep = "\t")
    rownames(data_gene) <- data_gene[, 1]
    colnames(data_gene) <- "Gene"
    gene_analysis <- data_gene[, 1]
    if(sum(!gene_analysis %in% colnames(data_exp)) > 0)
    {
        lost_genes = gene_analysis[!gene_analysis %in% colnames(data_exp)]
        message("[Error] 指定分析的基因在 exp 文件中缺失:", lost_genes)
        q()
    }
    data_exp <- data_exp[, gene_analysis, drop = F]
}else
{
    gene_analysis <- colnames(data_exp)
}

if(sum(!samples %in% rownames(data_exp)) > 0)
{
    lost_samples = samples[!samples %in% rownames(data_exp)]
    message("[Error] surv 中的样本在 exp 文件中缺失:", lost_samples)
    q()
}
data_exp_clean = data_exp[samples, , drop = F] 

# 根据提供基因注释信息，转换基因名
if(!is.null(gene_anno))
{
    message("read gene anno data")
    data_gene_anno <- read.table(gene_anno, header = F, stringsAsFactors = FALSE, check.names = F, sep = "\t")
    rownames(data_gene_anno) <- data_gene_anno[, 1]
    if(sum(!colnames(data_exp_clean) %in% data_gene_anno[, 1]) > 0)
    {
        unmatch_genes = colnames(data_exp_clean)[!colnames(data_exp_clean) %in% data_gene_anno[, 1]]
        message("[Error] exp 中的基因在 gene_anno 文件中缺失:", unmatch_genes)
        q()
    }
    colnames(data_exp_clean) <- data_gene_anno[colnames(data_exp_clean), 2]
}

# 特征名称临时替换，防止存在空格、-等特殊字符，导致报错
feature_names_input         = colnames(data_exp_clean)
feature_names_new           = paste0('gene', 1:ncol(data_exp_clean))
colnames(data_exp_clean)    = feature_names_new
names(feature_names_input)  = feature_names_new

##################
# COX分析（生存状态与feature）
##################
message("start COX analysis(Survive Status with Gene)")

# 准备输入数据
data_analysis_surv = cbind(data_surv, data_exp_clean, stringsAsFactors = FALSE)
data_analysis_surv <- data_analysis_surv[complete.cases(data_analysis_surv), ]
data_analysis_surv <- data_analysis_surv[apply(data_analysis_surv, 1, function(x){sum(x == '')}) == 0, ]
if(sum(data_analysis_surv$status == 0) == 0 || sum(data_analysis_surv$status == 1) == 0)
{
    message("[Error] 过滤缺失后，剩余样本生存状态 status 为单一类别，无法继续分析")
    q()
}

# 拟合
formula_surv <- as.formula( paste0('Surv(time, status) ~ ', paste(feature_names_new, collapse = ' + ') ))
fit_surv <- coxph(formula_surv, data = data_analysis_surv)

coef_surv <- summary(fit_surv)$coefficients[, c(1, 2, 5)]
confint_surv <- summary(fit_surv)$conf.int[, c(3,4)]

# 保存 Surv feature COX 分析结果
result_surv <- cbind(feature = feature_names_input[rownames(coef_surv)], coef_surv, confint_surv)
colnames(result_surv) <- c('Feature', 'coef', 'exp(coef)', 'pvalue', 'lower .95', 'upper .95')
surv_gene_file <- paste(output, "/survival.cox.gene.txt", sep = "")
write.table(result_surv, surv_gene_file, row.names = F, quote = F, sep = '\t')

##################
# 计算风险打分、基于风险打分分组
##################
message("start risk score calculate")

risk_score <- as.matrix(data_analysis_surv[, rownames(coef_surv), drop = F]) %*% coef_surv[, 1]
data_analysis_surv$risk_score <- risk_score
data_analysis_surv$risk_group <- "Low"
data_analysis_surv$risk_group[data_analysis_surv$risk_score >= median(risk_score)] = "High"
samples_sort_by_risk_score <- rownames(data_analysis_surv)[order(data_analysis_surv$risk_score)]

# 保存分组结果
result_risk_group <- data.frame(rownames(data_analysis_surv), data_analysis_surv[, c("time", "status", "risk_score", "risk_group")])
colnames(result_risk_group) <- c("Sample", "Time", "Status", "Risk Score", "Risk Group")
risk_group_file <- paste(output, "/sample.risk.result.txt", sep = "")
write.table(result_risk_group, risk_group_file, row.names = F, quote = F, sep = '\t')

##################
# 热图、risk score图
##################
data_for_heatmap_plot <- t(data_exp_clean[samples_sort_by_risk_score, ])
if(!is.null(do_log2)) data_for_heatmap_plot <- log2(data_for_heatmap_plot + 1)
rownames(data_for_heatmap_plot) <- feature_names_input[rownames(data_for_heatmap_plot)]

data_for_risk_score_dot_plot <- data.frame("sample" = c(1:length(samples_sort_by_risk_score)), "risk_score" = data_analysis_surv[samples_sort_by_risk_score, "risk_score"], "risk_group" = data_analysis_surv[samples_sort_by_risk_score, "risk_group"])
rownames(data_for_risk_score_dot_plot) <- samples_sort_by_risk_score

data_for_survive_status_plot <- data.frame("sample" = c(1:length(samples_sort_by_risk_score)), "time" = data_analysis_surv[samples_sort_by_risk_score, 'time'], "status" = data_analysis_surv[samples_sort_by_risk_score, 'status'])
data_for_survive_status_plot$status <- factor(data_for_survive_status_plot$status, levels = 1:0, labels = c("Dead", "Alive"))
rownames(data_for_survive_status_plot) <- samples_sort_by_risk_score

myheatcol = colorRampPalette(c("#CC3366", "white", "#00FFCC"))(100)
# 基因数大于50，字号调整
fontsize_col <- 10
if(ncol(data_exp_clean) > 50){fontsize_col <- 500/ncol(data_exp_clean)}

p1 = pheatmap(data_for_heatmap_plot,
    scale = 'row',
    margins = c(8,10),
    cluster_rows = T,
    cluster_cols = F,
    color = myheatcol,
    fontsize_col = fontsize_col,
    show_rownames = T,
    show_colnames = F
)

p2 = ggplot(data_for_risk_score_dot_plot) +
    geom_point(aes(x = sample, y = risk_score, color = risk_group), size = 3) +
    scale_color_manual(values = c('#CC3366', '#00FFCC')) +
    geom_vline(xintercept = length(samples_sort_by_risk_score)/2, lty = 3, col = "red", lwd = 1) +
    labs(x = '', y = 'risk score', title = '') +
    theme(panel.background = element_rect(fill = 'transparent', color = 'black'),
        panel.border = element_rect(fill = 'transparent', color = 'transparent'),
        panel.grid = element_blank(),
        axis.title.y = element_text(size = 10, face = 'bold'),
        legend.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 10, face = 'bold'),
        legend.position = 'right',
        legend.key.size = unit(0.8,'cm'))

p3 = ggplot(data_for_survive_status_plot) +
    geom_point(aes(x = sample, y = time, color = status), size = 3) +
    scale_color_manual(values = c("Dead" = '#CC3366', "Alive" = '#00FFCC'), name = "Vital_Status") +
    geom_vline(xintercept = length(samples_sort_by_risk_score)/2, lty = 3, col = "red", lwd = 1) +
    labs(x = '', y = 'Days', title = '') +
    theme(panel.background = element_rect(fill = 'transparent', color = 'black'),
        panel.border = element_rect(fill = 'transparent', color = 'transparent'),
        panel.grid = element_blank(),
        axis.title.y = element_text(size = 10, face = 'bold'),
        legend.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 10, face = 'bold'),
        legend.position = 'right',
        legend.key.size = unit(0.8,'cm'))

risk_status_heatmap_file <- paste(output, "/sample.risk.status.heatmap.pdf", sep = "")
pdf(risk_status_heatmap_file, width = 14)
print(p1)
# patchwork 等拼接工具只适合gg系列R包绘制的图
print(p2/p3)
dev.off()


##################
# Survive Logrank 绘图
##################
message("start survival logrank plot")

pdf_survival_logrank <- paste(output, "/survival.logrank.pdf", sep = "")
data_analysis_survival_logrank <- data_analysis_surv[, c("time", "status", "risk_group")]

fit_survival_logrank       <- survfit(Surv(time, status) ~ risk_group, data = data_analysis_survival_logrank)
surv_diff_survival_logrank <- survdiff(Surv(time, status) ~ risk_group, data = data_analysis_survival_logrank) # 差异分析
pvalue_survival_logrank    <- pchisq(surv_diff_survival_logrank$chisq, length(surv_diff_survival_logrank$n)-1, lower.tail = FALSE) # 获取P值
surv_info_survival_logrank <- surv_summary(fit_survival_logrank, data = data_analysis_survival_logrank) # 拟合详细结果 

# 保存 Survive Logrank 拟合结果
result_surv_logrank_info_low <- surv_info_survival_logrank[surv_info_survival_logrank$risk == 'Low', ]
result_surv_logrank_info_high <- surv_info_survival_logrank[surv_info_survival_logrank$risk == 'High', ]
surv_logrank_info_group_low_file <- paste(output, "/Low.risk.group.surv.logrank.info.txt", sep = "")
surv_logrank_info_group_high_file <- paste(output, "/High.risk.group.surv.logrank.info.txt", sep = "")
write.table(result_surv_logrank_info_low, surv_logrank_info_group_low_file, row.names = F, quote = F, sep = '\t')
write.table(result_surv_logrank_info_high, surv_logrank_info_group_high_file, row.names = F, quote = F, sep = '\t')

# 生存曲线绘图
pdf(pdf_survival_logrank, width = 7, height = 7)
p <- ggsurvplot(fit_survival_logrank, data = data_analysis_survival_logrank,
    fun = NULL,
    pval = TRUE, 
    conf.int = FALSE,
    risk.table = TRUE, # Add risk table
    surv.median.line = "hv", # Specify median survival
    ggtheme = theme_bw(), # Change ggplot2 theme
    add.all = FALSE,
    palette = "npg", #杂志nature的配色
    legend.title= "risk",
    legend = "top",
    xlab = "Time",
)
# 第一幅图要设置newpage = FALSE, 否则会多一个空白页
result = tryCatch(print(p, newpage = FALSE), error=function(e){cat("ERROR :",conditionMessage(e),"\n")}) # 注部分图可能无法绘制，报错
dev.off()


##################
# Survive ROC 分析绘图
##################
message("start survival ROC plot")

# 研究截止时间(默认 1 year 3 year 5 year 8 year)
cutoff_times <- unlist(strsplit(time_cutoff, ','))
time_labs <- unlist(strsplit(time_label, ','))

if(sum(is.na(cutoff_times)) > 0)
{
    message("[ERROR] 研究截止时间time_cutoff必须是数字形式")
    q()
}

if(length(cutoff_times) != length(time_labs))
{
    message("指定time_labs数量与time_cutoff数不一致")
    q()
}

# 颜色设置
mycol_ROC <- c(119,132,147,454,89,404,123,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)  
mycol_ROC <- colors()[mycol_ROC]
mycol_ROC <- mycol_ROC[1:length(cutoff_times)]
names(mycol_ROC) = cutoff_times

pdf_survival_roc <- paste(output, "/survival.roc.pdf", sep = "")
pdf(pdf_survival_roc)

for (i in 0:nrow(coef_surv))
{
    data_for_survivalROC_plot <- data.frame(sample = rownames(data_analysis_surv), time = data_analysis_surv$time, status = data_analysis_surv$status)
    if(i == 0)
    {
        # risk group ROC
        data_for_survivalROC_plot$marker_risk_value = data_analysis_surv$risk_group
        data_for_survivalROC_plot$marker_risk_value[data_for_survivalROC_plot$marker_risk_value == "High"] = 1
        data_for_survivalROC_plot$marker_risk_value[data_for_survivalROC_plot$marker_risk_value == "Low"] = 0
    }else
    {
        # 各基因 ROC
        data_for_survivalROC_plot$marker_risk_value = as.numeric(data_analysis_surv[, rownames(coef_surv)[i]])
    }
    # 记录 AUC
    val = c()
    for (j in 1:length(cutoff_times))
    {
        ### 选择fitting方法
        if (method == "KM")
        {
            ROC_opt= survivalROC(Stime = data_for_survivalROC_plot$time,
                    status = data_for_survivalROC_plot$status,
                    marker = data_for_survivalROC_plot$marker_risk_value,
                    predict.time = as.numeric(cutoff_times[j]), method = "KM")
        }else{

            nobs <- NROW(data_for_survivalROC_plot)
            ROC_opt= survivalROC(Stime = data_for_survivalROC_plot$time,
                    status = data_for_survivalROC_plot$status,
                    marker = data_for_survivalROC_plot$marker_risk_value,
                    predict.time = as.numeric(cutoff_times[j]), span = 0.25*nobs^(-0.20))
        }

        val[j] = round(ROC_opt$AUC, 3)

        if (j == 1)
        {
            plot(ROC_opt$FP, ROC_opt$TP, type="l", xlim=c(0, 1), ylim=c(0, 1),
                xlab="FP",
                ylab="TP",col = mycol_ROC[cutoff_times[j]])
            # 加对角线
            abline(0,1)
        }else{
            # 追加不同时间点的 ROC
            lines(ROC_opt$FP, ROC_opt$TP, type="l", xlim=c(0, 1), ylim=c(0, 1), col = mycol_ROC[cutoff_times[j]])
        }
    }
    # 记录图例文字
    label = c()
    for (k in 1:length(time_labs)){
        label[k] = paste("AUC at", time_labs[k], ":", val[k], sep = " ")
    }

    title(main = ifelse(i == 0, "Risk Group", feature_names_input[rownames(coef_surv)[i]]))
    legend("bottomright", label, lty=rep(1, length(time_labs)), col = mycol_ROC)
}
dev.off()


##################
# COX分析（生存状态与风险分组）
##################
message("start COX analysis(Survive Status with Risk Group)")

# 准备输入数据、去掉缺失
data_analysis_surv2 = data_analysis_surv[, c("time", "status", "risk_group")]
data_analysis_surv2$risk_group[data_analysis_surv2$risk_group == "High"] = 1
data_analysis_surv2$risk_group[data_analysis_surv2$risk_group == "Low"] = 0

# 拟合
formula_surv2 <- as.formula('Surv(time, status) ~ risk_group')
fit_surv2 <- coxph(formula_surv2, data = data_analysis_surv2)

coef_surv2 <- summary(fit_surv2)$coefficients[, c(1, 2, 5)]
confint_surv2 <- summary(fit_surv2)$conf.int[, c(3,4)]

# 保存分组生存分析统计结果
risk_groups <- c('Low', 'High')
os_rates_table <- matrix(NA, nrow = 2, ncol = length(cutoff_times)+1)
rownames(os_rates_table) <- risk_groups
colnames(os_rates_table) <- c('median', time_labs)

for (group in risk_groups)
{
    result_surv_logrank_info_group <- surv_info_survival_logrank[surv_info_survival_logrank$risk == group, ]
	for (i in 1:nrow(result_surv_logrank_info_group))
	{
		time <- as.numeric(result_surv_logrank_info_group[i, "time"])
		surv_rate <- as.numeric(result_surv_logrank_info_group[i, "surv"])
		if(surv_rate <= 0.5 && is.na(os_rates_table[group, 1])) os_rates_table[group, 1] = time
		for (j in 1:length(cutoff_times))
		{
			if(time >= cutoff_times[j] && is.na(os_rates_table[group, time_labs[j]])) os_rates_table[group, time_labs[j]] = paste(round(as.numeric(surv_rate*100), 1), "%", sep = "")
		}
	}
}

result_surv2 <- matrix(nrow = 2, ncol = 9)
rownames(result_surv2) <- risk_groups
colnames(result_surv2) <- c('Models and risk groups', 'Patient, n', 'Median OS Time', paste(time_labs, 'OS rates', sep = " "), 'HR(95%CI)', 'P value')

low_risk_sample <- rownames(data_analysis_surv)[data_analysis_surv$risk_group == "Low"]
high_risk_sample <- rownames(data_analysis_surv)[data_analysis_surv$risk_group == "High"]
low_risk_sample_num <- length(low_risk_sample)
high_risk_sample_num <- length(high_risk_sample)
low_risk_sample_percent <- paste(round(as.numeric(low_risk_sample_num*100/nrow(data_analysis_surv)), 0), "%", sep = "")
high_risk_sample_percent <- paste(round(as.numeric(high_risk_sample_num*100/nrow(data_analysis_surv)), 0), "%", sep = "")

hr = round(as.numeric(coef_surv2[2]), 3)
hr_low = round(as.numeric(confint_surv2[1]), 3)
hr_up = round(as.numeric(confint_surv2[2]), 3)
hazard_ratio <- paste0(hr, '(', hr_low, '-', hr_up, ')')
pvalue <- as.numeric(coef_surv2[3])

result_surv2["Low", ] = c("Low", paste0(low_risk_sample_num, '(', low_risk_sample_percent, ')'), os_rates_table["Low", ], "-", "")
result_surv2["High", ] = c("High", paste0(high_risk_sample_num, '(', high_risk_sample_percent, ')') , os_rates_table["High", ], hazard_ratio, pvalue)

surv_group_file <- paste(output, "/survival.group.txt", sep = "")
write.table(result_surv2, surv_group_file, row.names = F, quote = F, sep = '\t')

