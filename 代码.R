####TCGA数据分析代码####

####一####
####TCGA-UCEC 数据下载####
setwd('TCGA-UCEC')
setwd("TCGAdata")
library(tidyverse)#加载包
#安装TCGAbiolinks包 R中没有这个包 要通过外界的包下载需要的包
install.packages("BiocManager")#首先下载R里的BiocManager包
library(BiocManager)#加载包
#通过BiocManager下载我们需要的包
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("remotes")
BiocManager::install("ExperimentHub")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)#加载包
cancer_type = 'TCGA-UCEC'  #肿瘤类型，这里可修改癌症类型
#TCGA 肿瘤缩写：https://www.jianshu.com/p/3c0f74e85825
#首先找到我们需要的数据
expquery <- GDCquery(project = cancer_type,
                     data.category = 'Transcriptome Profiling',
                     data.type = 'Gene Expression Quantification',
                     workflow.type = 'STAR - Counts'
)
#下载找到的数据
GDCdownload(expquery,directory = 'GDCdata') #diectory是目录文件夹
#数据整理，格式更改
expquery2 <- GDCprepare(expquery,directory = 'GDCdata',summarizedExperiment = T)
save(expquery2,file = 'UCEC.gdc_2025.rda')#保存rda格式
#退出Rstudio
setwd('TCGA-UCEC')
setwd('TCGAdata')
library(tidyverse)#加载包
load('UCEC.gdc_2025.rda')#导入文件，rda格式文件也可直接从文件夹
load('gene_annotation_2022.rda')#导入基因注释文件
table(gene_annotation_2022$type)#table 分组jishu分组计数
#基因名称symbol ENSEMBL
#提取counts <- 以下三句无需掌握 tpms
#data <- read.table('E:/TCGAxunlian/TCGA-UCEC/TCGAdata/gene_annotation_2025.csv',header = TRUE,sep = ',',stringsAsFactors = FALSE)
#gene_annotation_2025 <- as.data.frame(data)#导入基因注释文件
#table(gene_annotation_2025$type)
#方法一
counts <- expquery2@assays@data@listData[['unstranded']]
colnames(counts) <- expquery2@colData@rownames
rownames(counts) <- expquery2@rowRanges@ranges@NAMES
#基因ID转换
counts <- counts %>%
  as.data.frame() %>%
  rownames_to_column('ENSEMBL') %>%
  inner_join(gene_annotation_2022,'ENSEMBL') %>%
  .[!duplicated(.$symbol),]

a <- rownames(counts)
b <- rownames(gene_annotation_2022)
identical(a,b)#行数不同，无法运行
counts$ENSEMBL <- as.character(gene_annotation_2022$symbol)
rownames(counts) <- counts$ENSEMBL   #将行名变为Gene Symbol,但是行数不同所以不能用
ncol(gene_annotation_2025)
nrow

#拆解以上代码
#方法二
counts1 <- expquery2@assays@data@listData[['unstranded']]
colnames(counts1) <- expquery2@colData@rownames
rownames(counts1) <- expquery2@rowRanges@ranges@NAMES
counts1 <- as.data.frame(counts1)
counts1 <- rownames_to_column(counts1,var = 'ENSEMBL')#赋予行名变成列 把列变成行名 column_to_rownames
counts1 <- inner_join(counts1,gene_annotation_2022,'ENSEMBL')

#先比对两个表格是否相似，表格行相同才能整合，之后去重复
a <- rownames(counts1)
b <- rownames(gene_annotation_2022)
identical(a,b)
counts1$ENSEMBL <- as.character(gene_annotation_2022$symbol) #新增Gene Symbol
counts1 <- counts1[!duplicated(counts1$symbol),] #去重复
#拆解完毕，counts=counts1
rownames(counts1) <- counts1$ENSEMBL   #将行名变为Gene Symbol
ncol(gene_annotation_2022)
nrow
counts <- counts1[which(counts1$type == "protein_coding"),] #只要编码RNA
counts <- counts[,-ncol(counts)]   #去除最后一列
write.table(counts, file = "COAD_counts_mRNA_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#方法三
#按照提供的程序比对counts和counts1是否相似
counts <- expquery2@assays@data@listData[['unstranded']]
colnames(counts) <- expquery2@colData@rownames
rownames(counts) <- expquery2@rowRanges@ranges@NAMES
counts <- counts %>%
  as.data.frame() %>%
  rownames_to_column('ENSEMBL') %>%
  inner_join(gene_annotation_2022,'ENSEMBL') %>%
  .[!duplicated(.$symbol),]

counts1 <- expquery2@assays@data@listData[['unstranded']]
colnames(counts1) <- expquery2@colData@rownames
rownames(counts1) <- expquery2@rowRanges@ranges@NAMES
counts1 <- as.data.frame(counts1)
counts1 <- rownames_to_column(counts1,var = 'ENSEMBL')#赋予行一个列名 把列变成行名 column_to_rownames
counts1 <- inner_join(counts1,gene_annotation_2022,'ENSEMBL')
counts1 <- counts1[!duplicated(counts1$symbol),] #去重复
#拆解完毕，counts=counts1
a <- c('a','b','a','b','c')
b <- c('a','b','b','a','c')
identical(a,b)#顺序也要相同
identical(colnames(counts),colnames(counts1))#列相同
identical(rownames(counts),rownames(counts1))#行相同
#继续跑
rownames(counts) <- NULL
counts <- counts %>% column_to_rownames('symbol')
#拆解
rownames(counts1) <- NULL
counts1 <- column_to_rownames(counts1,var = 'symbol')
#继续跑
#保留mRNA
#table
table(counts$type)#(注：可通过table（counts$type)查看基因类型数目
counts <- counts[counts$type == "protein_coding",] #只要编码RNA
#counts <- counts[counts$type == "lncRNA",] 

counts <- counts[,-c(1,ncol(counts))]   #去除第一列和最后一列
#ncol是数列数
ncol(counts)
nrow(counts)
#把TCGA barcode（就是列名）切割为16位字符，并去除重复样本
colnames(counts) <- substring(colnames(counts),1,16)#提取列名前16位
counts <- counts[,!duplicated(colnames(counts))]
table(substring(colnames(counts),14,16))#01A是肿瘤样本 11A是正常的样本
#保留01A （注：可用table(substring(colnames(counts)),14,16)）
counts01A <- counts[,substring(colnames(counts),14,16) == c('01A')]
#保存11A
counts11A <- counts[,substring(colnames(counts),14,16) == c('11A')]
table(substring(colnames(counts01A),14,16))
table(substring(colnames(counts11A),14,16))

####tpms####
#和counts基本一模一样
##expquery2的data-listdata的unstranded是counts tpm_unstranded是tpms##
tpms <- expquery2@assays@data@listData[['tpm_unstrand']]
colnames(tpms) <- expquery2@colData@rownames
rownames(tpms) <- expquery2@rowRanges@ranges@NAMES
tpms <- tpms %>%
  as.data.frame() %>%
  rownames_to_column('ENSEMBL') %>%
  inner_join(gene_annotation_2022,'ENSEMBL') %>%
  .[!duplicated(.$symbol),]
rownames(tpms) <- NULL
tpms <- tpms %>% column_to_rownames('symbol')
tpms <- tpms[tpms$type == "protein_coding",] #只要编码RNA
#可用table（tpms$tpye）查看基因类型
tpms <- tpms[,-c(1,ncol(tpms))]   #去除第一列和最后一列
#把TCGA barcode（就是列名）切割为16位字符，并去除重复样本
colnames(tpms) <- substring(colnames(tpms),1,16)#提取列名前16位
tpms <- tpms[,!duplicated(colnames(tpms))]
table(substring(colnames(tpms),14,16))#01A是肿瘤样本 11A是正常的样本
#保留01A （注：可用table(substring(colnames(counts)),14,16)）
tpms01A <- tpms[,substring(colnames(tpms),14,16) == c('01A')]
#保存11A
tpms11A <- tpms[,substring(colnames(tpms),14,16) == c('11A')]
table(substring(colnames(counts01A),14,16))
table(substring(colnames(counts11A),14,16))
#判断counts和tpms的行列名是否一致
identical(rownames(counts01A),rownames(counts11A))
identical(rownames(tpms01A),rownames(tpms11A))
identical(rownames(counts01A),rownames(tpms01A))
identical(rownames(counts11A),rownames(tpms11A))
identical(colnames(counts11A),colnames(tpms11A))
#保存counts和tpms数据
write.table(counts01A,'counts01A.txt',sep = '\t',row.names = T,col.names = NA,quote = F)
write.table(counts11A,'counts11A.txt',sep = '\t',row.names = T,col.names = NA,quote = F)
write.table(tpms01A,'tpms01A.txt',sep = '\t',row.names = T,col.names = NA,quote = F)
write.table(tpms11A,'tpms11A.txt',sep = '\t',row.names = T,col.names = NA,quote = F)
#cbind列合并和rbind行合并
#cbind之前需要确认两个数据框的行名
counts <- cbind(counts01A,counts11A)#前提是行相同才能列合并
tpms <- cbind(tpms01A,tpms11A)
write.table(counts,'counts.txt',sep = '\t',row.names = T,col.names = NA,quote = F)
write.table(tpms,'tpms.txt',sep = '\t',row.names = T,col.names = NA,quote = F)

####tpms_log2####
#counts是差异分析
range(tpms)#查看数据范围
range(tpms01A)
range(tpms11A)
tpms_log2 <- log2(tpms+1)#log2转换 加1是为了数值大于0
range(tpms_log2)
tpms01A_log2 <- log2(tpms01A+1)
range(tpms01A_log2)
tpms11A_log2 <- log2(tpms11A+1)
range(tpms11A_log2)
#保存log2转换后的数据
write.table(tpms_log2,'tpms_log2.txt',sep = '\t',row.names = T,col.names = NA,quote = F)
write.table(tpms01A_log2,'tpms01A_log2.txt',sep = '\t',row.names = T,col.names = NA,quote = F)
write.table(tpms11A_log2,'tpms11A_log2.txt',sep = '\t',row.names = T,col.names = NA,quote = F)
#表达谱整理完毕

####二####
####ESTIMATE####
#计算患者免疫评分与肿瘤纯度# (基质组分+免疫组分+肿瘤组分)=1 肿瘤纯度
setwd('TCGA-UCEC')
setwd('ESTIMATE') #设置工作目录
#安装包
library(utils) #这个包应该不用下载 直接加载试试
#rforge <- "http://r-forge.r-project.org"
#install.packages('estimate',repos=rforge,dependencies=TRUE)
library(estimate)
library(tidyverse)
#读取肿瘤患者01A表达谱
exp <- read.table('tpms01A_log2.txt',sep = '\t',row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#计算免疫评分 CD8 CD4的比值等
filterCommonGenes(input.f = 'tpms01A_log2.txt', #输入文件名
                  output.f = 'tpms01A_log2.gct', #输出文件名
                  id = 'GeneSymbol') #行名为gene symbol
estimateScore('tpms01A_log2.gct', #刚才的输出文件名
              'tpms01A_log2_estimate_score.txt', #新的输出文件名
              platform = 'affymetrix') #默认平台

#提取结果并整理
ESTIMATE_result <- read.table('tpms01A_log2_estimate_score.txt',sep = '\t',row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ESTIMATE_result <- ESTIMATE_result[,-1]
colnames(ESTIMATE_result) <- ESTIMATE_result[1,]
ESTIMATE_result <- as.data.frame(t(ESTIMATE_result[-1,]))
rownames(ESTIMATE_result) <- colnames(exp)
#保存结果
write.table(ESTIMATE_result,file = 'ESTIMATE_result.txt',sep = '\t',row.names = T,col.names = NA,quote = F)

####生存信息整理####
#xena官网：https://xenabrowser.net/datapages/
#下载生存信息
setwd('TCGA-UCEC')
setwd('Survival_data')
library(tidyverse)
#手动导入OS.txt取名survival
survival <- survival[,2:3]
survival$OS.time <- survival$OS.time/365
survival <- survival %>% rownames_to_column('sample')
#新建一列name
survival$name <- paste0(survival$sample,'A') #paste粘贴、连接
table(substring(survival$name,14,16))
rownames(survival) <- survival$name
survival <- survival[,2:3] #OS=1死亡时间发生0是生存
#合并PFI生存信息与表达谱
tpms01A_log2 <- read.table("tpms01A_log2.txt", sep = "\t",row.names = 1,check.names = F,header = T)
#取交集
a <- intersect(colnames(tpms01A_log2),rownames(survival))
table(substr(a,14,16))
exp_01A <- tpms01A_log2[,a]#用a的顺序展示tpms01A_log2
surv_01A <- survival[a,]#用a的顺序展示survival
#为了把exp数据和survival数据合并，需要行列转换
exp_01A <- exp_01A %>% t() %>% as.data.frame()
#合并前判断行名是否相同
identical(rownames(exp_01A),rownames(surv_01A))
#列合并数据
exp_surv_01A <- cbind(surv_01A,exp_01A)
##保存文件##
write.table(exp_surv_01A,"exp_surv_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#合并生存信息与ESTIMATE
ESTIMATE_result <- read.table("ESTIMATE_result.txt", sep = "\t",row.names = 1,check.names = F,header = T)
identical(rownames(ESTIMATE_result),rownames(surv_01A))
a <- intersect(rownames(ESTIMATE_result),rownames(surv_01A))
ESTIMATE_result <- ESTIMATE_result[a,]
ESTIMATE_result_surv_01A <- cbind(surv_01A,ESTIMATE_result)

##保存文件##
write.table(ESTIMATE_result_surv_01A,"ESTIMATE_result_surv_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####根据ESTIMATE_result高低组做生存分析####
setwd("TCGA-UCEC")
setwd("survival")
surv <- read.table("ESTIMATE_result_surv_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#把生存时间变成年为单位
surv$OS.time <- surv$OS.time/365

#方法一
#median 中位数
#ImmuneScore
#新增一列group 探究ImmuneScore免疫组分与生存的关系 用median(surv$ImmuneScore)求出ImmuneScore中位数
surv$group <- ifelse(surv$ImmuneScore > median(surv$ImmuneScore),"High","Low")
class(surv$group)

#将surv$group从字符转变成因子
surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)
#install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group, #将OS.time, OS与group关联
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.050, " < 0.050", paste0(" = ",round(pValue, 3))))

#install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间，扩散的阴影部分
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco，jama,lancet
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,15), # x轴长度
           break.time.by = 5, # x轴步长为5
           legend.title = "ImmuneScore", #标题
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间，扩散的阴影部分
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco，jama,lancet
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0, 15), # x轴长度
           break.time.by = 5, # x轴步长为5
           legend.title = "ImmuneScore", # 图例标题
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE,
           ggtheme = theme_minimal() + # 使用 theme_minimal 作为基础主题
             theme(legend.title = element_text(size = 20))) # 调整图例标题字号为14
dev.off()

#StromalScore
surv$group <- ifelse(surv$StromalScore > median(surv$StromalScore),"High","Low")
surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)
#install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
#install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,15), # x轴长度
           break.time.by = 5, # x轴步长为5
           legend.title = "StromalScore",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE,
           ggtheme = theme_minimal() + # 使用 theme_minimal 作为基础主题
             theme(legend.title = element_text(size = 20))) # 调整图例标题字号为14
dev.off()

#ESTIMATEScore
surv$group <- ifelse(surv$ESTIMATEScore > median(surv$ESTIMATEScore),"High","Low")
surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)
#install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.05, " < 0.05", paste0(" = ",round(pValue, 3))))
#install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,15), # x轴长度
           break.time.by = 5, # x轴步长为5
           legend.title = "ESTIMATEScore",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE,
           ggtheme = theme_minimal() + # 使用 theme_minimal 作为基础主题
             theme(legend.title = element_text(size = 20))) # 调整图例标题字号为14
dev.off()

####三####
####整理TCGA临床信息####
setwd("TCGA-UCEC")
setwd("clinical")
library(tidyverse)
load("UCEC.gdc_2025.rda")
#提取临床信息
clinical <- as.data.frame(expquery2@colData) %>%   
  .[!duplicated(.$sample),] #01A的definition是肿瘤，11A是正常组织。pathologic_stage是临床分析阶段，pathologic_t是T分析，pathologic_n是N分析，pathologic_m是M分析，MX不知道有没有转移，M0没有转移，M1原属转移，gender，vital_status生存状态，age,days_to_death(OS)

#提取需要的临床信息数据
clinical <-clinical[,c("gender","age_at_index","figo_stage",
                       "residual_disease","tumor_grade")]

class(clinical$gender)
class(clinical$age_at_index)
class(clinical$figo_stage)
class(clinical$residual_disease)
class(clinical$tumor_grade)

table(clinical$gender)
table(clinical$age_at_index)
table(clinical$figo_stage)
table(clinical$residual_disease)
table(clinical$tumor_grade)

#gsub()替换 使阶段只用数字区分
clinical$figo_stage <- gsub("A","",clinical$figo_stage)#ajcc_pathologic_stage这一列替换A为空白
clinical$figo_stage <- gsub("B","",clinical$figo_stage)
clinical$figo_stage <- gsub("C","",clinical$figo_stage)
clinical$figo_stage <- gsub("1","",clinical$figo_stage)
clinical$figo_stage <- gsub("2","",clinical$figo_stage)


#提取01A临床数据
rownames(clinical) <- substring(rownames(clinical),1,16)

##将基因表达谱和临床数据合并并保存
exp01A <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

clinical01A <- clinical[colnames(exp01A),] #用exp01A的行名排列方式排clinical01A  
exp01A <- exp01A %>% t() %>% as.data.frame()

identical(rownames(clinical01A),rownames(exp01A))

clinical.expr01A <- cbind(clinical01A,exp01A)

write.table(clinical.expr01A,"clinical.expr01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

##将ESTIMATE_result和临床数据合并并保存
ESTIMATE_result <- read.table("ESTIMATE_result.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

identical(rownames(clinical01A),rownames(ESTIMATE_result))

clinical.ESTIMATE_result01A <- cbind(clinical01A,ESTIMATE_result)
#保存表格用于后面作图
write.csv(clinical.ESTIMATE_result01A,file = "clinical.ESTIMATE_result01A.csv")
#解螺旋：https://www.helixlife.cn/class


####四####
####差异分析####
####免疫评分####
setwd("TCGA-UCEC")
setwd("Immune_DEG") #diff expression gene
library(BiocManager)
#BiocManager::install('DESeq2')
library(DESeq2)
library(tidyverse)
#TCGA差异分析用counts来做 因为是把01A患者分组做差异分析所以读取01A
counts_01A <- read.table("counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#因为是用免疫评分分组所以读取ESTIMATE_result
estimate <- read.table("ESTIMATE_result.txt", sep = "\t",row.names = 1,check.names = F,header = T)
#整理分组信息
x <- "ImmuneScore"
med <- as.numeric(median(estimate[,x])) #as.numeric作为数值
estimate <- as.data.frame(t(estimate)) #行列转换
identical(colnames(counts_01A),colnames(estimate))
#data.frame创造数据框
conditions=data.frame(sample=colnames(counts_01A), #第一列名字是sample
                      group=factor(ifelse(estimate[x,]>med,"high","low"),levels = c("low","high"))) %>% 
  column_to_rownames("sample")
#拆解上句长代码
#conditions=data.frame(sample=colnames(counts_01A),
#                      group=factor(ifelse(estimate[x,]>med,"high","low"),levels = c("low","high")))
#conditions <- column_to_rownames(conditions,"sample")

#差异分析准备工作
dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData = conditions,
  design = ~ group)

#开始差异分析
dds <- DESeq(dds)
#这句很重要
resultsNames(dds) #记住group里是high/low还是low/high
#提取结果
res <- results(dds)
save(res,file="DEG_ImmuneScore.Rda")

####热图绘制####
DEG <- as.data.frame(res)#第二列log2正数为高表达，负数为低表达和最后一列padj调节P值
#读取表达谱
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#添加上下调信息，第二列log2以1为界限划分高低表达，除去不符合标准的数据
logFC_cutoff <- 2
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
#新增一列
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)
#下载pheatmap包 
#install.packages("pheatmap")
library(pheatmap)
#提取差异基因表达谱
a <- filter(DEG,change == 'UP') #filter筛选函数
b <- filter(DEG,change == 'DOWN')
c <- rbind(a,b)
d <- rownames(c)
exp_diff <- exp[d,]
#设置分组信息
annotation_col <- conditions
#对exp_diff 列的顺序进行处理
a <- filter(annotation_col,group == 'high')
b <- filter(annotation_col,group == 'low')
exp_diff_high <- exp_diff[,rownames(a)]
exp_diff_low <- exp_diff[,rownames(b)]
exp_diff <- cbind(exp_diff_high,exp_diff_low)
#开始画图
color_breaks <- seq(1.5,-1.5,length.out = 80)
pheatmap(exp_diff,
         annotation_col=annotation_col,
         scale = 'row',
         main = "ImmuneScore", # 添加标题
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("blue", "white", "red"))(80),
         breaks = color_breaks,
         cluster_cols =F,#列聚类
         cluster_rows = T,
         fontsize = 10,#字体大小
         fontsize_row=12,
         fontsize_col=12)



#保存图片 调整大小
dev.off()#关闭画板

####基质评分####
#不用退出R更改路径方法
current_dir <- getwd()
current_dir
cat('E:/TCGA论文/TCGA-UCEC/Immune_DEG',current_dir,'\n')
parent_dir <- dirname(current_dir)
parent_dir
cat('E:/TCGA论文/TCGA-UCEC',parent_dir,'\n')
setwd("E:/TCGA论文/TCGA-UCEC")
setwd("Stromal_DEG")
library(BiocManager)
library(DESeq2)
library(tidyverse)
#TCGA差异分析用counts来做 因为是把01A患者分组做差异分析所以读取01A
counts_01A <- read.table("counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#因为是用免疫评分分组所以读取ESTIMATE_result
estimate <- read.table("ESTIMATE_result.txt", sep = "\t",row.names = 1,check.names = F,header = T)
#整理分组信息
x <- "StromalScore"
med <- as.numeric(median(estimate[,x]))
estimate <- as.data.frame(t(estimate))
identical(colnames(counts_01A),colnames(estimate))

conditions=data.frame(sample=colnames(counts_01A),
                      group=factor(ifelse(estimate[x,]>med,"high","low"),levels = c("low","high"))) %>% 
  column_to_rownames("sample")
#差异分析准备工作
dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData = conditions,
  design = ~ group)

#开始差异分析
dds <- DESeq(dds)
#这句很重要
resultsNames(dds)
#提取结果
res <- results(dds)
save(res,file="DEG_StromalScore.Rda")

####热图绘制####
DEG <- as.data.frame(res)
#读取表达谱
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#添加上下调信息
logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)

library(pheatmap)
#提取差异基因表达谱
a <- filter(DEG,change == 'UP')
b <- filter(DEG,change == 'DOWN')
c <- rbind(a,b)
d <- rownames(c)
exp_diff <- exp[d,]
#设置分组信息
annotation_col <- conditions
#对exp_diff 列的顺序进行处理
a <- filter(annotation_col,group == 'high')
b <- filter(annotation_col,group == 'low')
exp_diff_high <- exp_diff[,rownames(a)]
exp_diff_low <- exp_diff[,rownames(b)]
exp_diff <- cbind(exp_diff_high,exp_diff_low)
#开始画图
color_breaks <- seq(1.5,-1.5,length.out = 80)
pheatmap(exp_diff,
         annotation_col=annotation_col,
         scale = 'row',
         main = "StromalScore", # 添加标题
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("blue", "white", "red"))(80),
         breaks = color_breaks,
         cluster_cols =F,#列聚类
         cluster_rows = T,
         fontsize = 10,#字体大小
         fontsize_row=12,
         fontsize_col=12)

#保存图片 调整大小
dev.off()#关闭画板

####将两次差异分析的差异基因取交集####
current_dir <- getwd()
current_dir
cat('E:/TCGA论文/TCGA-UCEC/Stromal_DEG',current_dir,'\n')
parent_dir <- dirname(current_dir)
parent_dir
cat('E:/TCGA论文/TCGA-UCEC',parent_dir,'\n')
setwd("E:/TCGA论文/TCGA-UCEC")
setwd("Immune_Stromal_DEG")
#打开DEG_ImmuneScore.rda
DEG <- as.data.frame(res)
#添加上下调信息
logFC_cutoff <- 2
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)
#提取上下调基因
library(tidyverse)
a <- filter(DEG,change == 'UP')
b <- filter(DEG,change == 'DOWN')
write.csv(a, file = "Immune_up.csv")
write.csv(b, file = "Immune_down.csv")

#打开DEG_StromalScore.rda
DEG <- as.data.frame(res)
#添加上下调信息
logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)
#提取上下调基因
a <- filter(DEG,change == 'UP')
b <- filter(DEG,change == 'DOWN')
write.csv(a, file = "Stromal_up.csv")
write.csv(b, file = "Stromal_down.csv")
#仙桃画图

####用交集差异基因做富集分析####
setwd("TCGA-UCEC")
setwd("FUJI_Immune_Stromal_DEG")
library(tidyverse)
library("BiocManager")
#安装加载包
#BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)
#org.Hs.eg.db包主要注释人类基因:用于不同数据库ID间的转化
library(clusterProfiler)
library(DOSE)
#导入DEG_final.txt   UP交集共同基因和DOWN交集共同基因的集合
#导入immune或stromal差异分析结果 均可DEG_ImmuneScore/DEG_StromalScore
DEG <- as.data.frame(res)
DEG <- DEG[DEG_final$SYMBOL,]
DEG <- rownames_to_column(DEG,"SYMBOL")
genelist <- bitr(DEG$SYMBOL, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by="SYMBOL")

####GO####
#GO描述基因功能
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all", #三种形式 BP生物学功能 CC细胞层次 MF分子层次
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
save(ego,ego_res,file = "GO_DEG_final.Rda")

####KEGG####
#KEGG展示基因参与的通路
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa', #人
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kk_res <- kk@result
save(kk,kk_res,file = "KEGG_DEG_final.Rda")

#网络图
library(ggnewscale)
library(enrichplot)
#install.packages("ggnewscale")
List = DEG$log2FoldChange
names(List)= DEG$ENTREZID #用ENTREZID形式命名List
head(List) #查看数据的前几列
List = sort(List,decreasing = T) #sort函数是排序
#GO
ego <- DOSE::setReadable(ego,
                         OrgDb = 'org.Hs.eg.db',
                         keyType = 'ENTREZID')
cnetplot(ego, 
         foldChange = List,  # 添加基因的Fold Change信息
         circular = TRUE,    # 设置为圆形布局
         colorEdge = TRUE,   # 对边进行颜色编码
         node_label = "gene",  # 显示基因标签
         showCategory = 10,  # 显示前10个类别
         layout = "circle",
         categorySize = 'geneNum')  # 尝试使用圆形布局
#KEGG
kk <- DOSE::setReadable(kk,
                        OrgDb = 'org.Hs.eg.db',
                        keyType = 'ENTREZID')
cnetplot(kk, 
         foldChange = List,  # 添加基因的Fold Change信息
         circular = TRUE,    # 设置为圆形布局
         colorEdge = TRUE,   # 对边进行颜色编码
         node_label = "gene",  # 显示基因标签
         showCategory = 10,  # 显示前10个类别
         layout = "circle",
         categorySize = 'geneNum')  # 尝试使用圆形布局

####五####
####PPI####
#蛋白互作网络
#String：https://cn.string-db.org/cgi/input?sessionId=bXmYsv7CnUrH&input_page_active_form=multiple_identifiers
#制作网络-setting-最小的联系值调成最大
#Cytoscape 将String制作的网络图导入，进一步改进网络图导出图片以及右下角表格（重点是Degree蛋白关联度和Name）
#用仙桃学术把Cytoscape导出的表格提取Name和Degree，作一维条形图
####COX####
setwd("TCGA-UCEC")
setwd("COX")
setwd('..')#跳转上一级
#安装加载R包
#install.packages("survival")
#install.packages("forestplot")
library(survival)
library(forestplot)
library(tidyverse)
#读取文件
exp_surv_01A = read.table("exp_surv_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#手动读取DEG_final.txt 
#提取DEG_final
surv.expr <- cbind(exp_surv_01A[,1:2],exp_surv_01A[,DEG_final$SYMBOL])
#a <- exp_surv_01A[,1:2]
#b <- exp_surv_01A[,DEG_final$SYMBOL]
#Cox分析
#如何修改特定列的列名
#colnames(surv.expr)[ ] <- ""  #[]内填特定列数字 ""内填写修改的名字
Coxoutput <- NULL 

for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS.time,OS) ~ surv.expr[,i], data = surv.expr) # 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}
#HR是风险比，HR>1，是危险因素，HR<1，是保护因素
Coxoutput <- arrange(Coxoutput,pvalue) #arrange排序
###筛选top基因
a <- Coxoutput[Coxoutput$pvalue < 0.01,] # 取出p值小于0.05的基因
b <- Coxoutput[Coxoutput$HR < 2,] 
gene_sig <- intersect(a,b)
write.csv(gene_sig, file = "gene_sig.csv")
topgene <- gene_sig #为了下面不改topgene
#3. 绘制森林图
##3.1 输入表格的制作
#将小于0.001的数整合
topgene$pvalue <- round(topgene$pvalue,3)
topgene$pvalue <- ifelse(topgene$pvalue < 0.001,'< 0.001',topgene$pvalue)
tabletext <- cbind(c("Gene",topgene$gene),
                   c("HR",format(round(as.numeric(topgene$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(topgene$lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(topgene$upper),3),nsmall = 3)),
                   c("pvalue",topgene$pvalue))
##3.2 绘制森林图
#动态计算xticks的范围
#min_val <- min(topgene$lower,na.rm = T)
#max_val <- max(topgene$upper,na.rm = T)
#xticks_range <- seq(from = min_val,to = max_val)
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(topgene$HR)),
           lower=c(NA,as.numeric(topgene$lower)), 
           upper=c(NA,as.numeric(topgene$upper)),
           graph.pos=5,# 图在表中的列位置
           graphwidth = unit(.25,"npc"),# 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           xticks = c(0.5,1,1.5),# 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), # 在第一行上面画黑色实线
                           "2" = gpar(lwd=1.5, col="black"), # 在第一行标题行下画黑色实线
                           "50" = gpar(lwd=2, col="black")), # 在最后一行上画黑色实线
           lineheight = unit(.75,"cm"),# 固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
#保存图片 大小30*30
dev.off()

####六####
####CXCL13在肿瘤样本与正常样本中的表达差异####
####柱状图####
setwd("TCGA-UCEC")
setwd("CXCL13")
library(tidyverse)
#分析正常组织和肿瘤组织中CXCL13的表达差异
tpms01A_log2 <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tpms11A_log2 <- read.table("tpms11A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- "CXCL13"#以后修改这里即可
a <- tpms01A_log2[gene,]
b <- tpms11A_log2[gene,]
#t转换
a <- a %>% t() %>% as.data.frame()
b <- b %>% t() %>% as.data.frame()
write.csv(a, file = "CXCL13_01A.csv")
write.csv(b, file = "CXCL13_11A.csv")
#新建表格将01A和11A的数值合并
#用仙桃学术制作分组比较图

####配对图绘制####
tpms01A_log2 <- tpms01A_log2 %>% t() %>% as.data.frame()
tpms11A_log2 <- tpms11A_log2 %>% t() %>% as.data.frame()
rownames(tpms01A_log2) <- substring(rownames(tpms01A_log2),1,12)
rownames(tpms11A_log2) <- substring(rownames(tpms11A_log2),1,12)
a <- intersect(rownames(tpms01A_log2),rownames(tpms11A_log2))
tpms01A_log2 <- tpms01A_log2[a,]
tpms11A_log2 <- tpms11A_log2[a,]
peidui <- cbind(tpms11A_log2[,gene],tpms01A_log2[,gene])#11A放在前面
peidui <- as.data.frame(peidui)
write.csv(peidui,file = "peidui.csv")
#仙桃学术配对图

####根据CXCL13高低组做生存分析####
setwd('..')
#setwd("TCGA-UCEC")
setwd("survival")
surv <- read.table("exp_surv_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv$OS.time <- surv$OS.time/365

#median 中位数
#CXCL13
surv$group <- ifelse(surv$CXCL13 > median(surv$CXCL13),"High","Low")
class(surv$group)
surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)
#install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.05, " < 0.05", paste0(" = ",round(pValue, 3))))
#install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间，扩散的阴影部分
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco，jama,lancet
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0, 20), # x轴长度
           break.time.by = 5, # x轴步长为5
           legend.title = "CXCL13", # 图例标题
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE,
           ggtheme = theme_minimal() + # 使用 theme_minimal 作为基础主题
             theme(legend.title = element_text(size = 20))) # 调整图例标题字号为14
dev.off()


####不同分期CXCL13的表达####
setwd('..')
#setwd("TCGA-UCEC")
setwd("CXCL13")
library(tidyverse)
clinical.expr01A = read.table("clinical.expr01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- "CXCL13"
clinical_CXCL13 <- cbind(clinical.expr01A[,1:6],clinical.expr01A[,gene])
write.csv(clinical_CXCL13, file = "clinical_CXCL13.csv")
#仙桃学术制作分组比较图-箱型图

####七####
####CXCL13差异分析####
setwd("TCGA-UCEC")
setwd("CCL5_DEG")
library(DESeq2)
library(tidyverse)
counts_01A <- read.table("counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- read.table("tpms01A_log2.txt", sep = "\t",row.names = 1,check.names = F,header = T)
identical(colnames(counts_01A),colnames(exp))#习惯性判断以防万一
gene <- "CXCL13"#每次运行只改这个基因名
med=median(as.numeric(exp[gene,]))#取数值后取中位数

conditions=data.frame(sample=colnames(exp),
                      group=factor(ifelse(exp[gene,]>med,"high","low"),levels = c("low","high"))) %>% 
  column_to_rownames("sample")

dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData = conditions,
  design = ~ group)

dds <- DESeq(dds)

resultsNames(dds)
res <- results(dds)
save(res,file="DEG_CXCL13.Rda")

####GSEA####
#安装加载包
#GO KEGG 时已经装过 直接library即可
#BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db) #org.Hs.eg.db包主要注释基因:用于不同数据库ID间的转化
library(clusterProfiler)
DEG <- as.data.frame(res)%>% arrange(padj) #根据P值排序

DEG <- DEG %>% rownames_to_column("Gene")

geneList = DEG[,3]
names(geneList) = as.character(DEG[,'Gene'])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)#排序
head(geneList)

#GSEA基因集：https://zhuanlan.zhihu.com/p/504101161
#H C2 C5 C7
msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <- "h.all.v7.0.symbols.gmt"    
#读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

set.seed(1) #设置种子
gsea <-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
#####先HALLMARK跑####
#转换成数据框
gsea_result_df <- as.data.frame(gsea)
save(gsea,gsea_result_df,file = "GSEA_CXCL13_h.all.rda")
#绘图
#安装enrichplot
library(enrichplot)
#单个结果绘制
gseaplot2(gsea,14,color="red",pvalue_table = T)#’1‘表示表格第几列
#多个结果绘制
#A
gseaplot2(gsea, geneSetID = c(1,2,3,4,5,6,8,10), subplots = 1:3,pvalue_table = T, base_size = 30)
#B
gseaplot2(gsea, geneSetID = c(7,9,11,13,14), subplots = 1:3)
# 增加图形高度，为图例留出更多空间
gseaplot2(gsea, geneSetID = 1:7, subplots = 1:3,pvalue_table = T, base_size = 30)    # 调整基础字体大小
gseaplot2(gsea, geneSetID = 8:14, subplots = 1:3,pvalue_table = T, base_size = 30)
dev.off()

#####换C7跑####
msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <- "c7.all.v7.0.symbols.gmt"    
#读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

set.seed(1) #设置种子
gsea <-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
#转换成数据框
gsea_result_df <- as.data.frame(gsea)
save(gsea,gsea_result_df,file = "GSEA_CXCL13_c7.rda")
#绘图
gseaplot2(gsea,2285,color="red",pvalue_table = T)
#C前7个最显著的正富集通路
gseaplot2(gsea, geneSetID = c(377,357,349,218,83,73,63), subplots = 1:3,pvalue_table = T, base_size = 30)
#D前7个最显著的负富集通路
#gseaplot2(gsea,782,color="red",pvalue_table = T,base_size = 40)
gseaplot2(gsea, geneSetID = c(1251,1292,1360,1388,1401,1424,1448), subplots = 1:3,pvalue_table = T, base_size = 30)
#高峰是上调基因 每一个柱状图是加1 空格是减1
dev.off()

####八####
####cibersort####
setwd("TCGA-UCEC")
setwd("CIBERSORT")   
#install.packages('e1071')
#install.packages('parallel')
#BiocManager::install("preprocessCore", version = "3.20")
library(e1071)
library(parallel)
library(preprocessCore)
library(tidyverse)
source("CIBERSORT.R")   
sig_matrix <- "LM22.txt"   
mixture_file = 'tpms01A_log2.txt'   #肿瘤患者表达谱
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
res_cibersort <- res_cibersort[,1:22]   #去除后三列
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞
ciber.res <- as.data.frame(ciber.res)
write.table(ciber.res,"ciber.res.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####cibersort彩虹图 代码无需掌握####
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, # 
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.7, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   #关闭画板

#分组比较图
a <- ciber.res
#读取肿瘤患者01A表达谱
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
med=median(as.numeric(exp["CXCL13",]))
exp <- exp %>% t() %>% as.data.frame()
exp <- exp %>% mutate(group=factor(ifelse(exp$CXCL13>med,"high","low"),levels = c("low","high")))
class(exp$group)
identical(rownames(a),rownames(exp))
a$group <- exp$group
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=CIBERSORT,value = Fraction,-c(group,sample))
ggboxplot(b, x = "CIBERSORT", y = "Fraction",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()


####相关性热图####
#install.packages("ggstatsplot")
#install.packages("ggcorrplot")
#install.packages("corrplot")
library(ggstatsplot)
library(ggcorrplot)
library(corrplot)

cor<-sapply(ciber.res,function(x,y) cor(x,y,method="spearman"),ciber.res)#相关性
rownames(cor)<-colnames(ciber.res)

ggcorrplot(cor, 
           hc.order = TRUE, #使用hc.order进行
           type = "upper", #图片位置upper上方，lower下方
           outline.color = "white",#轮廓颜色
           lab = TRUE,#true为在图上添加相关系数
           ggtheme = ggplot2::theme_gray, #指ggplot2函数对象，默认值为thememinimal
           colors = c("#01468b", "white", "#ee0000"))

####基因与cibersort相关性散点图####
exp = read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- exp["CXCL13",]
ciber = read.table("ciber.res.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber <- ciber %>% t() %>% as.data.frame()
rownames(ciber) <- gsub(" ",".",rownames(ciber))#将空格变成点
identical(colnames(ciber),colnames(exp))
exp_ciber <- rbind(exp,ciber)
exp_ciber <- exp_ciber %>% t() %>% as.data.frame()

library(ggstatsplot)
library(ggside)
ggscatterstats(data = exp_ciber, #要分析的数据
               y = CXCL13, #设置Y轴
               x = B.cells.naive,#设置X轴，可以改名字
               type = "nonparametric", 
               margins = "both",#是否显示 边缘，默认为true                                      
               xfill = "#01468b", #x轴边缘图形的颜色
               yfill = "#ee0000", #y轴边缘图形的颜色
               marginal.type = "densigram")#在图片坐标轴边缘添加图形类型







####单细胞测序分析代码####

####一、安装包####
#安装这些包，加载
library(rjson)
library(limma)   
library(GEOquery) 
library(stringr)
library(survival)
library(glmnet)
library(survminer)
library(timeROC)
library(data.table)
library(ggpubr)
library(dplyr)
library(patchwork)
library(Matrix)
library(readr)
library(tibble)
library(ggplot2)
library(tidyverse) 
library(future)
library(pheatmap)
library(msigdbr)
library(clusterProfiler)
library(devtools)
library(Seurat)
library(glmGamPoi)
library(SingleR)
library(harmony)
library(DoubletFinder)
library(copykat)
library(GSVA)
library(AUCell)
library(monocle)
library(CellChat)
library(SCENIC)
library(RColorBrewer)
library(hdf5r)
install.packages('rjson')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
install.packages('limma')

install.packages('GEOquery')
install.packages('glmnet')
install.packages('timeROC')
install.packages('glmGamPoi')
install.packages('SingleR')
install.packages('harmony')
install.packages('DoubletFinder')
install.packages('copykat')
install.packages('AUCell')
install.packages('monocle')

install.packages('devtools')
library('devtools')
devtools::install_github('sqjin/CellChat')
library('CellChat')

install.packages('SCENIC')
install.packages('hdf5r')




options("repos"="https://mirrors.ustc.edu.cn/CRAN/")
if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")

cran_packages <- c('Matrix',
                   'tibble',
                   'dplyr',
                   'stringr',
                   'ggplot2',
                   'ggpubr',
                   "ggrepel",
                   "ggsci",
                   "gplots",
                   'factoextra',
                   'FactoMineR',
                   'devtools',
                   'cowplot',
                   'patchwork',
                   "pheatmap",
                   'basetheme',
                   'paletteer',
                   'AnnoProbe',
                   'ggthemes',
                   'VennDiagram',
                   'tinyarray') 

Biocductor_packages <- c('ReactomePA',
                         'COSG',
                         "Seurat",
                         'EnhancedVolcano',
                         "Seurat",
                         "TENxucecData",
                         "GSEABase",
                         "GSVA",
                         "clusterProfiler",
                         "org.Hs.eg.db",
                         "UpSetR",
                         "clustree",
                         "conos",
                         "cowplot",
                         "dorothea",
                         "entropy",
                         "future",
                         "msigdbr",
                         "pagoda2",
                         "scRNAseq",
                         "scRNAstat",
                         "tidyverse",
                         "viper",
                         "progeny",
                         "preprocesucecre",
                         "enrichplot")

for (pkg in cran_packages){
  if (! require(pkg,character.only=T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}


for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}

#前面的所有提示和报错都先不要管。主要看这里
for (pkg in c(Biocductor_packages,cran_packages)){
  require(pkg,character.only=T) 
}
library(limma)
install.packages('limma')
#没有任何提示就是成功了，如果有warning xx包不存在，用library检查一下。

#library报错，就单独安装。

####二、运行Seurat计算####
#加载
library(Seurat)

#工作目录
setwd("F:/4scRNA/1code/2Seurat")

#读入
ucec = readRDS("1ucec.rds")

#要把Seurat对象，当做一个数据库，包含多种数据
#其中，最重要的就是包含了多个表达矩阵和细胞注释信息（类似于临床信息）
#当需要用到表达矩阵或细胞注释信息时，会使用默认信息
#不想要使用默认信息，就必须指定或修改默认信息

#用@,$符号依次取下层，也可以[[]]

#assays 储存着表达矩阵
#counts 存储原始数据，是稀疏矩阵
#data 存储Normalize() 规范化的data
#scale.data 存储 ScaleData()缩放后的data
#SCT 储存SCT标准化之后的data

#meta.data中储存着细胞注释信息（类似于临床信息）

#active.assay 储存着默认的矩阵名

#active.ident 储存着默认的细胞注释信息（类似于临床信息）

dim(ucec)
#获取表达矩阵
a = as.matrix(GetAssayData(object = ucec@assays$RNA, layer = "counts")[1:20,1:20])
#或者用slot表示layer
b = as.matrix(GetAssayData(object = ucec@assays$RNA, layer = "data")[1:20,1:20])
c = as.matrix(GetAssayData(object = ucec@assays$SCT, layer = "counts")[1:20,1:20])
d = as.matrix(GetAssayData(object = ucec@assays$SCT, layer = "data")[1:20,1:20])

#修改默认的表达矩阵
DefaultAssay(ucec)
#DefaultAssay(ucec) = "RNA"
e = as.matrix(GetAssayData(object = ucec)[1:20,1:20])

table(d == e)


#获取细胞注释信息（类似于临床信息）
f = ucec@meta.data

table(ucec@meta.data$orig.ident)
colnames(ucec@meta.data)

#修改默认的细胞注释信息
Idents(ucec)
#Idents(ucec) = "nCount_RNA"


#总之，要把Seurat对象，当做一个数据库，包含多种数据
#其中，最重要的就是包含了多个表达矩阵和细胞注释信息（类似于临床信息）
#当需要用到表达矩阵或细胞注释信息时，会使用默认信息
#不想要使用默认信息，就必须指定或修改默认信息


####三、不同格式文件读取数据####
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)

#工作目录
setwd(".../3Read")


#####1.matrix.mtx、genes.tsv和barcodes.tsv####

#目录
setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

# 获取数据文件夹下的所有样本文件列表
samples <- list.files("seurat/")
samples

# 创建一个空的列表
seurat_list <- list()

#读取数据并创建Seurat对象
#R只能读取gz文件
#setwd("D:/Rstudio/生信幻想家单细胞测序学习/3Read/GSE173682/seurat/GSM5276933")
#R.utils::gzip("barcodes.tsv")  # 需要安装 R.utils 包
#R.utils::gzip("features.tsv")
#R.utils::gzip("matrix.mtx")
#setwd('..')
#setwd("D:/Rstudio/生信幻想家单细胞测序学习/3Read/GSE173682/seurat/GSM5276934")
#R.utils::gzip("barcodes.tsv")  # 需要安装 R.utils 包
#R.utils::gzip("features.tsv")
#R.utils::gzip("matrix.mtx")
#setwd('..')
#setwd("D:/Rstudio/生信幻想家单细胞测序学习/3Read/GSE173682/seurat/GSM5276935")
#R.utils::gzip("barcodes.tsv")  # 需要安装 R.utils 包
#R.utils::gzip("features.tsv")
#R.utils::gzip("matrix.mtx")
#setwd('..')
#setwd("D:/Rstudio/生信幻想家单细胞测序学习/3Read/GSE173682/seurat/GSM5276936")
#R.utils::gzip("barcodes.tsv")  # 需要安装 R.utils 包
#R.utils::gzip("features.tsv")
#R.utils::gzip("matrix.mtx")
#setwd('..')
#setwd("D:/Rstudio/生信幻想家单细胞测序学习/3Read/GSE173682/seurat/GSM5276937")
#R.utils::gzip("barcodes.tsv")  # 需要安装 R.utils 包
#R.utils::gzip("features.tsv")
#R.utils::gzip("matrix.mtx")
#setwd('..')
#setwd("D:/Rstudio/生信幻想家单细胞测序学习/3Read/GSE173682/seurat/GSM5276938")
#R.utils::gzip("barcodes.tsv")  # 需要安装 R.utils 包
#R.utils::gzip("features.tsv")
#R.utils::gzip("matrix.mtx")
#setwd('..')
#setwd("D:/Rstudio/生信幻想家单细胞测序学习/3Read/GSE173682/seurat/GSM5276939")
#R.utils::gzip("barcodes.tsv")  # 需要安装 R.utils 包
#R.utils::gzip("features.tsv")
#R.utils::gzip("matrix.mtx")
#setwd('..')
#setwd("D:/Rstudio/生信幻想家单细胞测序学习/3Read/GSE173682/seurat/GSM5276940")
#R.utils::gzip("barcodes.tsv")  # 需要安装 R.utils 包
#R.utils::gzip("features.tsv")
#R.utils::gzip("matrix.mtx")
#setwd('..')
#setwd("D:/Rstudio/生信幻想家单细胞测序学习/3Read/GSE173682/seurat/GSM5276941")
#R.utils::gzip("barcodes.tsv")  # 需要安装 R.utils 包
#R.utils::gzip("features.tsv")
#R.utils::gzip("matrix.mtx")
#setwd('..')
#setwd("D:/Rstudio/生信幻想家单细胞测序学习/3Read/GSE173682/seurat/GSM5276942")
#R.utils::gzip("barcodes.tsv")  # 需要安装 R.utils 包
#R.utils::gzip("features.tsv")
#R.utils::gzip("matrix.mtx")
#setwd('..')
#setwd("D:/Rstudio/生信幻想家单细胞测序学习/3Read/GSE173682/seurat/GSM5276943")
#R.utils::gzip("barcodes.tsv")  # 需要安装 R.utils 包
#R.utils::gzip("features.tsv")
#R.utils::gzip("matrix.mtx")
#setwd('..')
#删除，小于200个基因表达的细胞，小于3个细胞表达的基因
for (sample in samples) {
  #文件路径
  data.path <- paste0("seurat/", sample)
  
  #读取10x数据
  seurat_data <- Read10X(data.dir = data.path)
  
  #创建Seurat对象
  seurat_obj <- CreateSeuratObject(counts = seurat_data,project = sample,min.features = 200,min.cells = 3)
  
  #添加到列表中
  seurat_list <- append(seurat_list, seurat_obj)
}

#合并
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)

#将Layers融合
ucec = JoinLayers(seurat_combined)


##添加细胞注释信息

#获取细胞注释信息
abc789 = ucec@meta.data

#导出
write.table(data.frame(ID=rownames(abc789),abc789),file="meta.txt", sep="\t", quote=F, row.names = F,col.names = T)

#读入
meta = fread("meta.xlsx")

#将第一列转化为行名
meta <-  column_to_rownames(meta,"ID")

#添加meta.data信息，细胞名的顺序必须一致
ucec <- AddMetaData(object = ucec, 
                    metadata = meta,   
                    col.name = c("group1","Type2")) 

#当不一致时，根据seurat对象的名字进行排序
meta = meta[colnames(ucec),]


####四、QC####

library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)

#工作目录
setwd(".../3Read")


#####1.matrix.mtx、genes.tsv和barcodes.tsv####

#目录
setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

# 获取数据文件夹下的所有样本文件列表
samples <- list.files("seurat/")
samples

# 创建一个空的列表
seurat_list <- list()
#删除，小于200个基因表达的细胞，小于3个细胞表达的基因
for (sample in samples) {
  #文件路径
  data.path <- paste0("seurat/", sample)
  
  #读取10x数据
  seurat_data <- Read10X(data.dir = data.path)
  
  #创建Seurat对象
  seurat_obj <- CreateSeuratObject(counts = seurat_data,project = sample,min.features = 200,min.cells = 3)
  
  #添加到列表中
  seurat_list <- append(seurat_list, seurat_obj)
}

#合并
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)

#将Layers融合
ucec = JoinLayers(seurat_combined)





#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)

#目录
setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

#一定要用fread，否则其他函数要很久
#ucec <-fread("GSE173682_series_matrix.txt", sep="\t")
#ucec[1:5,1:5]

#直接转化为行名
#ucec <-  column_to_rownames(ucec,"V1")
#ucec[1:5,1:5]

#创建seurat对象
#ucec <- CreateSeuratObject(ucec, min.features = 200, min.cells = 3)

#将Layers融合
#ucec = JoinLayers(ucec)

#查看
table(ucec@meta.data$orig.ident)

#table(str_split(colnames(ucec),'-',simplify = T)[,2])
#增加或修改meta.dada信息
#ucec <- AddMetaData(object = ucec, 
#metadata = str_split(colnames(ucec),'-',simplify = T)[,2],   
#col.name = "orig.ident") 

#table(ucec@meta.data$orig.ident)

#%in% 判断前面一个向量内的元素是否在后面一个向量中，返回布尔值。
table(ucec@meta.data$orig.ident %in% c("GSM5276933","GSM5276934","GSM5276935",'GSM5276936','GSM5276937','GSM5276938','GSM5276939','GSM5276940','GSM5276941','GSM5276942','GSM5276943'))

dim(ucec)

#提取
ucec = subset(ucec,orig.ident %in% c("GSM5276933","GSM5276934","GSM5276935",'GSM5276936','GSM5276937','GSM5276938','GSM5276939','GSM5276940','GSM5276941','GSM5276942','GSM5276943'))

dim(ucec)

#查看
table(ucec@meta.data$orig.ident)

#质量控制
#线粒体基因比例
ucec[["percent.mt"]] <- PercentageFeatureSet(ucec, pattern = "^MT-")

#红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(ucec))
ucec[["percent.HB"]]<-PercentageFeatureSet(ucec, features=HB.genes) 

#查看相关性
FeatureScatter(ucec, "nCount_RNA", "percent.mt", group.by = "orig.ident")
FeatureScatter(ucec, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident")

#查看质控指标
#设置绘图元素
theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
group = "orig.ident"
#质控前小提琴图
plots = list()
for(i in c(1:length(plot.featrures))){
  plots[[i]] = VlnPlot(ucec, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)  
violin
#保存
ggsave("1vlnplot_before_qc.pdf", plot = violin, width = 14, height = 8) 
dim(ucec)

#设置质控指标
quantile(ucec$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(ucec$nFeature_RNA, seq(0.9, 1, 0.01))
#plots[[1]] + geom_hline(yintercept = 500) + geom_hline(yintercept = 4500)
quantile(ucec$nCount_RNA, seq(0.01, 0.1, 0.01))
quantile(ucec$nCount_RNA, seq(0.9, 1, 0.01))
#plots[[2]] + geom_hline(yintercept = 22000)
quantile(ucec$percent.mt, seq(0.9, 1, 0.01))
#plots[[3]] + geom_hline(yintercept = 20)
quantile(ucec$percent.HB, seq(0.9, 1, 0.01))
#plots[[4]] + geom_hline(yintercept = 1)

#设置质控标准
#基因
minGene=300
maxGene=10000
#counts
minUMI=600
#线粒体
pctMT=10
#血细胞
pctHB=1

#数据质控并绘制小提琴图
ucec <- subset(ucec, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
                 nCount_RNA > minUMI & percent.mt < pctMT & percent.HB < pctHB)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(ucec, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=3)    
violin

#保存
ggsave("2vlnplot_after_qc.pdf", plot = violin, width = 14, height = 8) 
dim(ucec)

#保存
saveRDS(ucec,"1ucec_qc.rds")

#读取
ucec = readRDS("1ucec_qc.rds")

####五、double####

#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
# 1. 安装Bioconductor管理器（如果你还没有安装的话）
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

# 2. 通过BiocManager安装关键的生物信息学依赖包
#BiocManager::install(c("Seurat", "ROCR", "matrixStats"))
# 注意：Seurat本身也有大量依赖，安装可能需要一些时间。

# 3. 安装CRAN上的依赖包
#install.packages(c("fields", "KernSmooth", "purrr", "RcppAnnoy"))

# 4. 安装开发工具包devtools，用于从GitHub安装
#install.packages("devtools")

# 5. 从GitHub加载devtools并安装DoubletFinder
#library(devtools)
#install_github('chris-mcginnis-ucsf/DoubletFinder')

# 安装依赖包
install.packages(c("remotes", "Matrix", "KernSmooth", "ROCR", "fields"))

# 从 GitHub 安装最新版本
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)

#目录
setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

#读取
ucec = readRDS("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682/1ucec_qc.rds")

#%>%（向右操作符，forward-pipe operator）是最常用的一种操作符，
#就是把左侧准备的数据或表达式，传递给右侧的函数调用或表达式进行运行，
#可以连续操作就像一个链条一样。

#数据标准化、寻找可变特征、数据缩放
ucec = ucec %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

#PCA、UMAP、寻找邻居、聚类
ucec = ucec %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>%  
  RunTSNE(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.1)

#首先获得最佳的pK值
#pK表示领域大小
sweep.res.list <- paramSweep(ucec, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats) #作图
pk_best = bcmvn %>% 
  dplyr::arrange(desc(BCmetric)) %>% 
  dplyr::pull(pK) %>% 
  .[1] %>% as.character() %>% as.numeric()

#然后估算出双细胞群中，homotypic doublets的比例
annotations <- ucec$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
print(homotypic.prop)

#双细胞占比为7%左右
nExp_poi <- round(0.07*nrow(ucec@meta.data))        
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 

#模拟出的artificial doublet数量。不同取值对识别结果影响不大，默认为0.25
ucec <- doubletFinder(ucec, PCs = 1:30, 
                      pN = 0.25, pK = pk_best, nExp = nExp_poi.adj, 
                      sct = FALSE)

#将列名改为"Double_ucecre"和"Is_Double"
colnames(ucec@meta.data)

colnames(ucec@meta.data)[length(colnames(ucec@meta.data))-1] <- "Double_ucecre"
colnames(ucec@meta.data)[length(colnames(ucec@meta.data))] <- "Is_Double"

#查看DoubletFinder分析结果
head(ucec@meta.data[, c("Double_ucecre", "Is_Double")])

#绘制DoubletFinder分类的tsne图
pdf(file="Double 3.pdf",width=7,height=6)
DimPlot(ucec,reduction = "tsne",label = F,group.by = "Is_Double")
dev.off()
pdf(file="Double 4.pdf",width=7,height=6)
DimPlot(ucec,reduction = "umap",label = F,group.by = "Is_Double")
dev.off()

#绘制双细胞分类的小提琴图
VlnPlot(ucec, group.by = "Is_Double", 
        features = c("nCount_RNA", "nFeature_RNA"), 
        pt.size = 0, ncol = 2)

#过滤非单细胞数据
#ucec <- subset(ucec, Is_Double == "Singlet")

#保存
saveRDS(ucec,"2ucec_double.rds")


####六、CellCycle####

#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

#目录
setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

#读取
ucec = readRDS("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682/2ucec_double.rds")

#细胞周期评分
ucec <- NormalizeData(ucec)

#获取G2M期相关基因
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(ucec))

#获取S期相关基因
s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(ucec))

#细胞周期阶段评分
ucec <- CellCycleucecring(ucec, g2m.features=g2m_genes, s.features=s_genes)

colnames(ucec@meta.data)
table(ucec$Phase)

#画图
DimPlot(ucec,group.by = "Phase",reduction = "tsne")

#保存
saveRDS(ucec,"3ucec_CellCycle.rds")

####七、Normalize####

#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

#目录
setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

#读取
ucec = readRDS("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682/3ucec_CellCycle.rds")

# 每个细胞最初包含相同数量的 RNA 分子的假设。
# 结果存储在ucec@assays$RNA@layers$data
ucec <- NormalizeData(ucec, normalization.method = "LogNormalize", scale.factor = 10000)

# 由于单细胞是稀疏矩阵，很多基因的表达值几乎为0
# 有助于突出单细胞数据集中的生物信号
# 寻找2000个高变基因
ucec <- FindVariableFeatures(ucec, selection.method = "vst", nfeatures = 2000)

# 改变每个基因的表达，使细胞间的平均表达为 0
# 缩放每个基因的表达，使细胞间的方差为 1
# 在下游分析中给予同等权重。
# 在这里去除细胞周期的影响
# 结果存储在ucec@assays$RNA@layers$scale.data
ucec <- ScaleData(ucec,vars.to.regress = c("S.ucecre", "G2M.ucecre"))
#ucec <- ScaleData(ucec,vars.to.regress = c("S.ucecre", "G2M.ucecre"),features = rownames(ucec))

# SCT，相当于替代了上述的三个函数NormalizeData，FindVariable，ScaleData
# 寻找3000个高变基因

# 在常规分析中，使用少量的PC既能关注到关键的生物学差异，
# 又能够不引入更多的技术差异，相当于一种保守性的做法。
# 它会失去一些生物差异信息但是同时又在常规手段中比较安全。
# 但SCT的归一化、标准化都做得不错，
# 多输入一些PCs能提取更多的生物差异，并且兼顾不引入技术误差。
# SCT认为:新增加的这1000个基因就包含了之前没有检测到的微弱的生物学差异。
# 而且，即使使用全部的全部的基因去做下游分析，得到的结果也是和SCT的结果类似

# 结果存储在ucec@assays$SCT中
ucec <- SCTransform(ucec, vars.to.regress = c("S.ucecre", "G2M.ucecre"))

#默认的矩阵
DefaultAssay(ucec)
DefaultAssay(ucec) = "RNA"
DefaultAssay(ucec)

#后续还是以SCT为例进行
DefaultAssay(ucec) = "SCT"
DefaultAssay(ucec)

#保存
saveRDS(ucec,"4ucec_Normalize.rds")

####八、harmony####

library(reticulate)
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ROGUE)
library(clustree)
library(harmony)

#目录
setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

#读取
ucec = readRDS("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682/4ucec_Normalize.rds")

#设置激活的矩阵及分组信息
DefaultAssay(ucec)
DefaultAssay(ucec) = "SCT"
DefaultAssay(ucec)

table(Idents(ucec))
Idents(ucec) = "orig.ident"
table(Idents(ucec))

cxcl13_gene <- "CXCL13"   # 确保大小写与表达矩阵行名一致
cxcl13_in_hv <- "CXCL13" %in% VariableFeatures(ucec)
print(cxcl13_in_hv)   # TRUE 表示 CXCL13 被选为高变基因
features = "CXCL13"

#10个变化最大的基因
top10 <- head(VariableFeatures(ucec), 10)
pdf(file="2.1.pdf",width=7,height=6)
LabelPoints(plot = VariableFeaturePlot(object = ucec), points = top10, repel = TRUE)
dev.off()
#画图
pdf(file="1.pdf",width=7,height=6)
VariableFeaturePlot(object = ucec)
dev.off()

library(ggrepel)

p  <- VariableFeaturePlot(ucec)
cxcl13_df <- p$data["CXCL13", , drop = FALSE]

## 当前对数坐标真实值
x0 <- cxcl13_df$gmean
y0 <- cxcl13_df$residual_variance

## 轻轻往左上推 0.15 轴单位（对数坐标下已很小）
deltaX <- -1.5
deltaY <- +20

p2 <- p +
  geom_text_repel(
    data        = cxcl13_df,
    aes(x = gmean, y = residual_variance, label = "CXCL13"),
    colour      = "black",
    fontface    = "bold",
    size        = 5,
    segment.colour = "black",
    segment.size   = 0.5,
    force       = 0,          # 关掉自动排斥
    nudge_x     = deltaX,
    nudge_y     = deltaY,
    direction   = "both",
    box.padding = 0.25,
    point.padding = 0.25
  )

pdf("2.pdf", width = 7, height = 6)
print(p2)
dev.off()

#PCA
ucec <- RunPCA(ucec, verbose = F)

#主成分分析图形
pdf(file="3.pdf",width=7,height=6)
DimPlot(object = ucec, reduction = "pca")
dev.off()

#绘制每个PCA成分的相关基因
pdf(file="4.pdf",width=10,height=9)
VizDimLoadings(object = ucec, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

#主成分分析热图
pdf(file="5.pdf",width=10,height=9)
DimHeatmap(object = ucec, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#选取合适的PC
#主成分累积贡献大于90%,选择拐点
pdf(file="6.pdf",width=7,height=6)
ElbowPlot(ucec, ndims = 50)
dev.off()

#确定与每个 PC 的百分比   
pct <- ucec [["pca"]]@stdev / sum( ucec [["pca"]]@stdev) * 100
pct

#计算每个 PC 的累计百分比
cumu <- cumsum(pct)
cumu

#设置PC
pcs = 1:40

#harmony
ucec <- RunHarmony(ucec, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)

table(ucec@meta.data$orig.ident)

#选取合适的分辨率
#从0.1-2的resolution结果均运行一遍
seq = seq(0.1,2,by=0.1)
ucec <- FindNeighbors(ucec,  dims = pcs) 
for (res in seq){
  ucec = FindClusters(ucec, resolution = res)
}

#画图
p1 = clustree(ucec,prefix = "SCT_snn_res.")+coord_flip()
p = p1+plot_layout(widths = c(3,1))
ggsave("SCT_sun_res.png", p, width = 30, height = 14)

#降维聚类
# ucec <- FindNeighbors(ucec, reduction = "harmony",  dims = pcs) %>% FindClusters(resolution = 1)
# ucec <- RunUMAP(ucec, reduction = "harmony",  dims = pcs) %>% RunTSNE(dims = pcs, reduction = "harmony")
#降维聚类
ucec <- FindNeighbors(ucec, reduction = "pca",  dims = pcs) %>% FindClusters(resolution = 1)
ucec <- RunUMAP(ucec, reduction = "pca",  dims = pcs) %>% RunTSNE(dims = pcs, reduction = "pca")

colnames(ucec@meta.data)
#画图
pdf(file="7.pdf",width=7,height=6)
DimPlot(ucec, reduction = "umap", label = T)
dev.off()

pdf(file="8.pdf",width=7,height=6)
DimPlot(ucec,reduction = "umap",label = F,group.by = "orig.ident")
dev.off()

pdf(file="9.pdf",width=7,height=6)
DimPlot(ucec,reduction = "umap",label = F,group.by = "Is_Double")
dev.off()

pdf(file="10.pdf",width=7,height=6)
DimPlot(ucec, reduction = "tsne", label = T)
dev.off()

pdf(file="11.pdf",width=7,height=6)
DimPlot(ucec,reduction = "tsne",label = F,group.by = "orig.ident")
dev.off()

pdf(file="12.pdf",width=7,height=6)
DimPlot(ucec,reduction = "tsne",label = F,group.by = "Is_Double")
dev.off()

#保存
saveRDS(ucec,"5ucec_UMPA.TSNE.rds")

####九、SingleR####

#加载
#install.packages('matrixStats')
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
#devtools::install_github("PaulingLiu/ROGUE")
library(ROGUE)
library(clustree)
library(harmony)
library(BiocManager)
#BiocManager::install("SingleR")
library(SingleR)
#SingleR version 2.4.1
#remove.packages("SingleR")
#BiocManager::install("SingleR")
#目录
setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

#读取
ucec = readRDS("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682/5ucec_UMPA.TSNE.rds")

#加载参考数据库
load("ref_Human_all.RData")

#获取表达矩阵
testdata = GetAssayData(object = ucec@assays$RNA, layer = "counts")  # 注意：layer 改为 slot
#testdata = GetAssayData(object = ucec@assays$SCT, layer = "data")

#获取clusters
clusters <- ucec@meta.data$seurat_clusters

#参考数据库的大类和小类细胞群
table(ref_Human_all@colData@listData[["label.main"]])
table(ref_Human_all@colData@listData[["label.fine"]])

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("scrapper")
#运行singleR
cellpred <- SingleR(test = testdata, ref = ref_Human_all, clusters = clusters, assay.type.test = "logcounts", 
                    labels = ref_Human_all@colData@listData[["label.main"]], assay.type.ref = "logcounts")
# 检查包是否在CRAN上
#available_packages <- available.packages()
#"scrapper" %in% rownames(available_packages)
# 搜索相关包
#apt <- available.packages()
#grep("scrap", rownames(apt), value = TRUE, ignore.case = TRUE)

#获取每个cluster的细胞类型
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

#添加到seurat.metadata对象中
ucec@meta.data$SingleR = "NA"
for(i in 1:nrow(celltype)){
  ucec@meta.data[which(ucec$seurat_clusters == celltype$ClusterID[i]),'SingleR'] <- celltype$celltype[i]
}

#install.packages('ComplexHeatmapComplexHeatmap')
library(ComplexHeatmap)
#画图
p = plotScoreHeatmap(cellpred)
ggsave("1.pdf", p, width = 12, height = 5)

p1 <- DimPlot(ucec, group.by = "SingleR", label = T,reduction = "tsne")
p2 <- DimPlot(ucec, group.by = "SingleR", label = T,reduction = "umap")
p <- p1 | p2
p
ggsave("2.pdf", p, width = 12, height = 5)

#保存
saveRDS(ucec,"6ucec_SingleR.rds")

####十、annotation####

#####1. 所有细胞分群#####
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ROGUE)
library(clustree)
library(harmony)
library(SingleR)
library(dplyr)
library(RColorBrewer)
#目录
setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

#读取
ucec = readRDS("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682/6ucec_SingleR.rds")

DefaultAssay(ucec)
#DefaultAssay(ucec) = "RNA"
table(Idents(ucec))
#Idents(ucec) = "nCount_RNA"

#寻找marker
ucec.markers1 <- FindAllMarkers(ucec, only.pos = TRUE,logfc.threshold = 1,)
write.csv(ucec.markers1,file="markers.1.SCT.csv")

DefaultAssay(ucec) = "RNA"
ucec.markers2 <- FindAllMarkers(ucec, only.pos = TRUE,logfc.threshold = 1)
write.csv(ucec.markers2,file="markers.2.RNA.csv")
DefaultAssay(ucec)
DefaultAssay(ucec) = "RNA"

#创建marker集合
markers <- c("PTPRC", #immune
             "EPCAM", #epithelial
             "MME","PECAM1") #stromal

#画图
p <- FeaturePlot(ucec, features = markers, ncol = 2)
p
ggsave("1.pdf", p, width = 10, height = 10)

p <- DotPlot(ucec, features = markers) + RotatedAxis()
p
ggsave("2.pdf", p, width = 14, height = 7)

p <- VlnPlot(ucec, features = markers, stack = T, flip = T) + NoLegend()
p
ggsave("3.pdf", p, width = 14, height = 6)

pdf(file="4.pdf",width=7,height=6)
DimPlot(ucec, reduction = "umap", label = T)
dev.off()

#第一次
#immune epithelial stromal
ucec$celltype.1 <- recode(ucec@meta.data$seurat_clusters,
                          "0" = "immune",
                          "1" = "immune",
                          "2" = "stromal",
                          "3" = "stromal",
                          "4" = "stromal",
                          "5" = "epithelial",
                          "6" = "epithelial",
                          "7" = "immune",
                          "8" = "immune",
                          "9" = "epithelial",
                          "10" = "epithelial",
                          "11" = "stromal",
                          "12" = "epithelial",
                          "13" = "epithelial",
                          "14" = "epithelial",
                          "15" = "epithelial",
                          "16" = "stromal",
                          "17" = "epithelial",
                          "18" = "epithelial",
                          "19" = "epithelial",
                          "20" = "epithelial",
                          "21" = "stromal",
                          "22" = "stromal",
                          '23' = 'stromal',
                          '24' = 'epithelial',
                          '25' = 'immune',
                          '26' = 'stromal',
                          '27' = 'stromal',
                          '28' = 'stromal',
                          '29' = 'epithelial',
                          '30' = 'immune',
                          '31' = 'immune',
                          '32' = 'epithelial',
                          '33' = 'epithelial',
                          '34' = 'immune')

#查看各个细胞群的亚型
table(ucec@meta.data$celltype.1)

#换个配色
Biocols = c('#AB3282', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
            '#E95C59', '#E59CC4','#E5D2DD' , '#23452F', '#BD956A', '#8C549C', '#585658',
            '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
            '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
            '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
            '#968175')

#画图
p = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.1",cols = Biocols)
p
ggsave("5.pdf", p, width = 7, height = 6)







#第二次
#Marker基因标记
markers <- c("FUT4", #Epithelial_cells
             "CD68","CD163","CD14",'ADGRE1', #Macrophage
             "ITGB1","CD34",'PECAM1','VWF', #Endothelial_cells
             "ACTA2","TAGLN","MYH11", #Smooth_muscle_cells 
             "DCN","LUM","FGF7",'PDGFRA', 'MME',#Fibroblasts 
             "PROM1","ALDH1A1","LGR5", #Tissue_stem_cells
             "PDCD1","CTLA4","CD8A","PTPRC","CD4","BTLA","IL2RA","IL7R","CCR7","CD28","CD27","SLAMF1","DPP4","CD7","CD2","CD3G","CD3E","CD3D",'ENTPD1','ITGAE',#T_cells
             'THY1','ENG','NT5E', #MSC 
             'NCAM1','FCGR3A','NKG7','GNLY') #NK_cell

#画图
p <- FeaturePlot(ucec, features = markers, ncol = 5)
ggsave("6.pdf", p, width = 25, height = 45)
p <- DotPlot(ucec, features = markers) + RotatedAxis()
p
ggsave("annotation 7.pdf", p, width = 14, height = 7)
p <- VlnPlot(ucec, features = markers, stack = T, flip = T) + NoLegend()
p
ggsave("8.pdf", p, width = 14, height = 13)

#Fibroblast  Endothelia  Mast  Luminal Basal/intermediate  Monolytic  T
ucec$celltype.main <- recode(ucec@meta.data$seurat_clusters,
                             "0" = "T_cells",
                             "1" = "Macrophage",
                             "2" = "Endothelial_cells",
                             "3" = "Smooth_muscle_cells",
                             "4" = "Fibroblasts",
                             "5" = "Tissue_stem_cells",
                             "6" = "Epithelial_cells",
                             "7" = "T_cells",
                             "8" = "Macrophage",
                             "9" = "Epithelial_cells",
                             "10" = "Tissue_stem_cells",
                             "11" = "Smooth_muscle_cells",
                             "12" = "Tissue_stem_cells",
                             "13" = "Tissue_stem_cells",
                             "14" = "Tissue_stem_cells",
                             "15" = "Epithelial_cells",
                             "16" = "Smooth_muscle_cells",
                             "17" = "Epithelial_cells",
                             "18" = "Epithelial_cells",
                             "19" = "Epithelial_cells",
                             "20" = "Epithelial_cells",
                             "21" = "Smooth_muscle_cells",
                             "22" = "MSC",
                             '23' = 'Endothelial_cells',
                             '24' = 'Epithelial_cells',
                             '25' = 'NK_cell',
                             '26' = 'Smooth_muscle_cells',
                             '27' = 'Smooth_muscle_cells',
                             '28' = 'Smooth_muscle_cells',
                             '29' = 'Tissue_stem_cells',
                             '30' = 'Macrophage',
                             '31' = 'T_cells',
                             '32' = 'Epithelial_cells',
                             '33' = 'Epithelial_cells',
                             '34' = 'T_cells')

#查看各个细胞群的亚型
table(ucec@meta.data$celltype.main)
#画图
p = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.main")
ggsave("9.pdf", p, width = 7, height = 6)

#####插入 CXCL13 表达统计图 ================####
celltype_expr <- AverageExpression(ucec, assays = "RNA", features = "CXCL13", group.by = "celltype.main")

pdf("CXCL13_expression_by_celltype.pdf", width = 8, height = 8)

# 增加底部边距以容纳倾斜的标签
par(mar = c(10, 4, 4, 2) + 0.1)  # 下边距增加到10

# 创建条形图但不显示x轴标签
bp <- barplot(celltype_expr$RNA[1, ], 
              col = "steelblue", 
              main = "CXCL13 Expression by Cell Type",
              xlab = "",  # 清空x轴标签
              xaxt = "n") # 不显示x轴

# 添加倾斜45度的x轴标签
text(x = bp, 
     y = par("usr")[3] - 0.05 * (par("usr")[4] - par("usr")[3]), # 调整标签位置
     labels = names(celltype_expr$RNA[1, ]),
     srt = 45,    # 旋转45度
     adj = 1,     # 右对齐
     xpd = TRUE,  # 允许在绘图区域外绘制
     cex = 0.8)   # 调整字体大小

dev.off()


p1 = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.main")
p2 = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.1")
ggsave("10.pdf", p1|p2, width = 15, height = 6)


#第三次
#Fibroblast  Endothelia  Mast  Luminal Basal/intermediate  Monolytic  T
ucec$celltype.main <- recode(ucec@meta.data$seurat_clusters,
                             "0" = "T_cells",
                             "1" = "Macrophage",
                             "2" = "Endothelial_cells",
                             "3" = "Smooth_muscle_cells",
                             "4" = "Fibroblasts",
                             "5" = "Tissue_stem_cells",
                             "6" = "Epithelial_cells",
                             "7" = "T_cells",
                             "8" = "Macrophage",
                             "9" = "Epithelial_cells",
                             "10" = "Tissue_stem_cells",
                             "11" = "Smooth_muscle_cells",
                             "12" = "Tissue_stem_cells",
                             "13" = "Tissue_stem_cells",
                             "14" = "Tissue_stem_cells",
                             "15" = "Epithelial_cells",
                             "16" = "Smooth_muscle_cells",
                             "17" = "Epithelial_cells",
                             "18" = "Epithelial_cells",
                             "19" = "Epithelial_cells",
                             "20" = "Epithelial_cells",
                             "21" = "Smooth_muscle_cells",
                             "22" = "MSC",
                             '23' = 'Endothelial_cells',
                             '24' = 'Epithelial_cells',
                             '25' = 'NK_cell',
                             '26' = 'Smooth_muscle_cells',
                             '27' = 'Smooth_muscle_cells',
                             '28' = 'Smooth_muscle_cells',
                             '29' = 'Tissue_stem_cells',
                             '30' = 'Macrophage',
                             '31' = 'T_cells',
                             '32' = 'Epithelial_cells',
                             '33' = 'Epithelial_cells',
                             '34' = 'T_cells')

#查看各个细胞群的亚型
table(ucec@meta.data$celltype.main)

#画图
p = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.main")
p
ggsave("11.pdf", p, width = 7, height = 6)

p1 = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.main")
p2 = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.1")
p1|p2
ggsave("12.pdf", p1|p2, width = 15, height = 6)

p1 = DimPlot(ucec, reduction = "tsne", label = T,group.by = "celltype.main")
p2 = DimPlot(ucec, reduction = "tsne", label = T,group.by = "celltype.1")
p1|p2
p1 = DimPlot(ucec, reduction = "tsne", label = TRUE, group.by = "celltype.main", repel = TRUE)
p2 = DimPlot(ucec, reduction = "tsne", label = TRUE, group.by = "celltype.1", repel = TRUE)
p1 | p2
ggsave("13.pdf", p1|p2, width = 15, height = 6)






#第四次
#Marker基因标记
markers <- c("FUT4", #Epithelial_cells
             "CD68","CD163","CD14",'ADGRE1', #Macrophage
             "ITGB1","CD34",'PECAM1','VWF', #Endothelial_cells
             "ACTA2","TAGLN","MYH11", #Smooth_muscle_cells 
             "DCN","LUM","FGF7",'PDGFRA', 'MME',#Fibroblasts 
             "PROM1","ALDH1A1","LGR5", #Tissue_stem_cells
             "PDCD1","CTLA4","CD8A","PTPRC","CD4","BTLA","IL2RA","IL7R","CCR7","CD28","CD27","SLAMF1","DPP4","CD7","CD2","CD3G","CD3E","CD3D",'ENTPD1','ITGAE',#T_cells
             'THY1','ENG','NT5E', #MSC 
             'NCAM1','FCGR3A','NKG7','GNLY') #NK_cell


#画图
p <- DotPlot(ucec, features = markers,group.by = "celltype.main") + RotatedAxis()
ggsave("14.pdf", p, width = 14, height = 6)

#查看各个细胞群的亚型
table(ucec@meta.data$celltype.main)









#分组
table(Idents(ucec))
Idents(ucec) = "celltype.main"
table(Idents(ucec))

#矩阵
DefaultAssay(ucec)

#寻找每个细胞类型的markers
ucec.markers3 <- FindAllMarkers(ucec, only.pos = TRUE,logfc.threshold = 1,min.pct = 0.3)
write.csv(ucec.markers3,file="markers.celltype.RNA.csv")

#提取前10的marker
top10 <- ucec.markers3 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#对这些基因进行scaledata
markers = as.data.frame(top10[,"gene"])
ucec <- ScaleData(ucec, features = as.character(unique(markers$gene)))

#画图
p = DoHeatmap(ucec,
              features = as.character(unique(markers$gene)),
              group.by = "celltype.main")
p
ggsave("15.pdf", p, width = 14, height = 12)

#每个细胞亚群抽最少的细胞类型
allCells = names(Idents(ucec))
allType = levels(Idents(ucec))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(ucec)== x ]
  cg=sample(cgCells,min(table(ucec@meta.data$celltype.main)))
  cg
}))

#提取
cg_sce = ucec[, allCells %in% choose_Cells]
table(Idents(cg_sce))

#画图
p = DoHeatmap(cg_sce,
              features = as.character(unique(markers$gene)),
              group.by = "celltype.main")
ggsave("16.pdf", p, width = 14, height = 12)




#堆叠柱状图
cell.prop<-as.data.frame(prop.table(table(ucec@meta.data$celltype.main, ucec@meta.data$orig.ident)))
colnames(cell.prop)<-c("cluster","group","proportion")

p = ggplot(cell.prop,aes(group,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))
p <- ggplot(cell.prop, aes(group, proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("") +
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.5, 'cm'),
    axis.text.x = element_text(angle = 45, hjust = 1)  # 倾斜45度并右对齐
  ) +
  guides(fill = guide_legend(title = NULL))
ggsave("17.pdf", p, width = 10, height = 9)

#保存
saveRDS(ucec,"7ucec_celltype.rds")
#####插入：CXCL13 高/低表达组的细胞组成差异分析####
# 根据CXCL13中位数分组
ucec$cxcl13_group <- ifelse(
  GetAssayData(ucec, assay = "RNA")["CXCL13", ] > 
    median(GetAssayData(ucec, assay = "RNA")["CXCL13", ], na.rm = TRUE),
  "High", "Low"
)

# 计算各细胞类型在两组中的比例
cell_prop <- prop.table(table(ucec$celltype.main, ucec$cxcl13_group), margin = 2)

# 绘制堆叠柱状图
pdf("CXCL13_high_low_cell_composition.pdf", width = 8, height = 8)
barplot(
  cell_prop, 
  col = brewer.pal(nrow(cell_prop), "Set2"), 
  legend.text = rownames(cell_prop),
  xlab = "CXCL13 Expression Group",
  ylab = "Proportion",
  main = "Cell Type Composition by CXCL13 Expression"
)
dev.off()

# 计算各细胞类型在两组中的比例
cell_prop <- prop.table(table(ucec$celltype.main, ucec$cxcl13_group), margin = 2)

# 绘制堆叠柱状图
pdf("CXCL13_high_low_cell_composition.pdf", width = 12, height = 8)  # 进一步增加宽度

# 设置更大的右边距以容纳图例
par(mar = c(5, 4, 4, 12) + 0.1, xpd = TRUE)

# 绘制柱状图（不显示图例）
barplot(
  cell_prop, 
  col = brewer.pal(nrow(cell_prop), "Set2"), 
  xlab = "CXCL13 Expression Group",
  ylab = "Proportion",
  main = "Cell Type Composition by CXCL13 Expression",
  legend = FALSE  # 不显示默认图例
)

# 在图表右侧外部添加图例，调整位置
legend("topright", 
       inset = c(-0.18, 0),  # 减小负值，使图例更靠近图表
       legend = rownames(cell_prop), 
       fill = brewer.pal(nrow(cell_prop), "Set2"),
       title = "Cell Types",
       cex = 0.8)  # 减小字体大小

dev.off()




#####在这里插入CXCL13与免疫细胞浸润相关性分析的代码####
# 首先需要计算免疫细胞比例矩阵（假设已通过某种方法获得）

## 直接从 Seurat 对象创建免疫细胞比例矩阵
cell_types <- unique(ucec$celltype.main)
cat("细胞类型:", cell_types, "\n")

# 创建二进制矩阵表示每个细胞属于哪种类型
immune_cell_matrix <- matrix(
  0, 
  nrow = ncol(ucec), 
  ncol = length(cell_types),
  dimnames = list(colnames(ucec), cell_types)
)

# 填充矩阵
for(i in 1:length(cell_types)) {
  cell_type <- cell_types[i]
  cells_of_type <- which(ucec$celltype.main == cell_type)
  immune_cell_matrix[cells_of_type, i] <- 1
}

# 提取 CXCL13 表达数据
cxcl13_expression <- as.numeric(GetAssayData(ucec, assay = "RNA", slot = "data")["CXCL13", ])

# 检查维度是否匹配
cat("CXCL13 长度:", length(cxcl13_expression), "\n")
cat("免疫细胞矩阵维度:", dim(immune_cell_matrix), "\n")

# 现在计算相关性
cor_data <- cor(cxcl13_expression, immune_cell_matrix, method = "spearman")
print(cor_data)


# 假设您有 p 值数据（这里需要您实际计算 p 值）
# 以下仅为示例，您需要替换为实际的 p 值计算
p_values <- c(0.001, 0.05, 0.1, 0.01, 0.0001, 0.5, 0.8, 0.6, 0.7)  # 示例 p 值

# 创建数据框
cor_df <- data.frame(
  CellType = colnames(immune_cell_matrix),
  Correlation = as.numeric(cor_data),
  PValue = p_values,
  FDR = p.adjust(p_values, method = "fdr")
)

# 添加显著性标记
cor_df$Significance <- ifelse(cor_df$FDR < 0.001, "***",
                              ifelse(cor_df$FDR < 0.01, "**",
                                     ifelse(cor_df$FDR < 0.05, "*", "")))

# 绘制带有显著性标记的条形图
ggplot(cor_df, aes(x = reorder(CellType, Correlation), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  coord_flip() +
  theme_minimal() +
  geom_text(aes(label = Significance), 
            vjust = 0.5, hjust = ifelse(cor_df$Correlation >= 0, -0.2, 1.2), 
            size = 5, color = "black") +
  geom_text(aes(label = round(Correlation, 3)), 
            vjust = 0.5, hjust = ifelse(cor_df$Correlation >= 0, 1.2, -0.2), 
            size = 3.5, color = "black") +
  labs(title = "Correlation between CXCL13 Expression and Cell Types",
       x = "Cell Type", y = "Spearman Correlation Coefficient") +
  theme(plot.title = element_text(hjust = 0.5))
# 保存为 PDF
pdf("cxcl13_celltype_correlation_visualization.pdf", width = 10, height = 8)

# 选择您喜欢的可视化方法，例如：
ggplot(cor_df, aes(x = reorder(CellType, Correlation), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  coord_flip() +
  theme_minimal() +
  geom_text(aes(label = round(Correlation, 3)), 
            hjust = ifelse(cor_df$Correlation >= 0, -0.1, 1.1), 
            size = 3.5) +
  labs(title = "Correlation between CXCL13 Expression and Cell Types",
       x = "Cell Type", y = "Spearman Correlation Coefficient") +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()

# 计算每个细胞类型的相关性和p值
cell_types <- colnames(immune_cell_matrix)
cor_results <- numeric(length(cell_types))
p_values <- numeric(length(cell_types))

for(i in 1:length(cell_types)) {
  cor_test <- cor.test(cxcl13_expression, immune_cell_matrix[, i], method = "spearman")
  cor_results[i] <- cor_test$estimate
  p_values[i] <- cor_test$p.value
}

# 创建结果数据框
results_df <- data.frame(
  CellType = cell_types,
  Correlation = cor_results,
  PValue = p_values,
  FDR = p.adjust(p_values, method = "fdr")
)

# 查看结果
print(results_df)

#提取细胞类型进行细分亚群定义
ucec.T = ucec[, Idents(ucec) %in% c("T_cells")]
table(Idents(ucec.T))   # 应该只剩 T 这一类
dim(ucec.T)             # 行是基因，列是 T 细胞数
# 3. 后续常用：重新降维 + 亚群再聚类
ucec.T <- NormalizeData(ucec.T)
ucec.T <- FindVariableFeatures(ucec.T)
ucec.T <- ScaleData(ucec.T)
ucec.T <- RunPCA(ucec.T, npcs = 30)
ucec.T <- FindNeighbors(ucec.T, dims = 1:20)
ucec.T <- FindClusters(ucec.T, resolution = 0.6)
ucec.T <- RunUMAP(ucec.T, dims = 1:20)

# 4. 查看 CXCL13 在 T 细胞亚群中的表达
FeaturePlot(ucec.T, features = "CXCL13", label = TRUE)
VlnPlot(ucec.T, features = "CXCL13", pt.size = 0)
# 绘制T细胞亚群的UMAP图
DimPlot(ucec.T, label = TRUE, repel = TRUE) + NoLegend()
# 使用T细胞特异性标记来鉴定亚群
T_subset_markers <- c("CD4", "CD8A", "CD8B", "FOXP3", 
                      "CCR7", "SELL", # naive/memory
                      "GZMB", "GZMK", "GZMA", # cytotoxic
                      "PDCD1", "CTLA4", "HAVCR2", "LAG3", # exhaustion
                      "CXCL13", "BCL6", # Tfh
                      "IL2RA", # Treg
                      "NKG7", "GNLY") # NK-like
# 绘制气泡图
DotPlot(ucec.T, features = T_subset_markers) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 寻找各亚群的差异基因
markers <- FindAllMarkers(ucec.T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 保存整个ucec.T对象
saveRDS(ucec.T, file = "ucec_T_cells_subset.rds")

#####2. T细胞分群####
# 1. 载入R包
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(plyr)
library(scales)

# 设置工作目录
setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

# 2. 读取T细胞子集数据
ucec.T <- readRDS("ucec_T_cells_subset.rds")

# 3. 重新分析流程（确保一致性）
DefaultAssay(ucec.T) <- "RNA"
ucec.T <- NormalizeData(ucec.T) %>% 
  FindVariableFeatures(nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = FALSE)

# 选择最佳维度（根据肘部法则）
ElbowPlot(ucec.T, ndims = 30)
ggsave("T_cell_elbow_plot.pdf", width = 6, height = 4)

dims.use <- 1:20
ucec.T <- FindNeighbors(ucec.T, dims = dims.use) %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dims = dims.use)


# 4. 定义T细胞标记基因
T.markers <- c(
  # CD4+ T细胞
  "CD4", "CD3D", 
  # CD8+ T细胞  
  "CD8A", "CD8B",
  # Naive/Central Memory
  "CCR7", "SELL","CD28",
  # 调节性T细胞
  "FOXP3", "IL2RA", "CTLA4",'ENTPD1',#后两个有抑制标志作用
  # 耗竭T细胞
  "PDCD1", "HAVCR2", "LAG3", "TIGIT","BTLA",
  # 细胞毒性/NK样
  "GZMB", "GZMA", "PRF1","GZMK","NKG7", "GNLY",
  # Tfh
  "CXCL13", "BCL6","CD27","SLAMF1",
  # 组织驻留记忆T细胞 (Trm)
  'ITGAE'
)

# 5. 查看初始聚类标记表达
dot_initial <- DotPlot(ucec.T, features = T.markers, group.by = "seurat_clusters") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  labs(title = "T Cell Markers by Initial Clusters")
ggsave("T_dotplot_raw_cluster.pdf", dot_initial, width = 12, height = 8)

# 6. 基于标记表达进行亚群注释（根据您的DotPlot结果调整）
new.ids <- c(
  "0" = "CD8+_T",
  "1" = "Naive_T/Tcm", 
  "2" = "Naive_T/Tcm",
  "3" = "CD8+_T",
  "4" = "Treg",
  "5" = "Treg", 
  "6" = "CD8+_T" ,   
  "7" = "Treg",
  "8" = "Tfh",
  "9" = "CD8+_T",
  "10" = "CD8+_T",
  "11" = "exhausted_CD8+_T",
  "12" = "NK-like_T",
  "13" = "Trm",
  "14" = "Trm",
  "15" = "Treg",
  "16" = "exhausted_CD8+_T"
)

# 应用注释
ucec.T$subT <- mapvalues(Idents(ucec.T), 
                         from = names(new.ids), 
                         to = new.ids)
Idents(ucec.T) <- "subT"

# 7. 高质量可视化

# 7.1 UMAP图
umap_plot <- DimPlot(ucec.T, reduction = "umap", label = TRUE, 
                     repel = TRUE, label.size = 4, pt.size = 0.5) +
  theme_minimal() +
  ggtitle("T Cell Subtypes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
ggsave("T_subtypes_umap.pdf", umap_plot, width = 8, height = 7)

# 7.2 气泡图（按最终亚群）
dot_final <- DotPlot(ucec.T, features = T.markers, group.by = "subT") +
  RotatedAxis() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  labs(title = "T Cell Subtype Markers")
ggsave("T_subtypes_dotplot.pdf", dot_final, width = 14, height = 8)

# 7.3 关键标记基因FeaturePlot
feature_genes <- c("CXCL13")
feature_plot <- FeaturePlot(ucec.T, features = feature_genes, 
                            ncol = 3, order = TRUE, pt.size = 0.3) &
  theme_minimal() &
  theme(legend.position = "bottom",
        plot.title = element_text(size = 12))
ggsave("T_markers_featureplot.pdf", feature_plot, width = 16, height = 8)

# 7.4 小提琴图
vln_plot <- VlnPlot(ucec.T, features = feature_genes, 
                    pt.size = 0, ncol = 3, fill.by = "ident") &
  theme_minimal() &
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave("T_markers_vlnplot.pdf", vln_plot, width = 16, height = 8)

# 7.5 亚群比例统计和可视化
cell_counts <- as.data.frame(table(ucec.T$subT))
colnames(cell_counts) <- c("Subtype", "Count")
cell_counts$Percentage <- round(cell_counts$Count / sum(cell_counts$Count) * 100, 1)

# 柱状图
bar_plot <- ggplot(cell_counts, aes(x = reorder(Subtype, -Count), y = Count, fill = Subtype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
            vjust = -0.3, size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "T Cell Subtype", y = "Cell Count", 
       title = "T Cell Subtype Distribution")
bar_plot <- ggplot(cell_counts, aes(x = reorder(Subtype, -Count), y = Count, fill = Subtype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
            vjust = -0.3, size = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(margin = margin(b = 20))  # 增加标题底部边距
  ) +
  labs(x = "T Cell Subtype", y = "Cell Count", 
       title = "T Cell Subtype Distribution") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  # 扩展y轴上限，为注释留出空间
bar_plot
ggsave("T_subtype_distribution.pdf", bar_plot, width = 10, height = 8)

# 8. 保存带注释的对象
saveRDS(ucec.T, file = "ucec_T_cells_subtyped.rds")

# 9. 综合报告图
composite_plot <- (umap_plot | dot_final) / 
  (feature_plot | vln_plot) +
  plot_annotation(tag_levels = 'A', 
                  title = "T Cell Subtype Analysis Summary") +
  plot_layout(heights = c(1, 2))
composite_plot <- (umap_plot | dot_final) / 
  (feature_plot | vln_plot) +
  plot_annotation(
    tag_levels = 'A', 
    title = "T Cell Subtype Analysis Summary"
  ) +
  plot_layout(heights = c(1, 2))
ggsave("T_cell_analysis_composite.pdf", composite_plot, width = 18, height = 16)

# 10. 输出统计信息
write.csv(cell_counts, "T_cell_subtype_statistics.csv", row.names = FALSE)

cat("分析完成！生成的文件：\n")
cat("- T_subtypes_umap.pdf: T细胞亚群UMAP图\n")
cat("- T_subtypes_dotplot.pdf: 标记基因气泡图\n") 
cat("- T_markers_featureplot.pdf: 关键基因表达分布图\n")
cat("- T_markers_vlnplot.pdf: 关键基因小提琴图\n")
cat("- T_subtype_distribution.pdf: 亚群分布柱状图\n")
cat("- T_cell_analysis_composite.pdf: 综合分析报告图\n")
cat("- ucec_T_cells_subtyped.rds: 带注释的Seurat对象\n")
cat("- T_cell_subtype_statistics.csv: 亚群统计信息\n")

####十一、monocle####
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ROGUE)
library(clustree)
library(harmony)
library(SingleR)
library(dplyr)
#library(BiocManager)
#BiocManager::install("monocle")
library(monocle)
library(tidyverse)

#目录
setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

#读入
ucec.T <- readRDS("ucec_T_cells_subtyped.rds")

#画图
DimPlot(ucec.T, group.by = "subT", label = T) + DimPlot(ucec.T, group.by = "seurat_clusters", label = T)

#命名
ucec.T$celltype <- ucec.T$subT
colnames(ucec.T@meta.data)

#获取表达矩阵
data <- GetAssayData(ucec.T, assay = "RNA", slot = "counts")
#获取细胞注释信息，表型信息
pd <- new('AnnotatedDataFrame', data = ucec.T@meta.data[,c(1,7,37,39)])
#获取基因名
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#创建一个新的细胞数据集对象
dim(ucec)
mycds <- newCellDataSet(as.matrix(data),
                        #若数据太大，转化为稀疏矩阵
                        #as(as.matrix(data),"sparseMatrix"),
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())

#数据预处理，估计size factor和离散度，类似归一化，标准化
#运行时间较长，可修改核心数
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=8)

### 选择排序基因，2000左右
disp_table <- dispersionTable(mycds)
order.genes <- subset(disp_table, mean_expression >= 0.005 & dispersion_empirical >= 
                        1 * dispersion_fit) %>% pull(gene_id) %>% as.character()

#标记这些基因
mycds <- setOrderingFilter(mycds, order.genes)

#设置排序基因
plot_ordering_genes(mycds)
p <- plot_ordering_genes(mycds)
ggsave("0OrderGenes.pdf", p, width = 8, height = 6)

### 降维排序
#residualModelFormulaStr 
#减去“无趣的”变异源的影响，以减少它们对集群的影响。
mycds <- reduceDimension(mycds, max_components = 2, reduction_method = 'DDRTree', 
                         residualModelFormulaStr = "~orig.ident")


mycds <- orderCells(mycds)


#若出现报错
#Error :
#  !nei() was deprecated in igraph 2.1.0 and is now defunct.
# please use.nei() instead.
#Run rlang::last trace() to see where the error occurred.

#需要将igraph降低至1.5.1版本
#首先卸载这个包  
#remove.packages("igraph")
#然后安装指定版本 
#install.packages("igraph")
#packageurl = 'https://cran.r-project.org/src/contrib/Archive/igraph/igraph_1.5.1.tar.gz'
#install.packages(packageurl, repos = NULL, type = 'source')
#如果仍然无法安装，看2.R包安装视频，本地安装。

#报错
#Error in if (class(projection) != "matrix") projection <- as.matrix(projection) :
#运行
#trace('project2MST', edit = T, where = asNamespace("monocle"))
#找到
#if (elass(projection) != "matrix")
#删除并保存

# 可以人工设置起始点，根
#mycds <- orderCells(mycds,root_state = 5)

# 结果可视化
# naive (Tn), central memory (Tcm), effector memory (Tem)
# regulatory (tregs)，T helper (Th)

# State
p1 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("1Trajectory_State.pdf", plot = p1, width = 10, height = 6.5)
# Pseudotime
p2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("2Trajectory_Pseudotime.pdf", plot = p2, width = 10, height = 6.5)
# Celltype
p3 <- plot_cell_trajectory(mycds, color_by = "celltype")
ggsave("3Trajectory_Celltype2.pdf", plot = p3, width = 10, height = 6.5)
# orig.ident
p4 <- plot_cell_trajectory(mycds, color_by = "orig.ident")
ggsave("4Trajectory_Sample.pdf", plot = p4, width = 10, height = 6.5)
# 树形图
p5 <- plot_complex_cell_trajectory(mycds, x = 1, y = 2,
                                   color_by = "celltype")
ggsave("5Trajectory_dendrogram.pdf", plot = p5, width = 10, height = 6.5)
# 细胞密度图
p6 <- ggplot(pData(mycds),aes(Pseudotime,colour = celltype,fill = celltype)) +
  geom_density(bw = 0.5, size = 1, alpha = 0.5)+theme_classic()
ggsave("6Trajectory_Density.pdf", plot = p6, width = 10, height = 6.5)
# 指定基因的表达变化
genes = c(order.genes)[1:4]
p1 = plot_genes_in_pseudotime(mycds[genes],color_by="State")
p2 = plot_genes_in_pseudotime(mycds[genes],color_by="celltype")
p3 = plot_genes_in_pseudotime(mycds[genes],color_by="Pseudotime")
ggsave("7Trajectory_Pseudotime.pdf", plot = p1|p2|p3, width = 10, height = 6.5)

p1 = plot_genes_jitter(mycds[genes],grouping="State",color_by="State")
p2 = plot_genes_violin(mycds[genes],grouping="State",color_by="State")
p3 = plot_genes_in_pseudotime(mycds[genes],color_by="State")
ggsave("8Trajectory_jitter.pdf", plot = p1|p2|p3, width = 10, height = 6.5)

pData(mycds)$SAT1 = log2(exprs(mycds)["SAT1",] +1)
p1 = plot_cell_trajectory(mycds,color_by = "SAT1") + 
  scale_color_continuous(type = "viridis")

pData(mycds)$KLRB1 = log2(exprs(mycds)["KLRB1",] +1)
p2 = plot_cell_trajectory(mycds,color_by = "KLRB1") + 
  scale_color_continuous(type = "viridis")

pData(mycds)$HSPA8 = log2(exprs(mycds)["HSPA8",] +1)
p3 = plot_cell_trajectory(mycds,color_by = "HSPA8") + 
  scale_color_continuous(type = "viridis")

pData(mycds)$G3BP2 = log2(exprs(mycds)["G3BP2",] +1)
p4 = plot_cell_trajectory(mycds,color_by = "G3BP2") + 
  scale_color_continuous(type = "viridis")

ggsave("9Trajectory_Expression.pdf",plot = p1|p2|p3|p4,width=14,height=6.5)

#保存至seurat对象
pdata <- Biobase::pData(mycds)
ucec <- AddMetaData(ucec, metadata = pdata[,c("Pseudotime","State")])
saveRDS(ucec, file = "10ucec.pseudotime.rds")

#寻找拟时差异基因，通过monocle方法
Time_diff <- differentialGeneTest(mycds, cores = 10,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff, "11Time_diff_all.csv", row.names = F)
#画图，按qval排序，选取前100个
Time_genes <- Time_diff[order(Time_diff$qval), "gene_short_name"][1:100]
#num_clusters按行聚类，聚为多少类
p = plot_pseudotime_heatmap(mycds[Time_genes,], num_clusters=3, 
                            show_rownames=T, return_heatmap=T)
ggsave("12Time_heatmap.pdf", p, width = 5, height = 10)
#保存
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(Time_diff_sig, "13Time_diff_sig.csv", row.names = F)

# 单细胞轨迹的“分支”分析
# 寻找和分叉点相关的基因
# BEAM分析,用于寻找以依赖于分支的方式调控的基因。
beam_res <- BEAM(mycds, branch_point = 1, cores = 10,
                 progenitor_method = "duplicate")
write.csv(beam_res, "14BEAM_all.csv", row.names = F)
#前100个
BEAM_genes <- beam_res[order(beam_res$qval), "gene_short_name"][1:100]
p <- plot_genes_branched_heatmap(
  mycds[BEAM_genes,],  branch_point = 1, num_clusters = 3, show_rownames = T, return_heatmap = T)
ggsave("15BEAM_heatmap.pdf", p$ph_res, width = 6.5, height = 10)
#保存
hp.genes <- p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig <- beam_res[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(BEAM_sig, "16BEAM_sig.csv", row.names = F)




####CXCL13鉴定####
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ROGUE)
library(clustree)
library(harmony)
library(SingleR)
library(dplyr)

#目录
setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

#读入
ucec.T <- readRDS("ucec_T_cells_subtyped.rds")
# 3. 后续常用：重新降维 + 亚群再聚类
ucec.T <- NormalizeData(ucec.T)
ucec.T <- FindVariableFeatures(ucec.T)
ucec.T <- ScaleData(ucec.T)
ucec.T <- RunPCA(ucec.T, npcs = 30)
ucec.T <- FindNeighbors(ucec.T, dims = 1:20)
ucec.T <- FindClusters(ucec.T, resolution = 0.6)
ucec.T <- RunUMAP(ucec.T, dims = 1:20)

# 4. 查看 CXCL13 在 T 细胞亚群中的表达
FeaturePlot(ucec.T, features = "CXCL13", label = TRUE)
VlnPlot(ucec.T, features = "CXCL13", pt.size = 0)
# 绘制T细胞亚群的UMAP图
DimPlot(ucec.T, label = TRUE, repel = TRUE) + NoLegend()
# 使用T细胞特异性标记来鉴定亚群
T_subset_markers <- c("CD4", "CD8A", "CD8B", "FOXP3", 
                      "CCR7", "SELL", # naive/memory
                      "GZMB", "GZMK", "GZMA", # cytotoxic
                      "PDCD1", "CTLA4", "HAVCR2", "LAG3", # exhaustion
                      "CXCL13", "BCL6", # Tfh
                      "IL2RA", # Treg
                      "NKG7", "GNLY") # NK-like
# 绘制气泡图
DotPlot(ucec.T, features = T_subset_markers) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 寻找各亚群的差异基因
markers <- FindAllMarkers(ucec.T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




####CXCL13如何影响TME细胞互作网络####

# 加载必要的包
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)

# 1. 创建CellChat对象（假设已有seurat_obj，且含有细胞注释信息'celltype'）
data.input <- GetAssayData(ucec.T, assay = "RNA", slot = "data") # 获取标准化表达矩阵
meta.data <- ucec.T@meta.data # 获取细胞元数据

# 创建CellChat对象
cellchat <- createCellChat(object = data.input, meta = meta.data, group.by = "subT")

# 2. 使用配体-受体数据库
CellChatDB <- CellChatDB.human # 使用人类数据库
cellchat@DB <- CellChatDB

#提取细胞通讯信号基因
cellchat <- subsetData(cellchat)
# --- 关键修改部分：跳过过表达分析，直接进行通信计算 --- #
# 我们不运行 identifyOverExpressedGenes 和 identifyOverExpressedInteractions
# 而是直接使用 projectData 将我们的数据投射到PPI网络上，为计算通信概率做准备
cellchat <- projectData(cellchat, PPI.human) # 对于人和小鼠，可以使用内置的PPI数据

cellchat <- createCellChat(object = data.input, meta = meta.data, group.by = "subT") # 或 "Mouse"
cellchat <- computeCommunProbPathway(cellchat)
错误于apply(prob, c(1, 2), by, group, sum): dim(X)的值必需是正数
# 3. 计算通信概率
# 使用 raw.use = FALSE 使用标准化后的数据
# 使用 population.size = TRUE 考虑细胞群体大小的影响
cellchat <- computeCommunProb(cellchat)

# 可选：如果结果中仍有太多微弱信号，可以按需过滤
# cellchat <- filterCommunication(cellchat, min.cells = 10)

# 4. 聚合计算通路水平的通信概率
cellchat <- computeCommunProbPathway(cellchat)

# 5. 聚合整个细胞通信网络
cellchat <- aggregateNet(cellchat)

# 6. 可视化与分析：重点关注CXCL13
# 方法A：直接提取CXCL13相关的所有相互作用
# 找到所有以CXCL13为配体的受体
cxcl13_ligand_pairs <- subset(CellChatDB$interaction, ligand == "CXCL13")
# 找到所有以CXCL13为受体的配体
cxcl13_receptor_pairs <- subset(CellChatDB$interaction, receptor == "CXCL13")

print("CXCL13作为配体的相互作用对：")
print(cxcl13_ligand_pairs[, c("ligand", "receptor", "pathway_name")])
print("CXCL13作为受体的相互作用对：")
print(cxcl13_receptor_pairs[, c("ligand", "receptor", "pathway_name")])

# 方法B：手动指定我们感兴趣的配体-受体对（例如CXCL13-CXCR5）进行可视化
pair_to_show <- c("CXCL13", "CXCR5")

# 绘制特定配体-受体对的通信网络图
netVisual_individual(
  cellchat,
  signaling = "CXCL", # 这里填写该对所属的通路名，在cxcl13_ligand_pairs中可查到，通常是"CXC"
  pairLR.use = data.frame(ligand = pair_to_show[1], receptor = pair_to_show[2]),
  layout = "circle",
  vertex.weight = NULL # 顶点大小不按细胞数量加权
)

# 绘制气泡图，显示哪些细胞群是CXCL13的发送者，哪些是CXCR5的接收者
netVisual_bubble(
  cellchat,
  sources.use = NULL, # 显示所有发送者
  targets.use = NULL, # 显示所有接收者
  signaling = "CXCL", # 指定通路
  pairLR.use = data.frame(ligand = pair_to_show[1], receptor = pair_to_show[2]),
  remove.isolate = FALSE # 即使不相互作用也显示
)

# 方法C：提取并查看CXCL13通信强度的数值结果
df.net <- subsetCommunication(cellchat) # 获取所有相互作用的详细数据框
# 筛选出涉及CXCL13的通信
cxcl13_comm <- df.net[df.net$ligand == "CXCL13" | df.net$receptor == "CXCL13", ]
head(cxcl13_comm[order(-cxcl13_comm$prob), ]) # 按通信概率从高到低排序

# 7. 保存结果
saveRDS(cellchat, file = "cellchat_analysis_with_CXCL13.rds")


