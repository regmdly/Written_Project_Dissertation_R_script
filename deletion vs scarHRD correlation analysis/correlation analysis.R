#load data:scar HRD score for 69 samples, copy number profile for 80 samples
scarHRD = read.csv('vortex_Dx_scarHRD_results.CSV', header = T)
CN = read.csv('copy_number_change.csv', header = T)

library(tidyverse)
score = select(scarHRD, HRD.sum, sample)
CN = select(CN, sample, gene, CN_status)

##############correlation analysis between scar HRD score and level of deletion 
#categorize LOH and homozygous deletion as 'deletion'
by_CN = group_by(CN, sample)
Deletion = count(by_CN, CN_status == 'loss of heterozygosity'|CN_status =="homozygous deletion")
colnames(Deletion) <- c("sample", "deletion","level_dele")

#identify samples that do not have deletion
t = filter(Deletion, deletion == TRUE)
f = filter(Deletion, deletion == FALSE)
f = f[! f$sample %in% t$sample,]
f[,3]= 0

t = select(t, sample, level_dele)
f = select(f, sample, level_dele)
level_dele = rbind(t,f)
#####another way to do this, assign gain or diploid as 0, assign LOH, homozygous deletion as 1. then sum up


#a new data frame combining the scarHRD score and level of deletion 
Combined = left_join(level_dele, score, by= 'sample')
Combined = drop_na(Combined)

# do correlation analysis between level of deletion versus scarHRD score############

# Install
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library("ggpubr")

#check normality by ploting or shapiro test
ggplot(data = Combined)+
  geom_histogram(mapping = aes(x = level_dele), binwidth = 1)

ggplot(data = Combined)+
  geom_histogram(mapping = aes(x = HRD.sum), binwidth = 2)

shapiro.test(Combined$level_dele)
shapiro.test(Combined$HRD.sum)

#Kendall rank correlation test and Spearman rank correlation coefficient
cor.test(Combined$level_dele, Combined$HRD.sum,  method="kendall")
cor.test(Combined$level_dele, Combined$HRD.sum,  method="spearman")

#Visualize data using scatter plots
ggscatter(Combined, x = "level_dele", y = "HRD.sum", 
          add = "reg.line", conf.int = TRUE,
          xlab = "level of deletion", ylab = "scarHRD score"
          )+
  stat_cor(method = "kendall", label.x = 11, label.y = 79)

#check how many samples involved in each data frame
n_distinct(CN$sample)
n_distinct(level_dele$sample)
n_distinct(Combined$sample)

#########################################################################
#correlation analysis between scar HRD score and level of deletion at individual gene level
#determining HRD positive/negative and presence/absence of deletion 
Combined2 = left_join(CN, score)
Combined2 = drop_na(Combined2)
Combined2 = mutate(Combined2, HRD.binary=NA)

for (i in 1:nrow(Combined2)) {
  if(Combined2$HRD.sum[i]>=42){
    Combined2$HRD.binary[i]='positive'
  } else if (Combined2$HRD.sum[i]<42){
    Combined2$HRD.binary[i]='negative'
  }
}

for (i in 1:nrow(Combined2)) {
  if(Combined2$CN_status[i]=='loss of heterozygosity'|Combined2$CN_status[i]=='homozygous deletion'){
    Combined2$CN_status[i]='deletion'
  } else if (Combined2$CN_status[i]=='gain'|Combined2$CN_status[i]=='diploid'){
    Combined2$CN_status[i]='no deletion'
  }
}

#conduct fisher exact test for each gene and output p.values
result= matrix(c(1,1,1),ncol = 3,nrow = 16)
colnames(result) = c("gene", "p.value_42",'estimates_odds_ratio_42')
k=1
for (i in unique(Combined2$gene)) {
  a=filter(Combined2, gene == i)
  b=xtabs(~CN_status+HRD.binary, data = a)
  c=fisher.test(b)
  result[k,3] = round(c$estimate, digits = 4)
  result[k,2] = round(c$p.value, digits = 4)
  result[k,1] = i
  k =k+1
}

#determining statistical significant according to significance level 0.05
result=as.data.frame(result)
result= mutate(result, p.value_binary=NA)
for (i in 1:nrow(result)) {
  if(result$p.value_42[i]>=0.05){
    result$p.value_binary[i]='not statistically significant'
  } else if (result$p.value_42[i]<0.05){
    result$p.value_binary[i]='statistically significant'
  }
}


#fisher exact test for ATM gene
ATM = filter(Combined2, gene =='ATM')
ATM = xtabs(~CN_status+HRD.binary, data = ATM)
fisher.test(ATM)

#########################################fisher exact test for each gene, using 63 as cutoff value for HRD positive/negative
Combined3 = mutate(Combined2, HRD.binary=NA)

for (i in 1:nrow(Combined3)) {
  if(Combined3$HRD.sum[i]>=63){
    Combined3$HRD.binary[i]='positive'
  } else if (Combined3$HRD.sum[i]<63){
    Combined3$HRD.binary[i]='negative'
  }
}

for (i in 1:nrow(Combined3)) {
  if(Combined3$CN_status[i]=='loss of heterozygosity'|Combined3$CN_status[i]=='homozygous deletion'){
    Combined3$CN_status[i]='deletion'
  } else if (Combined3$CN_status[i]=='gain'|Combined3$CN_status[i]=='diploid'){
    Combined3$CN_status[i]='no deletion'
  }
}

#conduct fisher exact test for each gene and output p.values
result2= matrix(c(1,1,1),ncol = 3,nrow = 16)
colnames(result2) = c("gene", "p.value_63", 'estimates_odds_ratio')
k=1
for (i in unique(Combined3$gene)) {
  a=filter(Combined3, gene == i)
  b=xtabs(~CN_status+HRD.binary, data = a)
  c=fisher.test(b)
  result2[k,3] = round(c$estimate, digits = 4)
  result2[k,2] = round(c$p.value, digits = 4)
  result2[k,1] = i
  k =k+1
}

#determining statistical significance according to significance level 0.05
result2=as.data.frame(result2)
result2= mutate(result2, p.value_binary=NA)
for (i in 1:nrow(result2)) {
  if(result2$p.value_63[i]>=0.05){
    result2$p.value_binary[i]='not statistically significant'
  } else if (result2$p.value_63[i]<0.05){
    result2$p.value_binary[i]='statistically significant'
  }
}


xtabs(~CN_status+HRD.binary+gene, data = Combined2)
xtabs(~CN_status+HRD.binary+gene, data = Combined3)



