library(tidyverse)
#load data:scar HRD score for 80 samples, copy number profile for 80 samples
scarHRD_80 = read.csv('scarHRD_80samples.csv', header = T)
CN = read.csv('copy_number_change.csv', header = T)

score_80 = select(scarHRD_80, HRD.sum, SampleID)
colnames(score_80)[2]='sample'
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

#a new data frame combining the scarHRD score and level of deletion 
Combined_80 = left_join(level_dele, score_80, by= 'sample')

# do correlation analysis between level of deletion versus scarHRD score############
library("ggpubr")

#check normality by ploting or shapiro test
ggplot(data = Combined_80)+
  geom_histogram(mapping = aes(x = level_dele), binwidth = 1)

ggplot(data = Combined_80)+
  geom_histogram(mapping = aes(x = HRD.sum), binwidth = 2)

shapiro.test(Combined$level_dele)
shapiro.test(Combined$HRD.sum)
#From the output, the p-value > 0.05 implying 
#that the distribution of the data are not significantly different from normal distribution. 
#In other words, we can assume the normality.

#Kendall rank correlation test and Spearman rank correlation coefficient
cor.test(Combined_80$level_dele, Combined_80$HRD.sum,  method="kendall")
cor.test(Combined_80$level_dele, Combined_80$HRD.sum,  method="spearman")
#The p-value is less than the significance level alpha = 0.05. 
#conclude that significantly correlated with a correlation coefficient of 0.42 and p-value of 1.675^{-7}

#Visualize data using scatter plots
ggscatter(Combined_80, x = "level_dele", y = "HRD.sum", 
          add = "reg.line", conf.int = TRUE,
          xlab = "level of deletion", ylab = "scarHRD score"
          )+
  stat_cor(method = "kendall", label.x = 11, label.y = 79)

#########################################################################
#correlation analysis between scar HRD score and level of deletion at individual gene level
#determining HRD positive/negative and presence/absence of deletion 
Combined_80_2 = left_join(CN, score_80)
Combined_80_2 = mutate(Combined_80_2, HRD.binary=NA)

for (i in 1:nrow(Combined_80_2)) {
  if(Combined_80_2$HRD.sum[i]>=42){
    Combined_80_2$HRD.binary[i]='positive'
  } else if (Combined_80_2$HRD.sum[i]<42){
    Combined_80_2$HRD.binary[i]='negative'
  }
}

for (i in 1:nrow(Combined_80_2)) {
  if(Combined_80_2$CN_status[i]=='loss of heterozygosity'|Combined_80_2$CN_status[i]=='homozygous deletion'){
    Combined_80_2$CN_status[i]='deletion'
  } else if (Combined_80_2$CN_status[i]=='gain'|Combined_80_2$CN_status[i]=='diploid'){
    Combined_80_2$CN_status[i]='no deletion'
  }
}

#conduct fisher exact test for each gene and output p.values
result_80= matrix(c(1,1,1),ncol = 3,nrow = 16)
colnames(result_80) = c("gene", "p.value_42",'estimates_odds_ratio_42')
k=1
for (i in unique(Combined_80_2$gene)) {
  a=filter(Combined_80_2, gene == i)
  b=xtabs(~CN_status+HRD.binary, data = a)
  b=b[,c(2,1)]
#this code is to exchange the two columns (HRD positive/negative) in contingency table (matrix)
#now  for example estimated odds ratio=10 means 
#that odds of LOH occurring in this gene for scar HRD>=42 subjects is 10 times that for subjects with scar HRD<42
  c=fisher.test(b)
  result_80[k,3] = round(c$estimate, digits = 4)
  result_80[k,2] = round(c$p.value, digits = 4)
  result_80[k,1] = i
  k =k+1
}

#determining statistical significant according to significance level 0.05
result_80=as.data.frame(result_80)
result_80= mutate(result_80, p.value_binary=NA)
for (i in 1:nrow(result_80)) {
  if(result_80$p.value_42[i]>=0.05){
    result_80$p.value_binary[i]='not statistically significant'
  } else if (result_80$p.value_42[i]<0.05){
    result_80$p.value_binary[i]='statistically significant'
  }
}

#fisher exact test for each gene, using 63 as cutoff value for HRD positive/negative
Combined_80_3 = mutate(Combined_80_2, HRD.binary=NA)

for (i in 1:nrow(Combined_80_3)) {
  if(Combined_80_3$HRD.sum[i]>=63){
    Combined_80_3$HRD.binary[i]='positive'
  } else if (Combined_80_3$HRD.sum[i]<63){
    Combined_80_3$HRD.binary[i]='negative'
  }
}

for (i in 1:nrow(Combined_80_3)) {
  if(Combined_80_3$CN_status[i]=='loss of heterozygosity'|Combined_80_3$CN_status[i]=='homozygous deletion'){
    Combined_80_3$CN_status[i]='deletion'
  } else if (Combined_80_3$CN_status[i]=='gain'|Combined_80_3$CN_status[i]=='diploid'){
    Combined_80_3$CN_status[i]='no deletion'
  }
}

#conduct fisher exact test for each gene and output p.values
result_80_2= matrix(c(1,1,1),ncol = 3,nrow = 16)
colnames(result_80_2) = c("gene", "p.value_63",'estimates_odds_ratio_63')
k=1
for (i in unique(Combined_80_3$gene)) {
  a=filter(Combined_80_3, gene == i)
  b=xtabs(~CN_status+HRD.binary, data = a)
  c=fisher.test(b)
  result_80_2[k,3] = round(c$estimate, digits = 4)
  result_80_2[k,2] = round(c$p.value, digits = 4)
  result_80_2[k,1] = i
  k =k+1
}
#error in fisher test using 63 as threshold value because no sample was HRD positive 

write.table(result_80, file = 'fisher_exact_test_gene_level_80samples.csv', sep = ',', row.names = FALSE)

##function test
xtabs(~CN_status+HRD.binary+gene, data = Combined_80_2)
ATM=filter(Combined_80_2, gene=='ATM')
ATM=xtabs(~CN_status+HRD.binary, data = ATM)

