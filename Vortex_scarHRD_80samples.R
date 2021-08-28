#######2021/7/27 calculate scarHRD from segmentation files using scarHRD package
#install scarHRD package from github
install.packages('devtools')
install.packages('usethis')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("copynumber")

#library
library(usethis)
library(devtools)
install_github('sztup/scarHRD',build_vignettes = TRUE)
library(scarHRD)
library(tidyverse)

SEG = read.csv("test_VOR114_seg.csv")
write.table(SEG, file = 'test_VOR114_seg.txt', sep = '\t')
a = scar_score('test_VOR114_seg.txt', reference = 'grch37', seqz = FALSE, chr.in.names = FALSE)

#######################################for loop
SEG_file = read.csv('ASCAT_seg.CSV')
Ploidy_value = read.csv('Vortex_ploidy_values.CSV')
NEW_SEG = left_join(SEG_file,Ploidy_value)
Sample_list = unique(SEG_file$SampleID)

for (i in Sample_list){
  C = filter(NEW_SEG, SampleID== i)
  write.table(C, file = paste('seg_file/',i,'.txt',sep = ""), sep = '\t')
  d=scar_score(paste('seg_file/',i,'.txt',sep = ""), reference = 'grch37', seqz = FALSE, chr.in.names = FALSE)
  a = rbind(a,d)
}
a = a[-1,]
result = cbind(a,Ploidy_value[,1])
colnames(result)[5] <- "SampleID"

write.table(result, file = 'scarHRD_80samples.csv', sep = ',', row.names = FALSE)

###############
scar=read.csv('scarHRD_80samples.csv')
colnames(scar)[5]=c('sample')
Diagnosis = read.csv('D:/R (UCL project)/Vortex study copy number/second/copy_number_change.csv')
Diagnosis = select(Diagnosis,sample,diagnosis)
Diagnosis=Diagnosis[!duplicated(Diagnosis),]
scar=left_join(Diagnosis,scar)

DIA=c(1)
for (i in unique(scar$diagnosis)) {
  a=filter(scar,diagnosis==i)
  if (nrow(a)>3){
    DIA=append(DIA,i)
  }
}
DIA=DIA[-1]

for (i in 1:nrow(scar)) {
  if(scar$diagnosis[i]%in%DIA){
    scar$diagnosis[i]=scar$diagnosis[i]
  }else{
    scar$diagnosis[i]='other'
  }
}

colnames(scar)[3:5]=c('LOH','TAI','LST')
scar$diagnosis <- factor(scar$diagnosis,levels = c('other',
                                                   "myxoid liposarcoma",
                                                   "leiomyosarcoma",
                                                   "malignant peripheral nerve sheath tumour",
                                                   "myxofibrosarcoma",
                                                   "undifferentiated pleomorphic sarcoma"
                                                   ))
ggplot(data = scar)+
  geom_point(mapping = aes(x = diagnosis, y = HRD.sum))+
  theme(axis.text.x = element_text(angle = 15,vjust =1, hjust=0.75, size = 12))

DP=ggplot(data = scar)+
  geom_boxplot(mapping = aes(x = diagnosis, y = HRD.sum))+
  theme(axis.text.y = element_text(angle = 0,vjust =0.5, hjust=1, size = 12))+
  coord_flip()+
  theme(axis.text.y = element_text(angle = 0,vjust =0.5, hjust=1, size = 12))

ggplot(data = scar,aes(x = diagnosis, y = HRD.sum))+
  geom_dotplot(binaxis = 'y',stackdir = 'center',binwidth = 1,dotsize = 0.8)+
  theme(axis.text.y = element_text(angle = 0,vjust =0.5, hjust=1, size = 12))+
  coord_flip()+
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),geom="crossbar", width=0.3)

ggplot(data = scar,aes(x = diagnosis, y = HRD.sum,fill=diagnosis))+
  geom_boxplot(width=0.35,fill='white')+
  geom_dotplot(binaxis = 'y',stackdir = 'center',binwidth = 1,dotsize = 0.85)+
  theme(axis.text.y = element_text(angle = 0,vjust =0.5, hjust=1, size = 12))+
  coord_flip()+
  theme(legend.position="none")

ggplot(data = scar,aes(x = diagnosis, y = HRD.sum,fill=diagnosis))+
  geom_boxplot(width=0.35,fill='white')+
  geom_dotplot(binaxis = 'y',stackdir = 'center',binwidth = 1,dotsize = 0.85)+
  theme(axis.text.y = element_text(angle = 0,vjust =0.5, hjust=1, size = 12))+
  theme(legend.position="none")

##########compare scarHRD scores between subtypes 
P_v= pairwise.wilcox.test(scar$HRD.sum, scar$diagnosis,
                     p.adjust.method = "BH")
P_v=P_v$p.value

################add P value to plot
P_v <- tibble::tribble(
  ~group1,     ~group2,   ~p.adj, ~p.adj.signif,
  "other",     "myxofibrosarcoma", 0.01383689, '*',
  "other",     "undifferentiated pleomorphic sarcoma", 0.01297034,'*',
  "myxoid liposarcoma",     "leiomyosarcoma", 0.016193307,'*',
  "myxoid liposarcoma",     "myxofibrosarcoma", 0.000118018,'***',
  "myxoid liposarcoma",     "undifferentiated pleomorphic sarcoma", 0.000118018,'***'
  )

library(tidyverse)
library(ggpubr)
library(rstatix)

P_v <- P_v %>%
  mutate(y.position = c(60,64,48,52,56))
DP+stat_pvalue_manual(P_v, label = "p.adj.signif", tip.length = 0.01)

###############################
#table: samples size for each STS subtype in this study
Diagnosis_table = data.frame()
Diagnosis_table[1:17,1]=unique(Diagnosis$diagnosis)
Diagnosis_table[1:17,2]=NA
colnames(Diagnosis_table)=c('diagnosis','n')
for (i in 1:nrow(Diagnosis_table)) {
  Diagnosis_table$n[i]=nrow(filter(Diagnosis,diagnosis==Diagnosis_table$diagnosis[i]))
}
write.csv(Diagnosis_table, file = 'sample size of STS subtypes.csv')
