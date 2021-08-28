library(tidyverse)
##################load pathogenic mutations
M = read.csv('D:/R (UCL project)/Vortex study MAF/Vortex study mutation calls MuTect2/variants_HR_related.csv')
M = select(M,Hugo_Symbol,Variant_Classification,Tumor_Sample_Barcode,IMPACT,SIFT,PolyPhen)
M = rename(M,sample=Tumor_Sample_Barcode,gene=Hugo_Symbol)
M = mutate(M,pathogenicity=NA)
for (i in 1:nrow(M)) {
  if(grepl("tolerated",M$SIFT[i])&grepl("benign",M$PolyPhen[i])){
    M$pathogenicity[i]='non pathogenic'
  }else{
    M$pathogenicity[i]='pathogenic'
  }
}

#load copy number status
CN = read.csv('copy_number_change.csv')
CN = select(CN, sample,gene,CN_status)
M_CN=left_join(M,CN)
#copy number status absent for four samples
M_CN= filter(M_CN, !is.na(CN_status))
M_CN=mutate(M_CN,Mutation_CN=NA)

###########################################################################
#categorize samples as bi-allelic mutated if they had homozygous deletion or
#LOH in combination with pathogenic mutations in at least one of 16 HR pathway genes
#as mono-allelic if they has no pathogenic mutation but LOH or pathogenic mutation without LOH
#as WT if they had no pathogenic mutation or LOH or homozygous deletion in any 16 HRD genes
for (i in 1:nrow(M_CN)) {
  if(M_CN$pathogenicity[i]=='pathogenic'& M_CN$CN_status[i]=='loss of heterozygosity'){
    M_CN$Mutation_CN[i]='bi-allelic'
  } else if(
    M_CN$pathogenicity[i]=='pathogenic'&!M_CN$CN_status[i]=='loss of heterozygosity'
  ){
    M_CN$Mutation_CN[i]='mono-allelic'
  }else{
    M_CN$Mutation_CN[i]='WT'
  }
}
Bi=filter(M_CN,Mutation_CN=='bi-allelic')
unique(Bi$sample)
HomoDele=filter(CN,CN_status=='homozygous deletion')
Bi=c(unique(Bi$sample),HomoDele$sample)#a list of samples that have at least one bi-allelic mutated HRD gene


Deletion = count(
  group_by(CN, sample), 
  CN_status == 'loss of heterozygosity'|CN_status =="homozygous deletion")
colnames(Deletion) <- c("sample", "deletion","level_dele")
Y_deletion = filter(Deletion, deletion == TRUE)
N_deletion = filter(Deletion, deletion == FALSE)
Non_LOH = N_deletion[! N_deletion$sample %in% Y_deletion$sample,]
Non_LOH=c(unique(Non_LOH$sample))
WT = Non_LOH[!Non_LOH%in%M_CN$sample]#a list of samples 
#that have neither LOH/homozygous deletion nor pathogenic mutations in any 16 HRD genes

Total_sample=unique(CN$sample)
Mono=Total_sample[!Total_sample%in%Bi&!Total_sample%in%WT]#a list of samples 
#that have either LOH or pathogenic mutations in at least one of 16 HRD gene

################################################### Mann whitney / Fisher's exact test
#load scar HRD score for 80 samples
scar=read.csv('scarHRD_80samples.csv')
scar=select(scar, HRD.sum,SampleID)
colnames(scar)=c('HRD_score','sample')
scar=mutate(scar,CN_M_status=NA)
##divide samples into two groups: bi-allelic mutated and non-biallelic (mono-allelic and WT)
for (i in 1:nrow(scar)) {
  if(scar$sample[i]%in%Bi){
    scar$CN_M_status[i]='bi-allelic'
  }else if(
    scar$sample[i]%in%Mono|scar$sample[i]%in%WT
  ){
    scar$CN_M_status[i]='no_bi-allelic'
  }
}

#Mann-Whitney test comparing two groups
library("ggpubr")
wilcox.test(HRD_score~CN_M_status,data=scar)

#Contingency table and Fisher's exact test
scar_binary=mutate(scar,HRD_binary_42=NA)
for (i in 1:nrow(scar_binary)) {
  if(scar_binary$HRD_score[i]>=42){
    scar_binary$HRD_binary_42[i]='positive'
  } else {
    scar_binary$HRD_binary_42[i]='negative'
  }
}
scar_binary=scar_binary[,2:4]

Con_table=xtabs(~HRD_binary_42+CN_M_status,scar_binary)
Con_table=Con_table[c(2,1),]
fisher.test(Con_table)



