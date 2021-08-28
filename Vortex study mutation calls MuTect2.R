library(maftools)
library(tidyverse)
#####################################
#load 80 MAF files for all samples and combined them into one dataframe
list = list.files("Dian_MAFs", pattern=".maf", full.names=T)

a = read.maf(maf = 'Dian_MAFs/VOR11.mutect2.vep.maf')
a = a@data
b = a[1,]
for (i in list){
  C = read.maf(i)
  C = C@data
  b = rbind(b,C,fill=TRUE)
}
combined_files = b[-1,]

#filter rule for only 16 HR pathway genes
HR_genes = c('ATM','ATRX','BLM','BRCA1','BRCA2','BRIP1','CHEK2','MRE11A','NBN','PALB2','POLD1','RAD50','RAD51',
             'RAD51B','RAD51C','RECQL4')

#filter DKFZBias, set it to 'blank'
DKFZ = c('strand','damage','damage,strand')

#filter rule for only mutations in exons that effect protein 
LOSS_OF_FUNCTION = c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins',
                     'Missense_Mutation','Nonsense_Mutation','Splice_Site','Silent')

#apply filter rules to the dataframe, one result of variants in 16 HR genes, another one of variants in non HR genes
# GET high confidence mutations and somatic mutations 
PASS = filter(combined_files, FILTER == 'PASS')
PASS_BLANK = filter(PASS,!(DKFZBias%in%DKFZ))
CODING_region = filter(PASS_BLANK,Variant_Classification%in%LOSS_OF_FUNCTION)
HR_related = filter(CODING_region, Hugo_Symbol%in%HR_genes)
non_HR_related = filter(CODING_region, !(Hugo_Symbol%in%HR_genes))


HR_related_pathogenic= mutate(HR_related,pathogenicity=NA)
for (i in 1:nrow(HR_related_pathogenic)) {
  if(grepl("tolerated",HR_related_pathogenic$SIFT[i])&grepl("benign",HR_related_pathogenic$PolyPhen[i])){
    HR_related_pathogenic$pathogenicity[i]='non pathogenic'
  }else{
    HR_related_pathogenic$pathogenicity[i]='pathogenic'
  }
}
HR_related_pathogenic=filter(HR_related_pathogenic,pathogenicity=='pathogenic')

#output dataframe
write.table(HR_related, file = 'variants_HR_related.csv', sep = ',', row.names = FALSE)

#Visualization
Result_MAF = read.maf(HR_related)
Result_MAF = read.maf(HR_related_pathogenic)

getSampleSummary(Result_MAF)
getGeneSummary(Result_MAF)
getClinicalData(Result_MAF)
getFields(Result_MAF)

NON_HR_MAF = read.maf(non_HR_related)

plotmafSummary(Result_MAF,titvRaw = FALSE)
oncoplot(maf = Result_MAF,
         legendFontSize= 1.7,
         showTumorSampleBarcodes = TRUE,barcodeSrt=30,SampleNamefontSize=1.2)

plotmafSummary(NON_HR_MAF,titvRaw = FALSE)
oncoplot(maf = NON_HR_MAF,
         legendFontSize= 1.7,
         showTumorSampleBarcodes = TRUE,barcodeSrt=75,SampleNamefontSize=1.1)



#####################################################
#functions that filter MAF files
#SUBSET = subsetMaf(combined_files, genes = HR_genes)%>%
#subsetMaf(query = "FILTER == 'PASS'")

gene = unique(combined_files$Hugo_Symbol)
write.table(gene, file = 'gene_MAF.csv', sep = ',', row.names = FALSE)

#####################################################
CN = read.csv('D:/R (UCL project)/Vortex study copy number/second/copy_number_change.csv')
CN = select(CN,sample,gene,CN_status)
CN_2 = filter(CN,CN_status=='loss of heterozygosity'|CN_status=='homozygous deletion')
CN_a=CN_2
CN_a[,1]=CN_2[,2]
CN_a[,2]=CN_2[,1]
colnames(CN_a)[1:2]=c('gene','sample')

#filter out somatic variants in samples that have no scarHRD scores and copy number data
#下一步，从MAF data中删除VOR26，VOR50，VOR95,无scarHRD score 和copy number data
f=unique(HR_related_pathogenic$Tumor_Sample_Barcode)
a=unique(CN$sample)
f[!f%in%a]#注意这里没有冒号和后面的 “a[!a%in%b]” 不一样
HR_related_pathogenic=filter(HR_related_pathogenic,!Tumor_Sample_Barcode%in%c('VOR26','VOR50','VOR95'))

#identify 6 samples that had no deletion (LOH or homozygous deletion) and manually add these samples into MAF file
#往MAF data中加入六个没有copy number deletion的样本，以便在图中将他们显示出来
a=unique(CN$sample)
b=unique(CN_2$sample)
a[!a%in%b]

d=HR_related_pathogenic[1:6,]
d[1:6,]=NA
d$Tumor_Sample_Barcode[1:6]=c(a[!a%in%b])
HR_related_pathogenic=rbind(d,HR_related_pathogenic)



#load scar HRD scores
scar = read.csv('D:/R (UCL project)/Vortex study scarHRD_80samples/scarHRD_80samples.csv')
scar = select(scar,HRD.sum,SampleID)
colnames(scar)=c('scarHRD','Tumor_Sample_Barcode')

#samples order according to scarHRD scores from the highest to the lowest
#最后sample顺序按scarHRD score由高到低排列
scar=scar[order(-scar$scarHRD),]

#combine maf, copy number data and scarHRD together
HR_related_plus_cn_scar = read.maf(maf = HR_related_pathogenic,
                              cnTable = CN_a,
                              clinicalData = scar)

col = c('Splice_Site'="#FF7F00",
        'Frame_Shift_Del'="#1F78B4" ,
        'Missense_Mutation'="#33A02C",
        'Nonsense_Mutation'="#FF7F00",
        'In_Frame_Del'="#FFFF99",
        'loss of heterozygosity'="#A6CEE3",
        'homozygous deletion' = "#6A3D9A",
        'Multi_Hit'="#FB9A99"
        )

oncoplot(maf = HR_related_plus_cn_scar,topBarData = 'scarHRD',
         showTumorSampleBarcodes = TRUE,barcodeSrt=75,
         SampleNamefontSize=0.9,
         sampleOrder=scar$Tumor_Sample_Barcode,
         removeNonMutated = FALSE,
         colors = col)

lolli=read.maf(HR_related_pathogenic[7:23,])

#lollipop plot of ATRX pathogenic somatic variants
lollipopPlot(maf = lolli, gene = 'ATRX',AACol='HGVSp_Short', refSeqID ='NM_000489')


############################################################# no HR related gene
n_HR_p= mutate(non_HR_related,pathogenicity=NA)
for (i in 1:nrow(n_HR_p)) {
  if(grepl("tolerated",n_HR_p$SIFT[i])&grepl("benign",n_HR_p$PolyPhen[i])){
    n_HR_p$pathogenicity[i]='non pathogenic'
  }else{
    n_HR_p$pathogenicity[i]='pathogenic'
  }
}
n_HR_p=filter(n_HR_p,pathogenicity=='pathogenic')

#output dataframe
write.table(n_HR_p, file = 'variants_non_HR_related_pathogenic.csv', sep = ',', row.names = FALSE)

#Mann whitney test
TP53=filter(n_HR_p,Hugo_Symbol=='TP53')
TP53=select(TP53,Tumor_Sample_Barcode)
s=mutate(scar,TP53=NA)
for (i in 1:nrow(s)) {
  if(s$Tumor_Sample_Barcode[i]%in%TP53$Tumor_Sample_Barcode){
    s$TP53[i]='positive'
  }else{
    s$TP53[i]='negative'
  }
}
library("ggpubr")
wilcox.test(scarHRD~TP53,data=s)

RB1=filter(n_HR_p,Hugo_Symbol=='RB1')
RB1=select(RB1,Tumor_Sample_Barcode)
s=mutate(scar,RB1=NA)
for (i in 1:nrow(s)) {
  if(s$Tumor_Sample_Barcode[i]%in%RB1$Tumor_Sample_Barcode){
    s$RB1[i]='positive'
  }else{
    s$RB1[i]='negative'
  }
}
library("ggpubr")
wilcox.test(scarHRD~RB1,data=s)


