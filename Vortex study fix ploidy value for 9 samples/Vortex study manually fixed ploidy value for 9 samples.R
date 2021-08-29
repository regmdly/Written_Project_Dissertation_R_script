## Ploidy value calculation
## Latest update: 16/07/2021
## Version 1.0.0
## Author: Dian Lyu

#load ploidy files
D=data.frame(IOD=0,GALLERY=0,sample=0)
for (i in list.files(path = "VOR_9_samples")) {
  C=read.csv(i)
  C$sample=gsub('[.txt]', '', i)
  D = rbind(D,C)
}
total = D[-1,]

#calculate mean of IOD normal
library(tidyverse)
IOD_normal = filter(total,GALLERY == 1 | GALLERY == 2)
IOD_normal = group_by(IOD_normal,sample)%>%
  summarise(IOD_mean=mean(IOD))

#create density plots for IOD tumor
IOD_tumor = filter(total, GALLERY == 0 )
for (i in unique(IOD_tumor$sample)) {
  a=filter(IOD_tumor, sample == i)
  plot(density(a$IOD))
}

#calculate IOD tumor with maximal frequency
result= c(1)
for (i in unique(IOD_tumor$sample)) {
  a=filter(IOD_tumor, sample == i)
  b=density(a$IOD)
  result[i] =b$x[which.max(b$y)]
}
result = result[-1]

#integrate 9 IOD tumour max value into one dataframe
IOD_max_tumor= data.frame(result)
IOD_max_tumor$sample <- rownames(IOD_max_tumor)
names(IOD_max_tumor) <- c("IOD_max_tumor", "sample")

#combine IOD normal mean and IOD tumor max
IOD = left_join(IOD_normal,IOD_max_tumor)

#calculate ploidy value for each sample
IOD= mutate(IOD, ploidy_value = (IOD_max_tumor/IOD_mean)*2)



