# VOR114 = read.csv("VOR114.txt")
# VOR114$sample = "VOR114"
# VOR125 = read.csv("VOR125.txt")
# VOR125$sample = "VOR125"
# VOR135 = read.csv("VOR135.txt")
# VOR135$sample = "VOR135"
# VOR231 = read.csv("VOR231.txt")
# VOR231$sample = "VOR231"
# VOR253 = read.csv("VOR253.txt")
# VOR253$sample = "VOR253"
# VOR267 = read.csv("VOR267.txt")
# VOR267$sample = "VOR267"
# VOR67 = read.csv("VOR67.txt")
# VOR67$sample = "VOR67"
# VOR73 = read.csv("VOR73.txt")
# VOR73$sample = "VOR73"
# VOR8 = read.csv("VOR8.txt")
# VOR8$sample = "VOR8"
#total = rbind(VOR114,VOR125,VOR135,VOR231,VOR253,VOR267,VOR67,VOR73,VOR8)

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



