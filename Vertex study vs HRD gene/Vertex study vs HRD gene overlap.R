## determine target HR related gene
## Latest update: 25/05/2021
## Version 1.0.0
## Author: Dian Lyu

Vortex_study = read.table("Vortex_study_V1_Exons_10bp_Targets.txt")
HRD_gene = read.csv("HRD_gene_list.csv")
HRD_gene[1,] %in% Vortex_study[,]

#select overlap genes between the two gene lists 
HRD_gene [,] %in% Vortex_study[,]
Overlap = HRD_gene[HRD_gene [,] %in% Vortex_study[,],]
Overlap= data.frame(Overlap)
(Overlap)

#select 
not_overlap = HRD_gene[HRD_gene [,] %in% Vortex_study[,]==F,]
not_overlap = data.frame(not_overlap)
