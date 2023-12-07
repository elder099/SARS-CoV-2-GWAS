library(tidyverse)


###Load in the datasets
#Fill-in-the-Blanks IGED
IGED <- read.csv("~/CDPH/GWAS/GWASv6.0/IGED_Fill_in_R.csv")
#View(IGED)

#CR-merged IGED
CR_match<- read.csv("~/CDPH/GWAS/GWASv6.0/IGED_CR_Merge_Mask.csv")
#View(CR_match)


###Merge "filled-in" IGED and CR-merged IGED
joined<- inner_join(IGED,CR_match,by=c("incidentid")) #Will create some dupes
#View(joined)


###Data Cleaning
#remove unnecessary columns
joined_trim <- joined[ c(3:5,8,12,24:dim(joined)[2]) ]
colnames(joined_trim)<-c("gisaid_accession","greek","lineage","incidentid","sample_id",
                         "phen1","phen2","phen3","dem1","dem2","dem3")

#Filter to only include those with real gisaid_accessions
IGED_GIS<-joined_trim[!is.na(joined_trim$gisaid_accession) & grepl("EPI_ISL_",joined_trim$gisaid_accession),]
#Remove ".2's" from gisaid #'s
IGED_GIS$gisaid_accession<-gsub("^(EPI_ISL_\\d+)[.]2","\\1",IGED_GIS$gisaid_accession)
#Remove NA's in dem2
IGED_GIS<-IGED_GIS[!is.na(IGED_GIS$dem2),]
#dim(IGED_GIS)

###Output data to file
write.csv(IGED_GIS,"~/CDPH/GWAS/GWASv6.0/IGED_Full_Phenotype.csv",row.names = FALSE)                    
