#####
#####MATCHING AND SUBSETTING ON DEMOGRAPHICS
#####
#install.packages("MatchIt")
library("MatchIt")
library(tidyverse)

#####Load in the data
#IGED<-read.csv("~/CDPH/GWAS/IGED_w_GISandHosp.csv")
IGED<-read.csv("~/CDPH/GWAS/GWASv6.0/IGED_Full_Phenotype.csv") #Rename this for standardization
#View(IGED)



###Perform the matching
set.seed(101) #Set seed for reproducibility
matched<-matchit(phen1 ~ dem2,data=IGED,method="nearest",ratio = 3)
#summary(matched)

#Turn data into data.frame
match_data<-match.data(matched)
#Output dataset for future use
write.csv(match_data,"~/CDPH/GWAS/GWASv6.0/Matched_Dem2_Metadata.csv",row.names = FALSE)






#####
#####GISAIDR DOWNLOADS
#####

library(devtools)
#devtools::install_github("Wytamma/GISAIDR")
library(GISAIDR)

credentials<-login(username="username",password="password",database="EpiCoV")


#####
#####Download matched data from GISAID
#####

###Data Cleaning
matched_gis<-match_data$gisaid_accession #Only take this column


#Deduplicate so we don't download dupes
length(unique(matched_gis))
matched_gis<-matched_gis[!duplicated(matched_gis)]
length(matched_gis)


#####
#####For loop to automate download of sequences
#####No limit for downloads

#empty dataframe to fill
df<-data.frame()
chunk<-3000             #Size of chunk to download each loop
acclen<-length(matched_gis) #Number of total accessions to download

for(i in 1:ceiling(acclen/chunk)){
  if(chunk*i<acclen){
    acclist<-matched_gis[(1+(i-1)*chunk):(1+i*chunk)] #before end of column
  } else{
    acclist<-tail(matched_gis,acclen%%chunk) #end of column
  }
  #Download the subsets
  df1<-download(
    credentials=credentials,
    list_of_accession_ids = acclist
  )
  #Concatenate the subsets
  df<-rbind(df,df1)
}

#View(df)


#The number of unique accession numbers matches
#length(unique(matched_gis))
#length(unique(df$accession_id))
#matched_gis[!(matched_gis %in% df$accession_id)] 
#dim(df)

#Remove duplicates
df_dedup<-df[!duplicated(df$accession_id),]
#dim(df_dedup)

export_fasta(df_dedup,
             out_file_name = '~/CDPH/GWAS/GWASv6.0/GISAID_Seqs.fasta',
             delimiter = '|')


