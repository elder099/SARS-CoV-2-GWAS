#install.packages('tidyverse')
library(tidyverse)

#install.packages("devtools") #Has many dependencies
library(devtools)

#install_github("caitiecollins/treeWAS", build_vignettes = TRUE)
library(treeWAS) 



#####
#####GENETIC DATASET
#####

#Load in aligned sequence file
print("Read DNA")
dna<-read.dna(file="~/CDPH/GWAS/GWASv6.0/MAFFT_aligned.fasta",format="fasta") #from "ape" package

print("Create Binary Matrix")
binary_matrix<-DNAbin2genind(dna,polyThres = 0.005)
#View(binary_matrix$tab)

#Create mutation matrix a dataframe
genetic_dataset<-binary_matrix$tab
genetic_dataset<-genetic_dataset[,genetic_dataset[1,]!=1] #Remove all reference snps
#Remove reference sequence entirely
genetic_dataset<-genetic_dataset[rownames(genetic_dataset)!="NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome",]
genetic_dataset<-genetic_dataset[,colMeans(genetic_dataset,na.rm = T)>0.005]  #Remove ultra-low frequency polymorphic snps
#View(genetic_dataset)

#Fix genetic_dataset rownames to be GISAID Acc #'s
rownames(genetic_dataset)<-
  gsub(".*[|](EPI_ISL_\\d*)[|]\\d.*","\\1",rownames(genetic_dataset))




#####
#####PHENOTYPE DATA
#####

#Load in Phenotype data
phen_df<-read.csv("~/CDPH/GWAS/GWASv6.0/Matched_Dem2_Metadata.csv")
#View(phen_df)
print(paste("Total Unique Phenotypes:",length(unique(phen_df$gisaid_accession)) ) )
#just phenotype dataframe with corresponding names
pheno<-as.factor(phen_df$phen1)
names(pheno)<-phen_df$gisaid_accession
#View(as.data.frame(pheno))

###Filtering so all names match -- deprecated
#Remove all rows from genetic dataset not in phenotype
#genetic_dataset<-genetic_dataset[rownames(genetic_dataset) %in% names(pheno),]
#Remove all rows from pheno not in genetic dataset
#pheno<-pheno[names(pheno) %in% rownames(genetic_dataset)]
#length(pheno)



#####
#####TREE BUILDING
#####

###Load in our Auspice tree
tree <- read.tree("~/CDPH/GWAS/GWASv6.0/Newicktree.nwk")
#head(tree$tip.label,30)

#Rename tip labels
tree$tip.label<-gsub("_new[\']*$","",tree$tip.label) #get rid of _new and other issues
tree$tip.label<-gsub("^\'","",tree$tip.label)   #Remove leading apostrophe
tree$tip.label<-gsub(".*[|](EPI_ISL_\\d*)[|]\\d.*","\\1",tree$tip.label) #Trim down to EPI_ISL #'s
#tree$tip.label[1]<-"NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome" #Fix ugly reference name

#Remove unwanted tip labels to match other datasets
rownames(genetic_dataset)[!(rownames(genetic_dataset) %in% tree$tip.label)] 
sum(tree$tip.label %in% rownames(genetic_dataset))
tree_filter<-drop.tip(tree,tree$tip.label[!(tree$tip.label %in% rownames(genetic_dataset))])  #Remove tips not in genetic dataset


###Perform QC so it'll work in TreeWAS
print(paste(
        "TreeWAS Input QC Check:",
        all(tree_filter$tip.label %in% rownames(genetic_dataset)),
        all(rownames(genetic_dataset) %in% tree_filter$tip.label),
        all(tree_filter$tip.label %in% rownames(genetic_dataset)),
        all(rownames(genetic_dataset) %in% names(pheno)),
        all(names(pheno) %in% rownames(genetic_dataset))
      ))

#Yay! It worked! All tips accounted for, no extras





#####
#####RUNNING TREEWAS
#####
#Perform the analysis
out <- treeWAS(snps = genetic_dataset,phen = pheno,tree = tree_filter,seed = 101,n.snps.sim=10*ncol(genetic_dataset),mem.lim=10,correct.prop = T)
#Simulate 10X snps and adjust for unbalanced case/control

print("List All:")
out$treeWAS

#Process TreeWAS output to grab all significant SNPs
term<-out$terminal
sim<-out$simultaneous
subs<-out$subsequent

termsig<-term$sig.snps
if(class(termsig)!="data.frame"){
  termsig<-data.frame()
}
simsig<-sim$sig.snps
if(class(simsig)!="data.frame"){
  simsig<-data.frame()
}
subssig<-subs$sig.snps
if(class(subssig)!="data.frame"){
  subssig<-data.frame()
}
#termsig$Score<-"Terminal"
#simsig$Score<-"Simultaneous"
#subssig$Score<-"Subsequent"

print("Output Each Score to File")
#Output significant SNPs
sigs<-bind_rows(list(Terminal=termsig,
                     Simultaneous=simsig,
                     Subsequent=subssig),.id="Score")
write.csv(sigs,"~/CDPH/GWAS/GWASv6.0/SigSNPs.csv")

#Output all SNPs (significant or not)
pvals<-as.data.frame(term$p.vals)
write.csv(pvals,"~/CDPH/GWAS/GWASv6.0/AllSNPs.csv")

