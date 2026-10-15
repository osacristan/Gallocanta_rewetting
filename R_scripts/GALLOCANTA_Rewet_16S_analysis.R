#########################################################
### A) Installing and loading required packages
#########################################################

getRversion()
#[1] ‘4.3.2’

## cleaning buffer
rm(list=ls(all=TRUE))

install.packages("here")
library(here)

setwd(here("data/"))

library(knitr)

library(fitdistrplus) # fit distributions

#install.packages("Rcpp", dependencies = TRUE)
#install.packages("fastmap")
#library(devtools)
#devtools::install_github('ggloor/CoDaSeq/CoDaSeq')
library(CoDaSeq)

library(zCompositions)
library(igraph)
if (!require("car")) {
  install.packages("car", dependencies = TRUE)
  library(car)
}
library(grDevices)

#if (!require("propr")) {
 # install.packages("propr", dependencies = TRUE)
  #library(propr)
#}

#devtools::install_github("tpq/propr")
library(propr)

library(vegan)

# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")

library(qiime2R)

library(ggplot2)
if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}
library(ggsignif)
library(RColorBrewer)
library(dendextend)
library(ggdendro)

#BiocManager::install("microbiome")
library(microbiome)
library(phyloseq)

#########################################################
### B) Reading in data
########################################################

## REWETTING ##
# Data: import qiime2 qza
qzaobj <- read_qza('gallocanta_rewetting_asvtable.qza')
names(qzaobj)
asvtable <- qzaobj$data

dim(asvtable) # [1] 6499   47
#filt <- asvtable[rowSums(asvtable)>0,]
#dim(filt) # [1] 6499   47 # No ASVs with 0 across all samples

taxonomy <- read_qza('gallocanta_rewetting_tax_silva138.qza')
head(taxonomy$data)
taxonomy<-parse_taxonomy(taxonomy$data)
head(taxonomy)
class(taxonomy)

metadata <- read.table("Gallocanta_metadataRewetting_R.txt", header=T, row.names=1, sep="\t")
View(metadata)
dim(metadata)
row.names(metadata)
str(metadata)


# Add a column with Time and Plot
metadata$TimePlot <- paste(metadata$Time,"_", metadata$Plot)
metadata$TreatTime <- paste(metadata$Treatment,"_", metadata$Time)
metadata$TypeTreat <- paste(metadata$Type,"_", metadata$Treatment)


tree <- read_qza('gallocanta_rewetting_align2tree/rooted_tree.qza')

metadata.f.av <- read.table("Gallocanta_metadataRewetting_R_filtAv.txt", header=T, row.names=1, sep="\t")
metadata.f.av$Sample <- rownames(metadata.f.av)
metadata.f.av$TreatTimeType <- paste(metadata.f.av$Treatment,metadata.f.av$Time,metadata.f.av$Type) 
metadata.f.av$TreatType <- paste(metadata.f.av$Treatment,metadata.f.av$Type) 
metadata.f.av$TreatTime <- paste(metadata.f.av$Treatment,metadata.f.av$Time)

metadata.DNA.f.av <- metadata.f.av[metadata.f.av$Type=="DNA",]
metadata.RNA.f.av <- metadata.f.av[metadata.f.av$Type=="cDNA",]

# order levels of metadata variables
# transform to a ordered factor vector. With gsub("\\s*", "") you remove extra spaces from the column avoiding the appearance of NAs due to a mismatch btw the levels of the vector and the ones you input.
metadata.DNA.f.av$Treatment<- factor(metadata.DNA.f.av$Treatment, levels = c("Dry","Rewetted"))
metadata.DNA.f.av$Time <- factor(metadata.DNA.f.av$Time, levels = c("T1","T3","T7","T11"))
metadata.DNA.f.av$Sample <- factor(metadata.DNA.f.av$Sample, levels = c("Dry T1","Dry T3","Dry T7","Dry T11","Rewetted T1","Rewetted T3","Rewetted T7","Rewetted T11"))
rownames(metadata.DNA.f.av) <- factor(rownames(metadata.DNA.f.av), levels = c("Dry T1","Dry T3","Dry T7","Dry T11","Rewetted T1","Rewetted T3","Rewetted T7","Rewetted T11"))


metadata.RNA.f.av$Treatment<- factor(metadata.RNA.f.av$Treatment, levels = c("Dry","Rewetted"))
metadata.RNA.f.av$Time <- factor(metadata.RNA.f.av$Time, levels = c("T1","T3","T7","T11"))
metadata.RNA.f.av$Sample <- factor(metadata.RNA.f.av$Sample, levels = c("cDry T1","cDry T3","cDry T7","cDry T11","cRewetted T1","cRewetted T3","cRewetted T7"))

# change row names metadata by DNA Sample data
rownames(metadata.RNA.f.av) <- metadata.DNA.f.av$Sample[1:7]


## DATA FILTRATION ##

# asvtable DNA
asvtableDNA <- asvtable[,19:47]
dim(asvtableDNA) #[1] 6499   29
asvtableDNA <- asvtableDNA[rowSums(asvtableDNA)>0,] # filter ASV with 0 across all samples
dim(asvtableDNA) #[1] 5861   29
colnames(asvtableDNA)
head(asvtableDNA)
asvDNA <- rownames(asvtableDNA)

# asvtable cDNA
asvtableRNA <- asvtable[,1:18]
dim(asvtableRNA) #[1] 6499   18
asvtableRNA <- asvtableRNA[rowSums(asvtableRNA)>0,] # filter ASV with 0 across all samples
dim(asvtableRNA) #[1] 2512   18
asvRNA <- rownames(asvtableRNA)

# taxonomy DNA
taxonDNA <- taxonomy[asvDNA,]
dim(taxonDNA) # [1] 5861    7

# taxonomy RNA
taxonRNA <- taxonomy[asvRNA,]
dim(taxonRNA) # [1] 2512    7


# Metadata DNA
metadata.DNA <- metadata[metadata$Type == "DNA",]

# Metadata RNA
metadata.RNA <- metadata[metadata$Type == "cDNA",]



# BLOCK1
## This is the entire dataset point
# filter the dataset:min reads per sample = 7000; min proportional abundance of a read in any sample=0.1%; min sample fraction of non-0 reads for each feature in all samples=5%)
f <- codaSeq.filter(asvtable, min.reads = 7000, min.prop=0.001, min.occurrence=0.05, samples.by.row=FALSE)
dim(f) # [1] 986  41 features x samples
class(f)
head(f)
colnames(f)

# filter metadata
metadata.f <- metadata[colnames(f),]
 

# check the order of samples
order_idx <- match(rownames(metadata.f), colnames(f))

# Reorder filtered data frame to have the sample names in the same order as the metadata data frame
f  <- f[ , order_idx]


# VERIFY SAMPLE ORDER!!!
all(rownames(metadata.f) == colnames(f)) #[1] TRUE

# NOW Change sample names for more informative labels
#colnames(f) <- metadata$SampleLocality

#write.table(f, "Gallocanta_Rewetting_ASV_0.1%.txt", sep="\t")

## DNA DATA
# filter the dataset
f.DNA <- f[,13:41]
dim(f.DNA) #[1] 986  29
f.DNA <- f.DNA[rowSums(f.DNA)>0,] # filter out ASVs with 0 counts across samples
dim(f.DNA) #[1] 982  29 features x samples

metadata.DNA.f <- metadata.DNA[colnames(f.DNA),]

## RNA DATA
# filter the dataset
f.RNA <- f[,1:12]
dim(f.RNA) #[1] 986  12
f.RNA <- f.RNA[rowSums(f.RNA)>0,] # filter out ASVs with 0 counts across samples
dim(f.RNA) #[1] 804  12 features x samples

metadata.RNA.f <- metadata.RNA[colnames(f.RNA),]

# ORDER levels of metadata variables
# transform to a ordered factor vector. With gsub("\\s*", "") you remove extra spaces from the column avoiding the appearance of NAs due to a mismatch btw the levels of the vector and the ones you input.
metadata.f$Time <- factor(metadata.f$Time, levels = c("T1","T3","T7","T11"))
metadata.f$Treatment <- factor(metadata.f$Treatment, levels = c("Dry","Rewetted"))
metadata.f$TreatTime <- factor(metadata.f$TreatTime,
                               levels = c("Dry _ T1","Dry _ T3","Dry _ T7","Dry _ T11",
                                          "Rewetted _ T3","Rewetted _ T7","Rewetted _ T11"))
metadata.f$Type <- factor(metadata.f$Type, levels = c("DNA","cDNA"))
metadata.f$TimePlot <- factor(metadata.f$TimePlot,
                                  levels = c("T1 _ D3","T1 _ D4","T1 _ D5","T3 _ D1","T7 _ D1","T7 _ D2",
                                             "T7 _ D3","T7 _ D4","T7 _ D5","T11 _ D1","T11 _ D3","T11 _ D4",
                                             "T11 _ D5",
                                             "T1 _ R6","T1 _ R7","T1 _ R8","T1 _ R9","T3 _ R6","T3 _ R7",
                                             "T3 _ R8","T3 _ R9","T7 _ R6","T7 _ R7","T7 _ R8","T7 _ R9",
                                             "T7 _ R10","T11 _ R6","T11 _ R7","T11 _ R8","T11 _ R9",
                                             "T11 _ R10"))


metadata.f$Sample <- rownames(metadata.f)
levels(factor(metadata.f$Sample))
metadata.f$Sample <- factor(metadata.f$Sample, levels = c("D3T1","D4T1","D5T1","D1T3","D1T7","D2T7","D3T7","D4T7","D5T7","D1TF","D3TF","D4TF",
                                                          "cD3T1","cD1T3","cD1T7","cD5T7","cD5TF",
                                                          "R6T1","R7T1","R9T1","R6T3","R7T3","R8T2","R9T3","R6T7","R7T7","R8T7","R9T7","R10T7","R6TF","R7TF","R8TF","R9TF","R10TF",
                                                          "cR7T1","cR8T1","cR6T3","cR7T2","cR8T3","cR9T3","cR9T7"))

metadata.f <- metadata.f[order(metadata.f$Sample),]
levels(metadata.f$Sample)

# order levels of metadata variables
# transform to a ordered factor vector. With gsub("\\s*", "") you remove extra spaces from the column avoiding the appearance of NAs due to a mismatch btw the levels of the vector and the ones you input.
metadata.DNA.f$Time <- factor(metadata.DNA.f$Time, levels = c("T1","T3","T7","T11"))
metadata.DNA.f$Treatment <- factor(metadata.DNA.f$Treatment, levels = c("Dry","Rewetted"))
metadata.DNA.f$TreatTime <- factor(metadata.DNA.f$TreatTime,
                                   levels = c("Dry _ T1","Dry _ T3","Dry _ T7","Dry _ T11",
                                              "Rewetted _ T3","Rewetted _ T7","Rewetted _ T11"))

metadata.DNA.f$TimePlot <- factor(metadata.DNA.f$TimePlot,
                                  levels = c("T1 _ D3","T1 _ D4","T1 _ D5","T3 _ D1","T7 _ D1","T7 _ D2",
                                             "T7 _ D3","T7 _ D4","T7 _ D5","T11 _ D1","T11 _ D3","T11 _ D4",
                                             "T1 _ R6","T1 _ R7","T1 _ R9","T3 _ R6","T3 _ R7","T3 _ R8",
                                             "T3 _ R9","T7 _ R6","T7 _ R7","T7 _ R8","T7 _ R9",
                                             "T7 _ R10","T11 _ R6","T11 _ R7","T11 _ R8","T11 _ R9","T11 _ R10"))


metadata.DNA.f$Sample <- rownames(metadata.DNA.f)
levels(factor(metadata.DNA.f$Sample))
metadata.DNA.f$Sample <- factor(metadata.DNA.f$Sample, levels = c("D3T1","D4T1","D5T1","D1T3","D1T7","D2T7","D3T7","D4T7","D5T7","D1TF","D3TF","D4TF",
                                                                 "R6T1","R7T1","R9T1","R6T3","R7T3","R8T2","R9T3","R6T7","R7T7","R8T7","R9T7","R10T7","R6TF","R7TF","R8TF","R9TF","R10TF"))

metadata.DNA.f <- metadata.DNA.f[order(metadata.DNA.f$Sample),]
levels(metadata.DNA.f$Sample)

# order levels of metadata variables
# transform to a ordered factor vector. With gsub("\\s*", "") you remove extra spaces from the column avoiding the appearance of NAs due to a mismatch btw the levels of the vector and the ones you input.
metadata.RNA.f$Time <- factor(metadata.RNA.f$Time, levels = c("T1","T3","T7","T11"))
metadata.RNA.f$Treatment <- factor(metadata.RNA.f$Treatment, levels = c("Dry","Rewetted"))
metadata.RNA.f$TreatTime <- factor(metadata.RNA.f$TreatTime,
                                   levels = c("Dry _ T1","Dry _ T3","Dry _ T7","Dry _ T11",
                                              "Rewetted _ T3","Rewetted _ T7"))
metadata.RNA.f$TimePlot <- factor(metadata.RNA.f$TimePlot,
                                  levels = c("T1 _ D3","T3 _ D1","T7 _ D1","T7 _ D5","T11 _ D5",
                                             "T1 _ R7","T1 _ R8","T3 _ R6","T3 _ R7","T3 _ R8","T3 _ R9","T7 _ R9"))


metadata.RNA.f$Sample <- rownames(metadata.RNA.f)
levels(factor(metadata.RNA.f$Sample))
metadata.RNA.f$Sample <- factor(metadata.RNA.f$Sample, levels = c("cD3T1","cD1T3","cD1T7","cD5T7","cD5TF",
                                                                  "cR7T1","cR8T1","cR6T3","cR7T2","cR8T3","cR9T3","cR9T7"))

metadata.RNA.f <- metadata.RNA.f[order(metadata.RNA.f$Sample),]
levels(metadata.RNA.f$Sample)


# Environmental data
envdata <- read.table("Environmental_data_rewetting.txt", header=T, sep="\t")

str(envdata)
rownames(envdata) <- gsub(" ","",paste(envdata$Plot,envdata$Time)) # gsub remove spaces in the combined variable
envdata$TimePlot <- paste(envdata$Time,"_",envdata$Plot)

envdata.DNA.f <- envdata[envdata$TimePlot %in% metadata.DNA.f$TimePlot,]
rownames(envdata.DNA.f)

envdata.RNA.f <- envdata[envdata$TimePlot %in% metadata.RNA.f$TimePlot,]
rownames(envdata.RNA.f)


#position of the levels of metadata in envdata 
match(envdata.DNA.f$TimePlot, metadata.DNA.f$TimePlot)

levels(metadata.DNA.f$TimePlot)
levels(factor(envdata.DNA.f$TimePlot))
dim(envdata.DNA.f)

match(envdata.RNA.f$TimePlot, metadata.RNA.f$TimePlot)

levels(metadata.RNA.f$TimePlot)
levels(factor(envdata.RNA.f$TimePlot))
dim(envdata.RNA.f)

###################
##### PHYLOSEQ #####
####################

#define the otutable
f2 = data.matrix(f) #need to be transformed to matrix of numeric/integers

# check the order of samples
order_idx2 <- match(rownames(metadata.f), colnames(f2))

# Reorder filtered data frame to have the sample names in the same order as the metadata data frame
f2  <- f2[ , order_idx2]

# VERIFY SAMPLE ORDER!!!
all(rownames(metadata.f) == colnames(f2)) #[1] TRUE

#match sample names of filtered asv table and metadata
colnames(f2) <- rownames(metadata.f)

OTU = otu_table(f2, taxa_are_rows = T)  
#very important: pay attention to the orientation of the table (taxa (OTUs) as rows or columns)


#define the taxonomy
tax2 = as.matrix(taxonomy)

TAX = tax_table(tax2)

#defining the map file
DAT <- sample_data(metadata.f)

#defining the phyloseq object
psq = phyloseq(OTU, TAX, DAT)
#Ok, done
psq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 986 taxa and 41 samples ]
# sample_data() Sample Data:       [ 41 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 986 taxa by 7 taxonomic ranks ]

# check the data 
ntaxa(psq)
nsamples(psq)
sample_names(psq)
rank_names(psq)
sample_variables(psq)
otu_table(psq)[1:5, 1:5]
tax_table(psq)[1:5, 1:5]
taxa_names(psq)[1:10]



# First lets remove out features with ambiguous phylum annotation.
psq.1 <- subset_taxa(psq, !is.na(Phylum) & !Phylum %in% c("Bacteria_unclassified", "Eukaryota_unclassified", "Archaea_unclassified", "Eukaryota","unknown"))
psq.1
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 979 taxa and 41 samples ]
# sample_data() Sample Data:       [ 41 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 979 taxa by 7 taxonomic ranks ]

# Transform data to dataframes
asv_f <- as.data.frame(otu_table(psq.1))
class(asv_f)
tax_f <- as.data.frame(tax_table(psq.1))
class(tax_f)

## Average replicates ## 
# TreatTimeType
# transpose
tasv_f <- as.data.frame(t(asv_f))
# add columns, merge, and delete
tasv_f$Treatment <- metadata.f$Treatment
tasv_f$Time <- metadata.f$Time
tasv_f$Type <- metadata.f$Type
tasv_f$TreatTimeType <- paste(tasv_f$Treatment,tasv_f$Time,tasv_f$Type)
tasv_f$Treatment <- NULL
tasv_f$Time <- NULL
tasv_f$Type <- NULL
# aggregate and delete col
agtasv_f<-aggregate(tasv_f, by=list(Treat_Time_Type=tasv_f$TreatTimeType), FUN =mean)
agtasv_f$TreatTimeType <- NULL
# col to rownames
agtasv_f_tbl <- column_to_rownames(as_tibble(agtasv_f), var="Treat_Time_Type")
# transpose again
asv_f_mean <- as.data.frame(t(agtasv_f_tbl))

# TreatType
# transpose
tasv_f <- as.data.frame(t(asv_f))
# add columns, merge, and delete
tasv_f$Treatment <- metadata.f$Treatment
tasv_f$Type <- metadata.f$Type
tasv_f$TreatType <- paste(tasv_f$Treatment,tasv_f$Type)
tasv_f$Treatment <- NULL
tasv_f$Type <- NULL
# aggregate and delete col
agtasv_f<-aggregate(tasv_f, by=list(Treat_Type=tasv_f$TreatType), FUN =mean)
agtasv_f$TreatType <- NULL
# col to rownames
agtasv_f_tbl <- column_to_rownames(as_tibble(agtasv_f), var="Treat_Type")
# transpose again
asv_f_meanTreatType <- as.data.frame(t(agtasv_f_tbl))


# check the order of samples. Unique levels cause we have 2 Dry T1 DNA and 2 Dry T1 cDNA (Rewetted T1 and Dry T1 are the same)
order_idx44 <- match(unique(metadata.f.av$TreatTimeType), colnames(asv_f_mean))

order_idx45 <- match(unique(metadata.f.av$TreatType), colnames(asv_f_meanTreatType))


# Reorder filtered data frame to have the sample names in the same order as the metadata data frame
asv_f_mean  <- asv_f_mean[,order_idx44]

asv_f_meanTreatType  <- asv_f_meanTreatType[,order_idx45]

# VERIFY SAMPLE ORDER!!!
all(unique(metadata.f.av$TreatTimeType) == colnames(asv_f_mean)) #[1] TRUE

all(unique(metadata.f.av$TreatType) == colnames(asv_f_meanTreatType)) #[1] TRUE


## DNA DATA ##
#define the otutable
f2.DNA = data.matrix(f.DNA) #need to be transformed to matrix of numeric/integers

# check the order of samples
order_idx3 <- match(rownames(metadata.DNA.f), colnames(f2.DNA))

# Reorder filtered data frame to have the sample names in the same order as the metadata data frame
f2.DNA  <- f2.DNA[ , order_idx3]

# VERIFY SAMPLE ORDER!!!
all(rownames(metadata.DNA.f) == colnames(f2.DNA)) #[1] TRUE

#match sample names of filtered asv table and metatada
colnames(f2.DNA) <- rownames(metadata.DNA.f)

OTU.DNA = otu_table(f2.DNA, taxa_are_rows = T)  
#very important: pay attention to the orientation of the table (taxa (OTUs) as rows or columns)

#define the taxonomy
tax2.DNA = as.matrix(taxonDNA)

TAX.DNA = tax_table(tax2.DNA)

#defining the map file
DAT.DNA<- sample_data(metadata.DNA.f)

#defining the phyloseq object
psq.DNA = phyloseq(OTU.DNA, TAX.DNA, DAT.DNA)
#Ok, done
psq.DNA
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 982 taxa and 29 samples ]
# sample_data() Sample Data:       [ 29 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 982 taxa by 7 taxonomic ranks ]

# First lets remove out features with ambiguous phylum annotation.
psq.1.DNA <- subset_taxa(psq.DNA, !is.na(Phylum) & !Phylum %in% c("Bacteria_unclassified", "Eukaryota_unclassified", "Archaea_unclassified", "Eukaryota","unknown"))
psq.1.DNA
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 975 taxa and 29 samples ]
# sample_data() Sample Data:       [ 29 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 975 taxa by 7 taxonomic ranks ]

# Transform data to dataframes
asv_f.DNA <- as.data.frame(otu_table(psq.1.DNA))
class(asv_f.DNA)
tax_f.DNA <- as.data.frame(tax_table(psq.1.DNA))
class(asv_f.DNA)


## Average replicates ##
# transpose
tasv_f.DNA <- as.data.frame(t(asv_f.DNA))
# add columns, merge, and delete
tasv_f.DNA$Treatment <- metadata.DNA.f$Treatment
tasv_f.DNA$Time <- metadata.DNA.f$Time
tasv_f.DNA$TreatTime <- paste(tasv_f.DNA$Treatment,tasv_f.DNA$Time)
tasv_f.DNA$Treatment <- NULL
tasv_f.DNA$Time <- NULL
# aggregate and delete col
agtasv_f.DNA<-aggregate(tasv_f.DNA, by=list(Treat_Time=tasv_f.DNA$TreatTime), FUN =mean)
agtasv_f.DNA$TreatTime <- NULL
# col to rownames
agtasv_f.DNA_tbl <- column_to_rownames(as_tibble(agtasv_f.DNA), var="Treat_Time")
# transpose again
asv_f.DNA_mean <- as.data.frame(t(agtasv_f.DNA_tbl))

# check the order of samples
order_idx4 <- match(unique(metadata.DNA.f.av$TreatTime), colnames(asv_f.DNA_mean))

# Reorder filtered data frame to have the sample names in the same order as the metadata data frame
asv_f.DNA_mean  <- asv_f.DNA_mean[,order_idx4]

# VERIFY SAMPLE ORDER!!!
all(unique(metadata.DNA.f.av$TreatTime) == colnames(asv_f.DNA_mean)) #[1] TRUE



## RNA DATA ##
#define the otutable
f2.RNA = data.matrix(f.RNA) #need to be transformed to matrix of numeric/integers

# check the order of samples
order_idx5 <- match(rownames(metadata.RNA.f), colnames(f2.RNA))

# Reorder filtered data frame to have the sample names in the same order as the metadata data frame
f2.RNA  <- f2.RNA[ , order_idx5]

# VERIFY SAMPLE ORDER!!!
all(rownames(metadata.RNA.f) == colnames(f2.RNA)) #[1] TRUE

#match sample names of filtered asv table and metadata
colnames(f2.RNA) <- rownames(metadata.RNA.f)

OTU.RNA = otu_table(f2.RNA, taxa_are_rows = T)  
#very important: pay attention to the orientation of the table (taxa (OTUs) as rows or columns)

#define the taxonomy
tax2.RNA = as.matrix(taxonRNA)

TAX.RNA = tax_table(tax2.RNA)

#defining the map file
DAT.RNA <- sample_data(metadata.RNA.f)

#defining the phyloseq object
psq.RNA = phyloseq(OTU.RNA, TAX.RNA, DAT.RNA)
#Ok, done
psq.RNA
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 804 taxa and 12 samples ]
# sample_data() Sample Data:       [ 12 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 804 taxa by 7 taxonomic ranks ]

# First lets remove out features with ambiguous phylum annotation.
psq.1.RNA <- subset_taxa(psq.RNA, !is.na(Phylum) & !Phylum %in% c("Bacteria_unclassified", "Eukaryota_unclassified", "Archaea_unclassified","Eukaryota","unknown"))
psq.1.RNA
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 798 taxa and 12 samples ]
# sample_data() Sample Data:       [ 12 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 798 taxa by 7 taxonomic ranks ]

# Transform data to dataframes
asv_f.RNA <- as.data.frame(otu_table(psq.1.RNA))
class(asv_f.RNA)
tax_f.RNA <- as.data.frame(tax_table(psq.1.RNA))
class(tax_f.RNA)

## Average replicates ##
# transpose
tasv_f.RNA <- as.data.frame(t(asv_f.RNA))
# add columns, merge, and delete
tasv_f.RNA$Treatment <- metadata.RNA.f$Treatment
tasv_f.RNA$Time <- metadata.RNA.f$Time
tasv_f.RNA$TreatTime <- paste(tasv_f.RNA$Treatment,tasv_f.RNA$Time)
tasv_f.RNA$Treatment <- NULL
tasv_f.RNA$Time <- NULL
# aggregate and delete col
agtasv_f.RNA<-aggregate(tasv_f.RNA, by=list(Treat_Time=tasv_f.RNA$TreatTime), FUN =mean)
agtasv_f.RNA$TreatTime <- NULL
# col to rownames
agtasv_f.RNA_tbl <- column_to_rownames(as_tibble(agtasv_f.RNA), var="Treat_Time")
# transpose again
asv_f.RNA_mean <- as.data.frame(t(agtasv_f.RNA_tbl))

# check the order of samples
order_idx6 <- match(unique(metadata.RNA.f.av$TreatTime), colnames(asv_f.RNA_mean))

# Reorder filtered data frame to have the sample names in the same order as the metadata data frame
asv_f.RNA_mean  <- asv_f.RNA_mean[,order_idx6]

# VERIFY SAMPLE ORDER!!!
all(unique(metadata.RNA.f.av$TreatTime) == colnames(asv_f.RNA_mean)) #[1] TRUE


#########################
### TAXONOMY BARPLOTS ###
#########################

taxasums_ord<-summarize_taxa(asv_f, tax_f)$Order


#average
taxasums_ord.av<-summarize_taxa(asv_f_mean, tax_f)$Order

#reorder cols
taxasums_ord.av<-taxasums_ord.av[,c(2,6,8,4,11,13,9,1,5,7,3,10,12)]

taxasums_ord.av2<-summarize_taxa(asv_f_meanTreatType, tax_f)$Order

#reorder cols
taxasums_ord.av2<-taxasums_ord.av2[,c(2,4,1,3)]

colorPhy.av = c("#A6CEE3","gray","yellow3","coral","cyan","green","pink","brown","orange","orange4","aquamarine","black","yellow","blue","red","goldenrod")
colorPhy.av2 = c("#A6CEE3","gray","yellow3","coral","cyan","green","pink","brown","orange","aquamarine","orange4","black","yellow","blue","red","goldenrod")


# create larger palettes by interpolating existing ones with constructor function colorRampPalette
# then getPalette(num of colors you want)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

bp3 <- taxa_barplot(taxasums_ord, metadata.f, "Treatment", ntoplot = 15)
bp3 + ylab("Relative Abundance (%)") + 
  scale_x_discrete(breaks=rownames(metadata.f),labels=metadata.f$Sample) + 
  scale_fill_manual(values = getPalette(16)) +
  theme(legend.position="bottom", text = element_text(size = 8),
        legend.title = element_text(size = 14, color="black", face = "bold"), legend.text = element_text(size=12),
        axis.title = element_text(size=14, color="black", face = "bold"), axis.text = element_text(size=12,color="black"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  )
ggsave("barplot_rewet_SILVA138_order.png", height=10, width=26, device="png", dpi=300) 


getPalette(16)
colorOrd.av = c("#A6CEE3","#63A8A0","#438EC0","cyan","#3BA432","#F16667","#B89B74","#E62F27","#FE982C","#F9A963","#ED8F47","#C3AAD2", "#7D54A5", "#B9A499", "#EAD27A","#B15928")

bp3.av <- taxa_barplot(taxasums_ord.av, metadata.f.av, "Treatment", ntoplot = 15)
#reorder samplesID in taxa_barplot object 
bp3.av$data$SampleID <- factor(bp3.av$data$SampleID,levels = c("Dry T1 DNA","Dry T3 DNA","Dry T7 DNA","Dry T11 DNA","Rewetted T3 DNA","Rewetted T7 DNA","Rewetted T11 DNA",
                                                               "Dry T1 cDNA","Dry T3 cDNA","Dry T7 cDNA","Dry T11 cDNA","Rewetted T3 cDNA","Rewetted T7 cDNA"))

bp3.av + ylab("Relative Abundance (%)") + 
  scale_x_discrete(breaks=unique(metadata.f.av$TreatTimeType),labels=unique(metadata.f.av$TreatTimeType)) + 
  scale_fill_manual(values = colorOrd.av) +
  theme(legend.position="bottom", text = element_text(size = 8),
        legend.title = element_text(size = 14, color="black", face = "bold"), legend.text = element_text(size=12),
        axis.title = element_text(size=14, color="black", face = "bold"), axis.text = element_text(size=12,color="black"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  )
ggsave("barplot_Rewet_SILVA138_order_av.svg", height=10, width=25, device="svg", dpi=300)

bp3.av2 <- taxa_barplot(taxasums_ord.av2, metadata.f.av, "Treatment", ntoplot = 15)
#reorder samplesID in taxa_barplot object 
bp3.av2$data$SampleID <- factor(bp3.av2$data$SampleID,levels = c("Dry DNA","Rewetted DNA","Dry cDNA","Rewetted cDNA"))

bp3.av2 + ylab("Relative Abundance (%)") + 
  scale_x_discrete(breaks=unique(metadata.f.av$TreatType),labels=unique(metadata.f.av$TreatType)) + 
  scale_fill_manual(values = getPalette(16)) +
  theme(legend.position="bottom", text = element_text(size = 8),
        legend.title = element_text(size = 12, color="black", face = "bold"), legend.text = element_text(size=10),
        axis.title = element_text(size=12, color="black", face = "bold"), axis.text = element_text(size=12,color="black"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  )
ggsave("barplot_Rewet_SILVA138_order_avDryRew.svg", height=10, width=20, device="svg", dpi=300)




## DNA DATA ##
taxasums_ord.DNA<-summarize_taxa(asv_f.DNA, tax_f.DNA)$Order


#average
taxasums_ord.DNA.av<-summarize_taxa(asv_f.DNA_mean, tax_f.DNA)$Order

#reorder cols
taxasums_ord.DNA.av<-taxasums_ord.DNA.av[,c(1,3,4,2,6,7,5)]


bp3.DNA <- taxa_barplot(taxasums_ord.DNA, metadata.DNA.f, "Treatment", ntoplot = 15)
bp3.DNA + ylab("Relative Abundance (%)") + 
          scale_x_discrete(breaks=rownames(metadata.DNA.f)) + 
          scale_fill_manual(values = getPalette(16)) +
          theme(legend.position="bottom", text = element_text(size = 8),
                legend.title = element_text(size = 14, color="black", face = "bold"), legend.text = element_text(size=12),
                axis.title = element_text(size=14, color="black", face = "bold"), axis.text = element_text(size=12,color="black"),
                panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  )
ggsave("barplot_Rewet_SILVA138_order_DNA.png", height=10, width=26, device="png", dpi=300) 


bp3.DNA.av <- taxa_barplot(taxasums_ord.DNA.av, metadata.DNA.f.av, "Treatment", ntoplot = 15)
#reorder samplesID in taxa_barplot object 
bp3.DNA.av$data$SampleID <- factor(bp3.DNA.av$data$SampleID,levels = c("Dry T1","Dry T3","Dry T7","Dry T11","Rewetted T1","Rewetted T3","Rewetted T7","Rewetted T11"))

bp3.DNA.av + ylab("Relative Abundance (%)") + 
  scale_x_discrete(breaks=rownames(metadata.DNA.f.av),labels=metadata.DNA.f.av$Sample) + 
  scale_fill_manual(values = getPalette(16)) +
  theme(legend.position="bottom", text = element_text(size = 8),
        legend.title = element_text(size = 14, color="black", face = "bold"), legend.text = element_text(size=12),
        axis.title = element_text(size=14, color="black", face = "bold"), axis.text = element_text(size=12,color="black"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  )
ggsave("barplot_Rewet_SILVA138_order_DNA_av.svg", height=10, width=25, device="svg", dpi=300) 



## RNA DATA ##
taxasums_ord.RNA<-summarize_taxa(asv_f.RNA, tax_f.RNA)$Order


#average
taxasums_ord.RNA.av<-summarize_taxa(asv_f.RNA_mean, tax_f.RNA)$Order

#reorder cols
taxasums_ord.RNA.av<-taxasums_ord.RNA.av[,c(1,3,4,2,5,6)]


bp3.RNA <- taxa_barplot(taxasums_ord.RNA, metadata.RNA.f, "Treatment", ntoplot = 15)
bp3.RNA + ylab("Relative Abundance (%)") + 
          scale_x_discrete(breaks=rownames(metadata.RNA.f)) + 
          scale_fill_manual(values = getPalette(16)) +
          theme(legend.position="bottom", text = element_text(size = 8),
                legend.title = element_text(size = 14, color="black", face = "bold"), legend.text = element_text(size=12),
                axis.title = element_text(size=14, color="black", face = "bold"), axis.text = element_text(size=12,color="black"),
                panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, size=1.5)
          )
ggsave("barplot_Rewet_SILVA138_order_RNA.png", height=10, width=26, device="png", dpi=300) 


bp3.RNA.av <- taxa_barplot(taxasums_ord.RNA.av, metadata.RNA.f.av, "Treatment", ntoplot = 15)
#reorder samplesID in taxa_barplot object 
bp3.RNA.av$data$SampleID <- factor(bp3.RNA.av$data$SampleID,levels = c("Dry T1","Dry T3","Dry T7","Dry T11","Rewetted T1","Rewetted T3","Rewetted T7"))

bp3.RNA.av + ylab("Relative Abundance (%)") + 
  scale_x_discrete(breaks=rownames(metadata.RNA.f.av),labels=metadata.RNA.f.av$Sample) + 
  scale_fill_manual(values = getPalette(16)) +
  theme(legend.position="bottom", text = element_text(size = 8),
        legend.title = element_text(size = 14, color="black", face = "bold"), legend.text = element_text(size=12),
        axis.title = element_text(size=14, color="black", face = "bold"), axis.text = element_text(size=12,color="black"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  )
ggsave("barplot_Rewet_SILVA138_order_RNA_av2.svg", height=10, width=24, device="svg", dpi=300) 

bp4.RNA.av <- taxa_barplot(taxasums_gen.RNA.av, metadata.RNA.f.av, "Treatment", ntoplot = 15)
#reorder samplesID in taxa_barplot object 
bp4.RNA.av$data$SampleID <- factor(bp4.RNA.av$data$SampleID,levels = c("Dry T1","Dry T3","Dry T7","Dry T11","Rewetted T1","Rewetted T3","Rewetted T7"))

bp4.RNA.av + ylab("Relative Abundance (%)") + 
  scale_x_discrete(breaks=rownames(metadata.RNA.f.av),labels=metadata.RNA.f.av$Sample) + 
  scale_fill_manual(values = getPalette(16)) +
  theme(legend.position="bottom", text = element_text(size = 8),
        legend.title = element_text(size = 14, color="black", face = "bold"), legend.text = element_text(size=12),
        axis.title = element_text(size=14, color="black", face = "bold"), axis.text = element_text(size=12,color="black"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  )
ggsave("barplot_Rewet_SILVA138_genus_RNA_av.png", height=14, width=38, device="png", dpi=300) 



# Percentage of each group by sample
ps <- as.data.frame(t(taxasums_ord)/rowSums(t(taxasums_ord))*100)
ps.t <- as.data.frame(t(ps))
ps.ord <- ps.t[order(-ps.t$D1T3),]
ps.ord[1:5,]

#write.table(ps.ord,"percentage_Rewet_ord.txt", sep="\t")

# Percentage of each group by sample
ps.DNA <- as.data.frame(t(taxasums_ord.DNA)/rowSums(t(taxasums_ord.DNA))*100)
ps.t.DNA <- as.data.frame(t(ps.DNA))
ps.ord.DNA <- ps.t.DNA[order(-ps.t.DNA$D1T3),]
ps.ord.DNA[1:5,]

#write.table(ps.ord.DNA,"percentage_Rewet_ord_DNA.txt", sep="\t")

# Percentage of each group by sample
ps.RNA <- as.data.frame(t(taxasums_ord.RNA)/rowSums(t(taxasums_ord.RNA))*100)
ps.t.RNA <- as.data.frame(t(ps.RNA))
ps.ord.RNA <- ps.t.RNA[order(-ps.t.RNA$cD3T1),]
ps.ord.RNA[1:5,]

#write.table(ps.ord.RNA,"percentage_Rewet_ord_RNA.txt", sep="\t")


# Percentage of each group by sample
ps.av <- as.data.frame(t(taxasums_ord.av)/rowSums(t(taxasums_ord.av))*100)
ps.t.av <- as.data.frame(t(ps.av))
ps.ord.av <- ps.t.av[order(-ps.t.av$`Dry T1 DNA`),]
ps.ord.av[1:5,]

#write.table(ps.ord.av,"percentage_Rewet_ord_av.txt", sep="\t")

# Percentage of each group by sample
ps.av2 <- as.data.frame(t(taxasums_ord.av2)/rowSums(t(taxasums_ord.av2))*100)
ps.t.av2 <- as.data.frame(t(ps.av2))
ps.ord.av2 <- ps.t.av2[order(-ps.t.av2$`Dry DNA`),]
ps.ord.av2[1:5,]

#write.table(ps.ord.av2,"percentage_Rewet_ord_avDvsR.txt", sep="\t")



#######################
### ALPHA DIVERSITY ###
#######################
# BLOCK2
# Diversity indices and ANOVAs (see below)
# data should be as samples x features (MARGIN=1): transpose filtered asvtable (f)
# input data: filtered asvtable 0.1% (f)

## SHANNON ##
H <- vegan::diversity(t(f), index="shannon", MARGIN = 1)
dfH <- as.data.frame(H)

#distribution
fit.gama=fitdist(dfH$H,"gamma")
fit.lognormal=fitdist(dfH$H,"lnorm")
fit.normal=fitdist(dfH$H,"norm")
par(mfrow=c(1,2))
cdfcomp (list(fit.gama,fit.lognormal,fit.normal))
qqcomp(list(fit.gama,fit.lognormal,fit.normal))
par(mfrow=c(1,1))
#Based on the graphical criterion, normal dist fits the best

gofstat(list(fit.lognormal,fit.gama,fit.normal), fitnames=c("Lognormal", "Gama","Normal"))
#Based on the statistics, (the lower AKAIKE value the better) normal dist fits the best
# Goodness-of-fit criteria
#                               Lognormal     Gama   Normal
# Akaike's Information Criterion  95.08256 93.10177 89.58179
# Bayesian Information Criterion  98.50970 96.52892 93.00893

# Boxplot
dfH$Type <- as.factor(metadata.f$Type)
dfH$Treatment <- as.factor(metadata.f$Treatment)
dfH$Time <- as.factor(metadata.f$Time)
dfH$TreatTime<- as.factor(metadata.f$TreatTime)

# Filter out RNA samples
dfH.DNA <- dfH[dfH$Type == "DNA",] 

dfH.DNA.summarized = dfH.DNA %>% group_by(Treatment) %>% summarize(Max.H=max(H)) 
dfH.DNA.summarized = dfH.DNA.summarized[order(factor(dfH.DNA.summarized$Treatment, levels = c("Dry","Rewetted"))),]

  
pH <- ggplot(dfH.DNA, aes(x=Treatment, y=H, color=Treatment)) + 
          geom_boxplot() +
          geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
          ylab("Shannon Diversity index") +
          xlab("Treatment") +
          scale_colour_manual(values=c("coral1","cyan3")) +
          theme(legend.position = "none", axis.text.x = element_text(hjust = 0.5,vjust= 0.5),
                panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, linewidth = 1.5), 
                axis.text=element_text(size=16, color="black"), 
                axis.title = element_text(size=18, color="black", face = "bold"),
                panel.grid =element_blank()) +
                geom_text(data=dfH.DNA.summarized, 
                    aes(x=Treatment,y=0.1+Max.H, label=hsd$groups[order(factor(rownames(hsd$groups), 
                                                                                  levels = c("Dry","Rewetted"))),2], vjust=0))
  
ggsave("Shannon_Treatment_DNA.png", height=4, width=6, device="png") # save a eps 4 inches by 6 inches


# Filter out DNA samples
dfH.RNA <- dfH[dfH$Type == "cDNA",] 

dfH.RNA.summarized = dfH.RNA %>% group_by(Treatment) %>% summarize(Max.H=max(H)) 
dfH.RNA.summarized = dfH.RNA.summarized[order(factor(dfH.RNA.summarized$Treatment, levels = c("Dry","Rewetted"))),]


pH2 <- ggplot(dfH.RNA, aes(x=Treatment, y=H, color=Treatment)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  ylab("Shannon Diversity index") +
  xlab("Treatment") +
  scale_colour_manual(values=c("coral1","cyan3")) +
  theme(legend.position = "none", axis.text.x = element_text(hjust = 0.5,vjust= 0.5),
        panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, linewidth = 1.5), 
        axis.text=element_text(size=16, color="black"), 
        axis.title = element_text(size=18, color="black", face = "bold"),
        panel.grid =element_blank()) +
  geom_text(data=dfH.DNA.summarized, 
            aes(x=Treatment,y=0.1+Max.H, label=wilcoxlet$Signif.group,vjust=0))

ggsave("Shannon_Treatment_RNA.png", height=4, width=6, device="png") # save a eps 4 inches by 6 inches


## SPECIES RICHNESS ##
S <- specnumber(t(f))
dfS <- as.data.frame(S)

# Boxplot
dfS$Type <- as.factor(metadata.f$Type)
dfS$Treatment <- as.factor(metadata.f$Treatment)
dfS$Time <- as.factor(metadata.f$Time)
dfS$TreatTime<- as.factor(metadata.f$TreatTime)

# Filter out RNA samples
dfS.DNA <- dfS[dfS$Type == "DNA",] 

dfS.DNA.summarized = dfS.DNA %>% group_by(Treatment) %>% summarize(Max.S=max(S)) 
dfS.DNA.summarized = dfS.DNA.summarized[order(factor(dfS.DNA.summarized$Treatment, levels = c("Dry","Rewetted"))),]


pS <- ggplot(dfS.DNA, aes(x=Treatment, y=S, color=Treatment)) + 
            geom_boxplot() +
            geom_jitter(shape=16, position=position_jitter(0.2), aes(color=Treatment), size=2) +
            ylab("Species richness") +
            xlab("Treatment") +
            scale_colour_manual(values=c("coral1","cyan3")) +
            theme(legend.position = "none", axis.text.x = element_text(hjust = 0.5,vjust= 0.5),
                  panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, linewidth = 1.5), 
                  axis.text=element_text(size=16, color="black"), 
                  axis.title = element_text(size=18, color="black", face = "bold"),
                  panel.grid =element_blank()) +
            geom_text(data=dfS.DNA.summarized, 
                  aes(x=Treatment,y=50+Max.S, label=hsd2$groups[order(factor(rownames(hsd2$groups), 
                                                                       levels = c("Dry","Rewetted"))),2], vjust=0))

            
ggsave("Richness_DNA_Treatment.png", height=4, width=6, device="png") # save a eps 4 inches by 6 inches


# Filter out DNA samples
dfS.RNA <- dfS[dfS$Type == "cDNA",] 

dfS.RNA.summarized = dfS.RNA %>% group_by(Treatment) %>% summarize(Max.S=max(S)) 
dfS.RNA.summarized = dfS.RNA.summarized[order(factor(dfS.RNA.summarized$Treatment, levels = c("Dry","Rewetted"))),]


pS2 <- ggplot(dfS.RNA, aes(x=Treatment, y=S, color=Treatment)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), aes(color=Treatment), size=2) +
  ylab("Species richness") +
  xlab("Treatment") +
  scale_colour_manual(values=c("coral1","cyan3")) +
  theme(legend.position = "none", axis.text.x = element_text(hjust = 0.5,vjust= 0.5),
        panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, linewidth = 1.5), 
        axis.text=element_text(size=16, color="black"), 
        axis.title = element_text(size=18, color="black", face = "bold"),
        panel.grid =element_blank()) +
  geom_text(data=dfS.RNA.summarized, 
            aes(x=Treatment,y=50+Max.S, label=wilcoxlet$Signif.group,vjust=0))


ggsave("Richness_RNA_Treatment.png", height=4, width=6, device="png") # save a eps 4 inches by 6 inches




# multiple plots in one page 
library(ggpubr)
ggarrange(pH , pS, #+ rremove("x.text")
          labels = c("A", "B"), label.x = c(0, 0), label.y = c(1, 1),
          align = "hv", # align plot vertically and horiz.
          font.label = list(size = 20, color = "black", face = "bold"),
          common.legend = FALSE, legend = NULL,
          ncol = 2, nrow = 1)
ggsave("AlphaDiversity_Rewet.svg", height=5, width=12, device="svg", dpi = 300)

ggarrange(pH2 , pS2, #+ rremove("x.text")
          labels = c("C", "D"), label.x = c(0, 0), label.y = c(1, 1),
          align = "hv", # align plot vertically and horiz.
          font.label = list(size = 20, color = "black", face = "bold"),
          common.legend = FALSE, legend = NULL,
          ncol = 2, nrow = 1)
ggsave("AlphaDiversity_RewetRNA.svg", height=5, width=12, device="svg", dpi = 300)


# ANOVAs
# use the Anova(mod, type) function to conduct the Type III SS ANOVA
# Type-III tests do not assume that terms higher-order to the term in
# > question are zero. For example, in a two-way design with factors A and B,
# > the type-III test for the A main effect tests whether the population
# > marginal means at the levels of A (i.e., averaged across the levels of B)
# > are the same. One can test this hypothesis whether or not A and B interact,
# > since the marginal means can be formed whether or not the profiles of means
# > for A within levels of B are parallel. Whether the hypothesis is of interest
# > in the presence of interaction is another matter, however. To compute
# > Type-III tests using incremental F-tests, one needs contrasts that are
# > orthogonal in the row-basis of the model matrix. In R, this means, e.g.,
# > using contr.sum, contr.helmert, or contr.poly (all of which will give you
# > the same SS), but not contr.treatment. Failing to be careful here will
# > result in testing hypotheses that are not reasonably construed, e.g., as
# > hypotheses concerning main effects.

##lme()
library(nlme)
##lmer()
library(lme4)
## Anova(), levenTest()
library(car)
#glm
library(glmmTMB)

## Multiple comparisons
library(multcomp)
library(emmeans)
library(multcompView)
library(agricolae)
library(RVAideMemoire)

## DNA ##
# DEPENDENT VARIABLE SHANNON
mod <- lm(H ~ Treatment*Time, data=dfH.DNA)
summary(mod)
Anova(mod, type = 3)

# Colinearity. Remove Time predictor (model with Time not significant)

mod.1 <- lm(H ~ Treatment, data=dfH.DNA)
summary(mod.1)
Anova(mod.1, type = 3, white.adjust="hc3")
# Analysis of Deviance Table (Type III tests)
# 
# Response: H
# Df       F    Pr(>F)    
# (Intercept)  1 592.624 < 2.2e-16 ***
#   Treatment    1  24.813 3.203e-05 ***
#   Residuals   27 

# Treatment significant

# Tukey test to study each pair of treatment :
hsd=HSD.test(mod.1, "Treatment", group=T)
# $groups
# H groups
# Rewetted 5.448507      a
# Dry      4.446210      b


# DEPENDENT VARIABLE RICHNESS
mod2 <- lm(S ~ Treatment*Time, data=dfS.DNA)
summary(mod2)
Anova(mod2, type = 3) 

# Colinearity. Remove Time predictor (model with Time not significant)

mod2.1 <- lm(S ~ Treatment, data=dfS.DNA)
summary(mod2.1)
Anova(mod2.1, type = 3) 

# Anova Table (Type III tests)
# 
# Response: S
# Sum Sq Df F value    Pr(>F)    
# (Intercept) 2125907  1   94.29 2.658e-10 ***
#   Treatment    231552  1   10.27  0.003459 ** 
#   Residuals    608753 27 

# Treatment significant

# Tukey test to study each pair of treatment :
hsd2=HSD.test(mod2.1, "Treatment", group=T)
# $groups
#           S (mean)        groups
# Rewetted 555.2857      a
# Dry      376.4667      b


## RNA ##
# DEPENDENT VARIABLE SHANNON
mod4 <- lm(H ~ Treatment*Time, data=dfH.RNA)
summary(mod4)
Anova(mod4, type = 3)

# Colinearity. Remove Time predictor (model with Time not significant)

# Not normally distributed data: Apply Kruskal-wallis test non-parametric
Kruskal <- kruskal.test(H ~ Treatment, data = dfH.RNA)
# Kruskal-Wallis rank sum test
# 
# data:  H by Treatment
# Kruskal-Wallis chi-squared = 0.53407, df = 1, p-value = 0.4649

# create a dataframe with letters meaning significance among groups based on pairwise wilcox test 
wilcoxlet <- data.frame("Treatmenet" =as.factor(c("Dry","Rewetted")), "Signif.group"= c("a","a"))

# DEPENDENT VARIABLE RICHNESS
mod5 <- lm(S ~ Treatment*Time, data=dfS.RNA)
summary(mod5)
Anova(mod5, type = 3) 

# Colinearity. Remove Time predictor (model with Time not significant)

mod5.1 <- lm(S ~ Treatment, data=dfS.RNA)
summary(mod5.1)
Anova(mod5.1, type = 3) 

# Not normally distributed data: Apply Kruskal-wallis test non-parametric
Kruskal2 <- kruskal.test(S ~ Treatment, data = dfS.RNA)
# Kruskal-Wallis rank sum test
# 
# data:  S by Treatment
# Kruskal-Wallis chi-squared = 0.0065934, df = 1, p-value = 0.9353


# create a dataframe with letters meaning significance among groups based on pairwise wilcox test 
wilcoxlet2 <- data.frame("Treatmenet" =as.factor(c("Dry","Rewetted")), "Signif.group"= c("a","a"))


### ASSUMPTIONS ###

# Normality test
residuals<-resid(mod)
hist(residuals)
shapiro.test(residuals) #Ok normal shannon DNA

residuals.1<-resid(mod.1)
hist(residuals.1)
shapiro.test(residuals.1) #Ok normal shannon DNA

residuals2.1<-resid(mod2.1)
hist(residuals2.1)
shapiro.test(residuals2.1) #Ok Normal Richness DNA 

residuals4.1<-resid(mod4.1)
hist(residuals4.1)
shapiro.test(residuals4.1) #Not Normal shannon DNA 

residuals5.1<-resid(mod5.1)
hist(residuals5.1)
shapiro.test(residuals5.1) #Not Normal richness DNA 


# Homogeneity of variances
leveneTest(H ~ Treatment, data=dfH.DNA) # heterocedasticity shannon.
leveneTest(S ~ Treatment, data=dfS.DNA) # homocedasticity richness



# if there's no homocedasticity you can correct for that in the Anova() by using the argument white.adjust="hc3"



## RAREFACTION CURVES
# Extract the min read counts of a sample
raremax <- min(rowSums(t(f))) # [1] 12938
raremax.DNA <- min(rowSums(t(f.DNA))) # [1] 12938
raremax.RNA <- min(rowSums(t(f.RNA))) # [1] 15398

# Rarefy species richness at that sampling depth
Srare <- rarefy(t(f), raremax)
Srare.DNA <- rarefy(t(f.DNA), raremax.DNA)
Srare.RNA <- rarefy(t(f.RNA), raremax.RNA)

# plot btw observed and rarefied ASV richness 
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
# almost fits in a y=x equation

# parameters for plotting. Combination of colors, line types and thickness
col <- rainbow(15)
lty <- c("solid", "dashed", "dotdash")
lwd <- c(1,2)
pars <- expand.grid(col = col, lty = lty, lwd = lwd, stringsAsFactors = FALSE)
head(pars)

# plot rarefaction curves of each sample
out <- with(pars[1:87, ],
            rarecurve(t(f), step = 20, sample = raremax, col = col,
                      lty = lty, lwd = lwd, label = FALSE))

out.DNA <- with(pars[1:87, ],
                rarecurve(t(f.DNA), step = 20, sample = raremax.DNA, col = col,
                          lty = lty, lwd = lwd, label = FALSE))

out.RNA <- with(pars[1:87, ],
                rarecurve(t(f.RNA), step = 20, sample = raremax.RNA, col = col,
                          lty = lty, lwd = lwd, label = FALSE))

# plot rarecurves as a fucntion of sample 
Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
Smax <- sapply(out, max)
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Size",
     ylab = "Species", type = "n")
abline(v = raremax)
for (i in seq_along(out)) {
  N <- attr(out[[i]], "Subsample")
  with(pars, lines(N, out[[i]], col = col[i], lty = lty[i], lwd = lwd[i]))
}

Nmax.DNA <- sapply(out.DNA, function(x) max(attr(x, "Subsample")))
Smax.DNA <- sapply(out.DNA, max)
plot(c(1, max(Nmax.DNA)), c(1, max(Smax.DNA)), xlab = "Sample Size",
     ylab = "Species", type = "n")
abline(v = raremax.DNA)
for (i in seq_along(out.DNA)) {
  N <- attr(out.DNA[[i]], "Subsample")
  with(pars, lines(N, out.DNA[[i]], col = col[i], lty = lty[i], lwd = lwd[i]))
}

Nmax.RNA <- sapply(out.RNA, function(x) max(attr(x, "Subsample")))
Smax.RNA <- sapply(out.RNA, max)

# plot rarecurves as a function of a categorical variable Type
col <- c('blue', 'red')
set.seed(3)
grp <- factor(metadata.f$Type) 
cols <- col[grp]

plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Size",
     ylab = "ASV", type = "n")
legend("bottomright", legend=c("DNA","RNA"), title ="Type",
       col = col, lwd=2, seg.len = 0.5, cex = 0.65)
#abline(v = raremax)
for (i in seq_along(out)) {
  N <- attr(out[[i]], "Subsample")
  lines(N, out[[i]], col = cols[i], lwd = 2)
}

# plot rarecurves as a function of a categorical variable Treatment
col <- c('red','cyan')
set.seed(3)
grp <- factor(metadata.DNA.f$Treatment) 
cols <- col[grp]

plot(c(1, max(Nmax.DNA)), c(1, max(Smax.DNA)), xlab = "Sample Size",
     ylab = "ASV", type = "n")
legend("bottomright", legend=c("Dry","Rewetted"), title ="Treatment",
       col = col, lwd=2, seg.len = 0.5, cex = 0.65)
#abline(v = raremax.DNA)
for (i in seq_along(out.DNA)) {
  N <- attr(out.DNA[[i]], "Subsample")
  lines(N, out.DNA[[i]], col = cols[i], lwd = 2)
}

# plot rarecurves as a function of a categorical variable Treatment
col <- c('red','cyan')
set.seed(3)
grp <- factor(metadata.RNA.f$Treatment) 
cols <- col[grp]

plot(c(1, max(Nmax.RNA)), c(1, max(Smax.RNA)), xlab = "Sample Size",
     ylab = "ASV", type = "n")
legend("bottomright", legend=c("Dry","Rewetted"), title ="Treatment",
       col = col, lwd=2, seg.len = 0.5, cex = 0.65)
#abline(v = raremax.RNA)
for (i in seq_along(out.RNA)) {
  N <- attr(out.RNA[[i]], "Subsample")
  lines(N, out.RNA[[i]], col = cols[i], lwd = 2)
}


########################
#### BETA DIVERSITY ####
########################

##### PCA #####
# Compositional biplot generated from the CLR transformed ASV values with 0 replacement. 
# The biplot is shown with scale=1, projecting the distances between samples onto the 2-D image as best as possible. 
# The distance between ASVs (red) is a measure of their compositional association, but this is only a crude measure relative to rho

# replace 0 values with an estimate: Count Zero Multiplicative method.
# with the default method Geometric Bayesian Multiplicative (GBM) I got an error probably due to taxa with a high percentage of zeros; 
# Unique label (numeric or character) used to denote count zeros in X (default label=0);
# output dataset in proportions (as default)
# according to cmultRepl the input should be samples x features (transpose)
f.n0 <- cmultRepl(t(f), method="CZM", label=0, output = "prop", z.warning=0.999)
dim(f.n0) # [1]  41 986

f.n0.DNA <- cmultRepl(t(f.DNA), method="CZM", label=0, output = "prop", z.warning=0.999)
dim(f.n0.DNA) # [1]  29 982

# Generate the centered log-ratio transformed data for plotting later
# samples by row (1)
# use for alpha diversity, ordination and dendrogram
f.clr <- apply(f.n0, 1, function(x) log(x) - mean(log(x)))
dim(f.clr) # [1] 986  41

f.clr.DNA <- apply(f.n0.DNA, 1, function(x) log(x) - mean(log(x)))
dim(f.clr.DNA) # [1] 982  29

## ALL DATA ##
# perform a singular value decomposition (SVD)
# samples x features (transpose f.clr)
pcx <- prcomp(t(f.clr))
summary(pcx)
names(pcx)
pcx$x[1:5,1:3] # PCA values of samples
pcx$rotation[1:5,1:3] # PCA loadings of variables (ASVs)

# dataframe of PC values from samples creating additional cols 
# check the order of samples
order_idx50 <- match(rownames(metadata.f), rownames(pcx$x))

# Reorder filtered data frame to have the sample names in the same order as the metadata data frame
pcx$x  <- pcx$x[order_idx50,]

# VERIFY SAMPLE ORDER!!!
all(rownames(metadata.f) == rownames(pcx$x)) #[1] TRUE

PCvalues <- data.frame(Sample = rownames(metadata.f), Type = metadata.f$Type, Treatment = metadata.f$Treatment, Time = metadata.f$Time, TreatTime = metadata.f$TreatTime, pcx$x)
class(PCvalues$Type)

# variance explained by each principal component
var_explained <- pcx$sdev^2/sum(pcx$sdev^2)
var_explained[1:5]

## DNA DATA ##
# perform a singular value decomposition (SVD)
# samples x features (transpose f.clr)
pcx.DNA <- prcomp(t(f.clr.DNA))
summary(pcx.DNA)
names(pcx.DNA)
pcx.DNA$x[1:5,1:3] # PCA values of samples
pcx.DNA$rotation[1:5,1:3] # PCA loadings of variables (ASVs)

# dataframe of PC values from samples creating additional cols 
# check the order of samples
order_idx51 <- match(rownames(metadata.DNA.f), rownames(pcx.DNA$x))

# Reorder filtered data frame to have the sample names in the same order as the metadata data frame
pcx.DNA$x  <- pcx.DNA$x[order_idx51,]

# VERIFY SAMPLE ORDER!!!
all(rownames(metadata.DNA.f) == rownames(pcx.DNA$x)) #[1] TRUE

PCvalues.DNA <- data.frame(Sample = rownames(metadata.DNA.f), Treatment = metadata.DNA.f$Treatment, Time = metadata.DNA.f$Time, TreatTime = metadata.DNA.f$TreatTime, pcx.DNA$x)
class(PCvalues.DNA$Treatment)

# variance explained by each principal component
var_explained.DNA <- pcx.DNA$sdev^2/sum(pcx.DNA$sdev^2)
var_explained.DNA[1:5]


# BLOCK3
# Plot of the relationship between the ASVs from the SVD. Here each ASV is colored red if it was identified as having an effect size greater than 1. 
# When we plot the individual ASVs separately as in `R_Block_3` on a loadings plot 
# the resulting plot contains qualitative information on both differential relative abundance of ASVs and 
# on association between ASVs. The distance and direction of ASV from the center of the plot is proportional 
# to the standard deviation of the CLR value of that ASV in the dataset. 
# In theory, ASVs that have a short link between them, (i.e., are nearly the same direction and distance from the origin) 
# should be compositionally associated.

# calculation of SVD in R_Block_2
# plot(pcx$rotation[,1], pcx$rotation[,2],
#      main="loadings", xlab=xlab, cex=0.5,
#      ylab=ylab, col=rgb(0,0,0,0.3), pch=19)
# points(pcx$rotation[,1][high.e],
#        pcx$rotation[,2][high.e], cex=0.8,
#        col=rgb(1,0,0,0.5))
# points(pcx$rotation[c("39306","39235"),1],
#        pcx$rotation[c("39306","39235"),2],
#        cex=0.5, col=rgb(0,0,1,1), pch=19)

pcx$rotation %>% 
  as.data.frame %>%
  rownames_to_column("ASV") %>%
  ggplot(aes(x=PC1,y=PC2)) + 
  geom_point(size=1, alpha=0.3, show.legend = FALSE) +
  #geom_point(data=PC_h.e.LF, aes(x=PC1,y=PC2), size=3, alpha=0.3, color= 'blue1', show.legend = FALSE) +
  #geom_point(data=PC_h.e.CC, aes(x=PC1,y=PC2), size=3, alpha=0.3, color= 'brown', show.legend = FALSE) +
  geom_segment(data = PCvalues, aes(x = 0, y = 0, xend = (PC1/500), yend = (PC2/500), color = Treatment), 
               alpha=0.7, arrow = arrow(length = unit(0.25, "picas"))) +
  #annotate("text", x = (PCvalues$PC1/1000), y = (PCvalues$PC2/1000), label = PCvalues$Sample, 
  #size = 3, colour ="black", fontface = "bold") +
  ggtitle("PCA of ASVs (Aitchison distance of CLR vectors)") +
  theme_bw(base_size=20) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="right", legend.key.size = unit(0.8,"line"), plot.title = element_text(size=12),
        legend.title = element_text(size = 10), legend.text = element_text(size=8),
        axis.title = element_text(size=10), axis.text = element_text(size=10))

ggsave("PCA_ASV_clr_Treatment.svg", units = "cm", height=10, width=15, device="svg", dpi=300) 


#DNA DATA
pcx.DNA$rotation %>% 
  as.data.frame %>%
  rownames_to_column("ASV") %>%
  ggplot(aes(x=PC1,y=PC2)) + 
  geom_point(size=1, alpha=0.3, show.legend = FALSE) +
  #geom_point(data=PC_h.e.LF, aes(x=PC1,y=PC2), size=3, alpha=0.3, color= 'blue1', show.legend = FALSE) +
  #geom_point(data=PC_h.e.CC, aes(x=PC1,y=PC2), size=3, alpha=0.3, color= 'brown', show.legend = FALSE) +
  geom_segment(data = PCvalues.DNA, aes(x = 0, y = 0, xend = (PC1/500), yend = (PC2/500), color = Treatment), 
               alpha=0.7, arrow = arrow(length = unit(0.25, "picas"))) +
  #annotate("text", x = (PCvalues$PC1/1000), y = (PCvalues$PC2/1000), label = PCvalues$Sample, 
  #size = 3, colour ="black", fontface = "bold") +
  ggtitle("PCA of ASVs (Aitchison distance of CLR vectors)") +
  theme_bw(base_size=20) + 
  labs(x=paste0("PC1: ",round(var_explained.DNA[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained.DNA[2]*100,1),"%")) +
  theme(legend.position="right", legend.key.size = unit(0.8,"line"), plot.title = element_text(size=12),
        legend.title = element_text(size = 10), legend.text = element_text(size=8),
        axis.title = element_text(size=10), axis.text = element_text(size=10))

ggsave("PCA_ASV_clr_Treatment_DNA.svg", units = "cm", height=10, width=15, device="svg", dpi=300) 


# BLOCK4
# Similarly, we can plot the individual samples in Figure 4 (`R_Block_4`) 
# Here the distances between points is proportional to the Euclidian distance of the CLR vectors of the samples. 
# This is referred to as the Aitchison distance.

# Plot of the individual samples from the SVD output.
# plot(pcx$x[,1], pcx$x[,2],
#      main="Samples", xlab=xlab, ylab=ylab,
#      col=factor(metadata$Century),
#     cex=2, pch=19)
pcx_samples <- pcx$x %>% 
  as.data.frame %>%
  rownames_to_column("PlotTime") %>%
  separate(PlotTime,c("Plot","Time"), "T") # avoid extra spaces in the cols with " _ " not "_"

pcx_samples <- pcx_samples %>%
  mutate(Type= case_when(
    grepl("^D.*", pcx_samples$Plot) ~ "DNA",
    grepl("^cD.*", pcx_samples$Plot) ~ "RNA",
    grepl("^R.*", pcx_samples$Plot) ~ "DNA",
    grepl("^cR.*", pcx_samples$Plot) ~ "RNA"))

pcx_samples <- pcx_samples %>% # be aware that Rewetted T1 samples are Dry samples
  mutate(Treatment= case_when(
    Plot == "cD1" ~ "Dry",
    Plot == "cD3" ~ "Dry",
    Plot == "cD5" ~ "Dry",
    Plot == "cR6" ~ "Rewetted",
    Plot == "cR7" & Time == "1" ~ "Dry",
    Plot == "cR7" & Time == "2" ~ "Rewetted",
    Plot == "cR8" & Time == "1" ~ "Dry",
    Plot == "cR8" & Time == "3" ~ "Rewetted",
    Plot == "cR9" ~ "Rewetted",
    Plot == "D1" ~ "Dry",
    Plot == "D2" ~ "Dry",
    Plot == "D3" ~ "Dry",
    Plot == "D4" ~ "Dry",
    Plot == "D5" ~ "Dry",
    Plot == "R6" & Time == "1" ~ "Dry",
    Plot == "R6" & Time == "3" ~ "Rewetted",
    Plot == "R6" & Time == "7" ~ "Rewetted",
    Plot == "R6" & Time == "F" ~ "Rewetted",
    Plot == "R7" & Time == "1" ~ "Dry",
    Plot == "R7" & Time == "3" ~ "Rewetted",
    Plot == "R7" & Time == "7" ~ "Rewetted",
    Plot == "R7" & Time == "F" ~ "Rewetted",
    Plot == "R8" ~ "Rewetted",
    Plot == "R9" & Time == "1" ~ "Dry",
    Plot == "R9" & Time == "3" ~ "Rewetted",
    Plot == "R9" & Time == "7" ~ "Rewetted",
    Plot == "R9" & Time == "F" ~ "Rewetted",
    Plot == "R10" ~ "Rewetted"
  ))

pcx_samples <- pcx_samples %>%
  mutate(TreatTime= case_when(
    Time == "1" & Treatment =="Dry" ~ "DryT1",
    Time == "3" & Treatment =="Dry" ~ "DryT3",
    Time == "7" & Treatment =="Dry" ~ "DryT7",
    Time == "F" & Treatment =="Dry" ~ "DryT11",
    Time == "2" & Treatment =="Rewetted" ~ "RewettedT3",
    Time == "3" & Treatment =="Rewetted" ~ "RewettedT3",
    Time == "7" & Treatment =="Rewetted" ~ "RewettedT7",
    Time == "F" & Treatment =="Rewetted" ~ "RewettedT11"
  ))


# transform to a ordered factor vector. With gsub("\\s*", "") you remove extra spaces from the column avoiding the appearance of NAs due to a mismatch btw the levels of the vector and the ones you input.
pcx_samples$Treatment <- factor(pcx_samples$Treatment, levels = c("Dry","Rewetted"))

pPCA <-  ggplot(pcx_samples, aes(x=PC1,y=PC2)) + geom_point(aes(color=Type), size=10, alpha=0.6) +
  ggtitle("PCA of samples (Aitchison distance of CLR vectors)") +
  geom_text(aes(label = Treatment, color = Type), hjust = "bottom", nudge_x = 0, size=4, vjust = "bottom", nudge_y = 1 ) +
  theme_bw(base_size=20) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="bottom", legend.key.size = unit(0.3,"line"),
        legend.title = element_text(size = 18, color="black", face = "bold"), legend.text = element_text(size=16),
        axis.title = element_text(size=18, color="black", face = "bold"), axis.text = element_text(size=16,color="black"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        panel.grid =element_blank()) +
  scale_color_manual(name="Treatment",values=c("gold","darkorchid1"))


ggsave("PCA_SAMPLE_clr_Type.svg", height=10, width=15, device="svg", dpi=300) 

# multiple plots in one page 
library(ggpubr)
ggarrange(pS , pH, pPCA, #+ rremove("x.text")
          labels = c("A", "B", "C"), label.x = c(0, 0, 0), label.y = c(1, 1, 1),
          align = "hv", # align plot vertically and horiz.
          font.label = list(size = 30, color = "black", face = "bold"),
          common.legend = FALSE, legend = NULL,
          ncol = 1, nrow = 3)
ggsave("AlphaBetaDiversity_Rewet.svg", height=18, width=10, device="svg", dpi = 300) 


##DNA DATA##

pcx_samples.DNA <- pcx.DNA$x %>% 
  as.data.frame %>%
  rownames_to_column("PlotTime") %>%
  separate(PlotTime,c("Plot","Time"), "T") # avoid extra spaces in the cols with " _ " not "_"

pcx_samples.DNA <- pcx_samples.DNA %>%
                    mutate(Treatment= case_when(
                      Plot == "D1" ~ "Dry",
                      Plot == "D2" ~ "Dry",
                      Plot == "D3" ~ "Dry",
                      Plot == "D4" ~ "Dry",
                      Plot == "D5" ~ "Dry",
                      Plot == "R6" & Time == "1" ~ "Dry",
                      Plot == "R6" & Time == "3" ~ "Rewetted",
                      Plot == "R6" & Time == "7" ~ "Rewetted",
                      Plot == "R6" & Time == "F" ~ "Rewetted",
                      Plot == "R7" & Time == "1" ~ "Dry",
                      Plot == "R7" & Time == "3" ~ "Rewetted",
                      Plot == "R7" & Time == "7" ~ "Rewetted",
                      Plot == "R7" & Time == "F" ~ "Rewetted",
                      Plot == "R8" ~ "Rewetted",
                      Plot == "R9" & Time == "1" ~ "Dry",
                      Plot == "R9" & Time == "3" ~ "Rewetted",
                      Plot == "R9" & Time == "7" ~ "Rewetted",
                      Plot == "R9" & Time == "F" ~ "Rewetted",
                      Plot == "R10" ~ "Rewetted"
                    ))

pcx_samples.DNA <- pcx_samples.DNA %>%
                    mutate(Time= case_when(
                      Time == "1" ~ "T1",
                      Time == "2" ~ "T3",
                      Time == "3" ~ "T3",
                      Time == "7" ~ "T7",
                      Time == "F" ~ "T11"
                    ))
pcx_samples.DNA <- pcx_samples.DNA %>%
                    mutate(TreatTime= case_when(
                      Time == "T1" & Treatment =="Dry" ~ "DryT1",
                      Time == "T3" & Treatment =="Dry" ~ "DryT3",
                      Time == "T7" & Treatment =="Dry" ~ "DryT7",
                      Time == "T11" & Treatment =="Dry" ~ "DryT11",
                      Time == "T3" & Treatment =="Rewetted" ~ "RewettedT3",
                      Time == "T7" & Treatment =="Rewetted" ~ "RewettedT7",
                      Time == "T11" & Treatment =="Rewetted" ~ "RewettedT11"
                    ))

# transform to a ordered factor vector. With gsub("\\s*", "") you remove extra spaces from the column avoiding the appearance of NAs due to a mismatch btw the levels of the vector and the ones you input.
pcx_samples.DNA$Treatment <- factor(pcx_samples.DNA$Treatment, levels = c("Dry","Rewetted"))
pcx_samples.DNA$Time <- factor(pcx_samples.DNA$Time, levels = c("T1","T3","T7","T11"))

# add rownames
rownames(pcx_samples.DNA) <- gsub(" ","",paste(pcx_samples.DNA$Plot,pcx_samples.DNA$Time))

# check the order of samples
order_idx7 <- match(rownames(pcx_samples.DNA), rownames(envdata.DNA.f))

# Reorder filtered data frame to have the sample names in the same order as the metadata data frame
envdata.DNA.f  <- envdata.DNA.f[order_idx7,]

# VERIFY SAMPLE ORDER!!!
all(rownames(pcx_samples.DNA) == rownames(envdata.DNA.f)) #[1] TRUE


pPCA.DNA <-  ggplot(pcx_samples.DNA, aes(x=PC1,y=PC2)) + geom_point(aes(color=Time), size=10, alpha=0.6) +
            ggtitle("PCA of samples (Aitchison distance of CLR vectors)") +
            geom_text(aes(label = Treatment, color = Time), hjust = "bottom", nudge_x = 0, size=4, vjust = "bottom", nudge_y = 1 ) +
            theme_bw(base_size=20) + 
            labs(x=paste0("PC1: ",round(var_explained.DNA[1]*100,1),"%"),
                 y=paste0("PC2: ",round(var_explained.DNA[2]*100,1),"%")) +
            theme(legend.position="bottom", legend.key.size = unit(0.3,"line"),
                  legend.title = element_text(size = 18, color="black", face = "bold"), legend.text = element_text(size=16),
                  axis.title = element_text(size=18, color="black", face = "bold"), axis.text = element_text(size=16,color="black"),
                  panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                  panel.grid =element_blank()) +
            scale_color_manual(name="Time",values=c("orange","darkolivegreen3","blue","red4")) 

ggsave("PCA_SAMPLE_clr_Time_DNA.svg", height=10, width=15, device="svg", dpi=300)


# environmental parameters (according to section 3.2 draft)
set.seed(12)
envfit_result <- envfit(pcx_samples.DNA[,3:4], envdata.DNA.f[,c(7, 9:12,15,17)], perm = 999)       
# ***VECTORS
# 
#                           PC1      PC2     r2 Pr(>r)    
# emissions_CO2_mmolm2day  0.47990  0.87732 0.4019  0.004 **
# Soil._temp              -0.96924 -0.24612 0.2842  0.020 * 
# pH                       0.37134 -0.92850 0.1564  0.104   
# Water_activity          -0.56673 -0.82391 0.3277  0.010 **
# Conductivity             0.59194  0.80598 0.3101  0.009 **
# Water_gW.gDW             0.91151  0.41127 0.0146  0.832   
# AFDW_.DW                 0.97904  0.20368 0.1278  0.170       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 999

# To get coordinates of the arrows those direction values should be multiplied
#by square root of r2 values that are stored in vectors$r
envfit_data <- as.data.frame(envfit_result$vectors$arrows*sqrt(envfit_result$vectors$r))
envfit_data$Environ<-rownames(envfit_data)


pPCA.DNA +  geom_text(data=envfit_data,aes(x=PC1*50, y=PC2*50, label = rownames(envfit_data)),nudge_x = -0.5, nudge_y = 1) +
            geom_segment(data=envfit_data,aes(x = 0, y = 0, xend = PC1*50, yend = PC2*50), 
                          arrow = arrow(length = unit(0.3, "cm")), color = "red", linewidth = 1.2)
ggsave("PCA_SAMPLE_clr_Time_DNA_EnvironVar.svg", height=10, width=15, device="svg", dpi=300)

envfit_plot <- ggplot(envfit_data, aes(x = PC1, y = PC2)) +
  geom_text(aes(label = rownames(envfit_data), x = PC1*2, y = PC2*2)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), 
               arrow = arrow(length = unit(0.3, "cm")), color = "red") +
              xlim(c(-1, 1)) + ylim(c(-1, 1)) +
              theme_minimal()
   



# BLOCK5
# Similarly, we can calculate an Aitchison distance matrix as in `R_Block_5` to generate the dendrogram 
# with sample names colored depending on if they fall within the specific century. 
# Dendrogram plot of the individual samples from distance matrix output.


# Calculate Aitchison distance
dist.clr <- dist(t(f.clr), method = "euclidean")
dist.clr.DNA <- dist(t(f.clr.DNA), method = "euclidean") #DNA data


#################
### PERMANOVA ###
#################
library(RVAideMemoire)
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# PERMANOVA with adonis and adonis2 functions from vegan package
# one-way PERMANOVA
adonis2(formula=dist.clr ~ Type, data=metadata.f, permutations = 999, sqrt.dist = FALSE, add = FALSE, by = "terms")
# Type significant

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist.clr ~ Type, data = metadata.f, permutations = 999, sqrt.dist = FALSE, add = FALSE, by = "terms")
# Df SumOfSqs      R2      F Pr(>F)    
# Type      1     8960 0.06076 2.5228  0.001 ***
#   Residual 39   138515 0.93924                  
# Total    40   147475 1.00000 

# three-way PERMANOVA
adonis2(formula=dist.clr ~ Type*Treatment*Time, data=metadata.f, permutations = 999, sqrt.dist = FALSE, add = FALSE, by = "terms")
# Interaction Type and Treatment significant

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist.clr ~ Type * Treatment * Time, data = metadata.f, permutations = 999, sqrt.dist = FALSE, add = FALSE, by = "terms")

#                     Df SumOfSqs      R2      F Pr(>F)    
# Type                 1     8960 0.06076 2.6923  0.001 ***
# Treatment            1     8790 0.05960 2.6412  0.001 ***
# Time                 3    11403 0.07732 1.1421  0.172    
# Type:Treatment       1     6724 0.04559 2.0203  0.008 ** 
# Type:Time            3     9967 0.06759 0.9983  0.472    
# Treatment:Time       2     5580 0.03784 0.8384  0.795    
# Type:Treatment:Time  1     2864 0.01942 0.8606  0.681    
# Residual            28    93186 0.63188                  
# Total               40   147475 1.00000  

# pair-wise comparisons 
pairwise.perm.manova(dist.clr, metadata.f$Type, nperm=999, p.method='fdr', R2 = TRUE)
pairwise.adonis2(dist.clr ~ TypeTreat, data=metadata.f, nperm=999)
# $parent_call
# [1] "dist.clr ~ TypeTreat , strata = Null , permutations 999"
# 
# $`DNA _ Dry_vs_cDNA _ Dry`
# Df SumOfSqs      R2      F Pr(>F)    
# TypeTreat  1     9700 0.12224 2.7853  0.001 ***
#   Residual  20    69655 0.87776                  
# Total     21    79355 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $`DNA _ Dry_vs_DNA _ Rewetted`
# Df SumOfSqs      R2      F Pr(>F)    
# TypeTreat  1    11536 0.10989 3.3335  0.001 ***
#   Residual  27    93439 0.89011                  
# Total     28   104975 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $`DNA _ Dry_vs_cDNA _ Rewetted`
# Df SumOfSqs      R2      F Pr(>F)   
# TypeTreat  1    10278 0.15095 3.2001  0.002 **
#   Residual  18    57814 0.84905                 
# Total     19    68092 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $`cDNA _ Dry_vs_DNA _ Rewetted`
# Df SumOfSqs      R2      F Pr(>F)
# TypeTreat  1     3673 0.05331 1.0698  0.296
# Residual  19    65237 0.94669              
# Total     20    68910 1.00000              
# 
# $`cDNA _ Dry_vs_cDNA _ Rewetted`
# Df SumOfSqs      R2      F Pr(>F)
# TypeTreat  1     3928 0.11713 1.3266  0.141
# Residual  10    29611 0.88287              
# Total     11    33539 1.00000              
# 
# $`DNA _ Rewetted_vs_cDNA _ Rewetted`
# Df SumOfSqs     R2      F Pr(>F)   
# TypeTreat  1     6505 0.1086 2.0712  0.009 **
#   Residual  17    53395 0.8914                 
# Total     18    59901 1.0000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1          

## Rewetted diff cDNA vs DNA
## Dry diff cDNA vs DNA
## DNA dry vs rewetted
## DNA dry vs cDNA rewetted


## DNA DATA

# two-way PERMANOVA 
adonis2(formula=dist.clr.DNA ~ Treatment*Time, data=metadata.DNA.f, permutations = 999, sqrt.dist = FALSE, add = FALSE, by = "terms")
# Time significant

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist.clr.DNA ~ Treatment * Time, data = metadata.DNA.f, permutations = 999, sqrt.dist = FALSE, add = FALSE, by = "terms")
# Df SumOfSqs      R2      F Pr(>F)  
# Treatment       1     3529 0.03690 1.0901  0.270  
# Time            3    12309 0.12870 1.2674  0.038 *
#   Treatment:Time  2     8576 0.08968 1.3246  0.053 .
# Residual       22    71220 0.74472                
# Total          28    95634 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1          

# pair-wise comparisons NO SIGNIFICANT PAIR

pw <- pairwise.adonis2(dist.clr.DNA ~ Time, data=metadata.DNA.f, nperm=999)

# $parent_call
# [1] "dist.clr.DNA ~ Time , strata = Null , permutations 999"
# 
# $T1_vs_T3
# Df SumOfSqs      R2      F Pr(>F)
# Time      1     3463 0.09796 0.9774  0.446
# Residual  9    31890 0.90204              
# Total    10    35353 1.00000              
# 
# $T1_vs_T7
# Df SumOfSqs      R2     F Pr(>F)
# Time      1     3325 0.06348 0.949  0.539
# Residual 14    49049 0.93652             
# Total    15    52374 1.00000             
# 
# $T1_vs_T11
# Df SumOfSqs      R2      F Pr(>F)
# Time      1     2283 0.05669 0.7212   0.88
# Residual 12    37989 0.94331              
# Total    13    40272 1.00000              
# 
# $T3_vs_T7
# Df SumOfSqs      R2     F Pr(>F)
# Time      1     3164 0.06246 0.866  0.754
# Residual 13    47500 0.93754             
# Total    14    50665 1.00000             
# 
# $T3_vs_T11
# Df SumOfSqs      R2      F Pr(>F)
# Time      1     3974 0.09832 1.1995  0.214
# Residual 11    36440 0.90168              
# Total    12    40414 1.00000              
# 
# $T7_vs_T11
# Df SumOfSqs     R2      F Pr(>F)
# Time      1     3991 0.0693 1.1914   0.18
# Residual 16    53600 0.9307              
# Total    17    57590 1.0000         


###############################
##### REDUNDANCY ANALYSIS #####
###############################

# if(!"devtools" %in% installed.packages()){
#   install.packages("devtools")
# }
# devtools::install_github("gmteunisse/fantaxtic")
#devtools::install_github("gmteunisse/ggnested")

library(fantaxtic)

library(vegan)

# EXPLANATORY VARIABLES

envDNA <- envdata.DNA.f[,c(7, 9:12,15,17)]
envRNA <- envdata.RNA.f[,c(7, 9:12,15,17)]

# check for correlations between variables, as the constrained ordination methods we will be using are 
# highly sensitive to collinearities in the explanatory matrix. 

# We can visually look for correlations between variables:
heatmap(abs(cor(envDNA)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topleft", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))

# the environmental data are all in different units and will therefore need to be standardized prior to performing any ordinations. 
# We can once again use the decostand() function to standardize the environmental variables.

# Scale and center variables
envDNA.z <- decostand(envDNA, method = "standardize")
envRNA.z <- decostand(envRNA, method = "standardize")

# Variables are now centered around a mean of 0
round(apply(envDNA.z, 2, mean), 1)
round(apply(envRNA.z, 2, mean), 1)

# and scaled to have a standard deviation of 1
apply(envDNA.z, 2, sd)
apply(envRNA.z, 2, sd)



# RESPONSE VARIABLE
# phylo object with DNA data 
psq.1.DNA
otu_table(psq.1.DNA)

# phylo object with RNA data 
psq.1.RNA


# Top 20 most abundant Orders over all samples by median abundance
# fantaxtic::top_taxa(psq.1.DNA, tax_level = "Order", n_taxa = 20, FUN = median, include_na_taxa = TRUE)
# This function identifies the top n taxa by some metric (e.g. mean, median, variance, etc.) in a phyloseq object. 
# It outputs a table with the top taxa, as well as a phyloseq object in which all other taxa have been merged into a single taxon.
# Top 5 most abundant Orders PER SAMPLE using absolute abundances
ps5ordDNA <- fantaxtic::top_taxa(psq.1.DNA, tax_level = "Order", n_taxa = 5, by_proportion = FALSE, include_na_taxa = TRUE, grouping = "sample_id")
ps5ordRNA <- fantaxtic::top_taxa(psq.1.RNA, tax_level = "Order", n_taxa = 10, by_proportion = FALSE, include_na_taxa = TRUE, grouping = "sample_id")

# Save data frames 
ord5DNA <- as.data.frame(otu_table(ps5ordDNA$ps_obj))
taxord5 <- as.data.frame(tax_table(ps5ordDNA$ps_obj))
ord5DNAtop <- ps5ordDNA$top_taxa

ord5RNA <- as.data.frame(otu_table(ps5ordRNA$ps_obj))
taxord5RNA <- as.data.frame(tax_table(ps5ordRNA$ps_obj))
ord5RNAtop <- ps5ordRNA$top_taxa

# Create a variable PhyOrd
taxord5$PhyOrd <- paste(taxord5$Phylum,";",taxord5$Order)
taxord5RNA$PhyOrd <- paste(taxord5RNA$Phylum,";",taxord5RNA$Order)

# Check rownames order
all(rownames(ord5DNA) == rownames(taxord5)) #TRUE
all(rownames(ord5RNA) == rownames(taxord5RNA)) #TRUE

# Change rownames of otu table
rownames(ord5DNA) <- taxord5$PhyOrd
rownames(ord5RNA) <- taxord5RNA$PhyOrd

# Transpose to get samples x Orders
dfOrd5DNA <- as.data.frame(t(ord5DNA))
dfOrd5RNA <- as.data.frame(t(ord5RNA))

# change row names to match env data
dfOrd5DNA <- dfOrd5DNA %>% 
            rownames_to_column() %>%
            separate(rowname,c("Plot","Time"), "T") # avoid extra spaces in the cols with " _ " not "_"

dfOrd5DNA <- dfOrd5DNA %>%
            mutate(Time= case_when(
              Time == "1" ~ "T1",
              Time == "2" ~ "T3",
              Time == "3" ~ "T3",
              Time == "7" ~ "T7",
              Time == "F" ~ "T11"
            ))

dfOrd5DNA$PlotTime <- paste0(dfOrd5DNA$Plot,dfOrd5DNA$Time)
          
rownames(dfOrd5DNA) <- dfOrd5DNA$PlotTime

dfOrd5DNA$PlotTime <- NULL
dfOrd5DNA$Plot <- NULL
dfOrd5DNA$Time <- NULL

dfOrd5RNA <- dfOrd5RNA %>% 
  rownames_to_column() %>%
  separate(rowname,c("Plot","Time"), "T") # avoid extra spaces in the cols with " _ " not "_"

dfOrd5RNA <- dfOrd5RNA %>%
  mutate(Time= case_when(
    Time == "1" ~ "T1",
    Time == "2" ~ "T3",
    Time == "3" ~ "T3",
    Time == "7" ~ "T7",
    Time == "F" ~ "T11"
  ))

dfOrd5RNA <- dfOrd5RNA %>%
  mutate(Plot= case_when(
    Plot == "cD3" ~ "D3",
    Plot == "cD1" ~ "D1",
    Plot == "cD5" ~ "D5",
    Plot == "cR6" ~ "R6",
    Plot == "cR7" ~ "R7",
    Plot == "cR8" ~ "R8",
    Plot == "cR9" ~ "R9"
  ))

dfOrd5RNA$PlotTime <- paste0(dfOrd5RNA$Plot,dfOrd5RNA$Time)

rownames(dfOrd5RNA) <- dfOrd5RNA$PlotTime

dfOrd5RNA$PlotTime <- NULL
dfOrd5RNA$Plot <- NULL
dfOrd5RNA$Time <- NULL


### FUNCTION ###
# # Change colnames
# # Define a function to extract the desired parts from the original column names
# clean_column_names <- function(colname) {
#   parts <- unlist(strsplit(colname, "\\.\\."))  # Split the column name by ".."
#   cleaned <- paste(parts[c(2, 4)], collapse = "_")  # Extract the 2nd and 4th parts and join with "_"
#   return(cleaned)
# }
# 
# # Apply the function to get the new column names
# new_colnames <- lapply(colnames(ordDNA), clean_column_names)
# 
# # Assign the new column names to the dataframe
# colnames(ordDNA) <- unlist(new_colnames)
################


# Order rows as env data
order_idx60 <- match(rownames(envDNA),rownames(dfOrd5DNA))
dfOrd5DNA <- dfOrd5DNA[order_idx60,]
all(rownames(envDNA) == rownames(dfOrd5DNA)) # TRUE

order_idx61 <- match(rownames(envRNA),rownames(dfOrd5RNA))
dfOrd5RNA <- dfOrd5RNA[order_idx61,]
all(rownames(envRNA) == rownames(dfOrd5RNA)) # TRUE

# Transform abundance data 
# replace 0 values with an estimate: NO 0 values in the dataframe

# Generate the centered log-ratio transformed data
# samples by row (MARGIN=1)
# Option in decostand(method = "rclr") that deal with zeroes without adding a pseudocount
dfOrd5DNA.clr <- decostand(dfOrd5DNA, method = "clr", MARGIN = 1)
dim(dfOrd5DNA.clr) # [1] 29 17

dfOrd5RNA.clr <- decostand(dfOrd5RNA, method = "rclr", MARGIN = 1)
dim(dfOrd5RNA.clr) # [1] 12 23

# Check samples order
all(rownames(dfOrd5DNA.clr) == rownames(envDNA.z)) #[1] TRUE
all(rownames(dfOrd5RNA.clr) == rownames(envRNA.z)) #[1] TRUE


# Run the Redundancy analysis DNA
rda_model <- rda(dfOrd5DNA.clr ~ ., 
                 data = envDNA.z)

summary(rda_model)

# Partitioning of variance:
#   Inertia Proportion
# Total           6.797     1.0000
# Constrained     2.303     0.3389
# Unconstrained   4.494     0.6611

# variance of Y explained by X (33.9%)
# unexplained variance (66.1%)
# The included environmental variables explain 33.9% of the variation in bacterial community composition across sites.

# Eigenvalues, and their contribution to the variance 
# 
# Importance of components:
#                       RDA1    RDA2    
# Eigenvalue            1.5078 0.44305 
# Proportion Explained  0.2218 0.06518 
# Cumulative Proportion 0.2218 0.28700 


# Forward selection of variables:
fwd.sel <- ordiR2step(rda(dfOrd5DNA.clr ~ 1, data = envDNA.z), # lower model limit (simple!)
                        scope = formula(rda_model), # upper model limit (the "full" model)
                        direction = "forward",
                        R2scope = TRUE, # can't surpass the "full" model's R2
                        pstep = 1000,
                        trace = FALSE) # change to TRUE to see the selection process!

# Check new model with foeward-selected variables
fwd.sel$call

# rda(formula = dfOrd5DNA.clr ~ AFDW_.DW, data = envDNA.z)

# Write new model
rda_model.signif <- rda(formula = dfOrd5DNA.clr ~ AFDW_.DW, data = envDNA.z)

# check the adjusted R2 LOWER THAN THE MODEL WITH ALL ENV VARIABLES
RsquareAdj(rda_model.signif)
# $r.squared
# [1] 0.100336
# 
# $adj.r.squared
# [1] 0.06701506

# Significance Testing
# Full Model
anova.cca(rda_model, step =1000)
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = dfOrd5DNA.clr ~ emissions_CO2_mmolm2day + Soil._temp + pH + Water_activity + Conductivity + Water_gW.gDW + AFDW_.DW, data = envDNA.z)
#         Df Variance      F Pr(>F)  
# Model     7   2.3033 1.5376  0.066 .
# Residual 21   4.4940   

# Each variable
anova.cca(rda_model, step =1000, by = "term")
# Permutation test for rda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = dfOrd5DNA.clr ~ emissions_CO2_mmolm2day + Soil._temp + pH + Water_activity + Conductivity + Water_gW.gDW + AFDW_.DW, data = envDNA.z)
#                         Df Variance      F Pr(>F)  
# emissions_CO2_mmolm2day  1   0.5296 2.4746  0.049 *
# Soil._temp               1   0.4586 2.1430  0.095 .
# pH                       1   0.1949 0.9105  0.475  
# Water_activity           1   0.1973 0.9220  0.408  
# Conductivity             1   0.4897 2.2882  0.060 .
# Water_gW.gDW             1   0.0975 0.4558  0.806  
# AFDW_.DW                 1   0.3357 1.5689  0.164  
# Residual                21   4.4940 

# Each canonical axis
anova.cca(rda_model, step =1000, by = "axis")
# Permutation test for rda under reduced model
# Forward tests for axes
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = dfOrd5DNA.clr ~ emissions_CO2_mmolm2day + Soil._temp + pH + Water_activity + Conductivity + Water_gW.gDW + AFDW_.DW, data = envDNA.z)
# Df Variance      F Pr(>F)  
# RDA1      1   1.5078 7.0456  0.058 .
# RDA2      1   0.4431 2.0703  0.832  
# RDA3      1   0.1558 0.7280  1.000  
# RDA4      1   0.0896 0.4186  1.000  
# RDA5      1   0.0645 0.3012  1.000  
# RDA6      1   0.0297 0.1386  1.000  
# RDA7      1   0.0130 0.0607  1.000  
# Residual 21   4.4940                

anova.cca(rda_model.signif, step =1000)
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = dfOrd5DNA.clr ~ AFDW_.DW, data = envDNA.z)
#         Df Variance      F Pr(>F)  
# Model     1   0.6820 3.0112  0.025 *
#   Residual 27   6.1153 

# RDA plot
# Type 1 scaling. Distances among objects reflect their similarities
ordiplot(rda_model, scaling = 1, type = "text")
# Sites (numbers) that are closer together have more similar communities.
# Species that are closer together occupy more sites in common.

# Type 2 scaling. Angles between objects reflect their correlation
ordiplot(rda_model, scaling = 2, type = "text")
# Longer arrows mean this variable strongly drives the variation in the community matrix.
# Arrows pointing in opposite directions have a negative relationship.
# Arrows pointing in the same direction have a positive relationship.

# Custom plots
## extract % explained by the first 2 axes
perc <- round(100*(summary(rda_model)$cont$importance[2, 1:2]), 2)

## extract scores - these are coordinates in the RDA space
sc_si <- scores(rda_model, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(rda_model, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(rda_model, display="bp", choices=c(1, 2), scaling=1)

sc_si2 <- scores(rda_model, display="sites", choices=c(1,2), scaling=2)
sc_sp2 <- scores(rda_model, display="species", choices=c(1,2), scaling=2)
sc_bp2 <- scores(rda_model, display="bp", choices=c(1, 2), scaling=2)

orders_abrev <- c("Actinomarinales","Cyanobacteriales","Rubrobacterales","Proteo-Uncult","Phormidesmiles","Other","Nitriliruptorales","Rhodobacterales",
                  "Cytophagales","Gemmatimonadales","Chitinophagales","BD2-11","Sphingomonadales","Tistrellales","Longimicrobiales","Ardenticatenales","Thermoanaerobaculales")
env_abrev <- c("CO2","Soil Temp.","pH","Water act.","Cond.","Water gW/gDW","AFDW %DW")

## Custom triplot, step by step

# Set up a blank plot with scaling, axes, and labels
plot(rda_model,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = TRUE,
     # set axis limits
     xlim = c(-2.5,2.5), 
     ylim = c(-2,2),
     # label the plot (title, and axes)
     main = "RDA - scaling 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
    )
# # add points for site scores
# points(sc_si, 
#        pch = 22, # set shape (here, square with a fill colour)
#        col = "black", # outline colour
#        bg = "steelblue", # fill colour
#        cex = 1.5) # size
# # add text labels for species abbreviations
# text(sc_si + c(0.05, 0.09), # adjust text coordinates to avoid overlap with points 
#      labels = rownames(sc_si), 
#      col = "black", 
#      font = 1, # not bold
#      cex = 1)
# add points for species scores
points(sc_sp, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black",
       bg = "black", 
       cex = 1.5)
# add text labels for species abbreviations
text(sc_sp + c(0.05, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = orders_abrev, 
     col = "black", 
     font = 1, # not bold
     cex = 1)
# add arrows for effects of the explanatory variables
arrows(0,0, # start them from (0,0)
       sc_bp[,1]*10, sc_bp[,2]*10, # end them at the score value
       col = "red", 
       lwd = 2)
# add text labels for arrows
text(x = sc_bp[,1]*10 -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2]*10 - 0.03, 
     labels = env_abrev, 
     col = "red", 
     cex = 1, 
     font = 1)


# Set up a blank plot with scaling, axes, and labels
plot(rda_model,
     scaling = 2, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = TRUE,
     # set axis limits
     xlim = c(-0.9,0.9), 
     ylim = c(-0.9,0.9),
     # label the plot (title, and axes)
     main = "RDA - scaling 2",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)

# add points for species scores
points(sc_sp2, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black",
       bg = "black", 
       cex = 1.5)
# add text labels for species abbreviations
text(sc_sp2 + c(0.05, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = orders_abrev, 
     col = "black", 
     font = 1, # not bold
     cex = 1)
# add arrows for effects of the explanatory variables
arrows(0,0, # start them from (0,0)
       sc_bp2[,1]*1.2, sc_bp2[,2]*1.2, # end them at the score value
       col = "red", 
       lwd = 2)
# add text labels for arrows
text(x = sc_bp2[,1]*1.2 -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp2[,2]*1.2 - 0.03, 
     labels = env_abrev, 
     col = "red", 
     cex = 1, 
     font = 1)


# Run the Redundancy analysis RNA
rda_modelR <- rda(dfOrd5RNA.clr ~ ., 
                 data = envRNA.z)

summary(rda_modelR)

# Partitioning of variance:
#   Inertia Proportion
# Total           31.93      1.000
# Constrained     17.31      0.542
# Unconstrained   14.62      0.458

# variance of Y explained by X (54.2%)
# unexplained variance (45.8%)
# The included environmental variables explain 54.2% of the variation in bacterial community composition across sites.

# Eigenvalues, and their contribution to the variance 
# 
# Importance of components:
#   RDA1   RDA2    
# Eigenvalue            6.4455 4.3564 
# Proportion Explained  0.2018 0.1364
# Cumulative Proportion 0.2018 0.3383


# Forward selection of variables:
fwd.selR <- ordiR2step(rda(dfOrd5RNA.clr ~ 1, data = envRNA.z), # lower model limit (simple!)
                      scope = formula(rda_modelR), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE) # change to TRUE to see the selection process!

# Check new model with foeward-selected variables
fwd.selR$call

# rda(formula = dfOrd5RNA.clr ~ 1, data = envRNA.z)

# Write new model
rda_model.signifR <- rda(formula = dfOrd5RNA.clr ~ 1, data = envRNA.z)

# check the adjusted R2 LOWER THAN THE MODEL WITH ALL ENV VARIABLES
RsquareAdj(rda_model.signifR)


# Significance Testing
# Full Model
anova.cca(rda_modelR, step =1000)

# Non signif  

# Each variable
anova.cca(rda_modelR, step =1000, by = "term")
# non-signif

# Each canonical axis
anova.cca(rda_modelR, step =1000, by = "axis")
# non-signif              



# Custom plots
## extract % explained by the first 2 axes
perc <- round(100*(summary(rda_modelR)$cont$importance[2, 1:2]), 2)

## extract scores - these are coordinates in the RDA space
sc_siR <- scores(rda_modelR, display="sites", choices=c(1,2), scaling=1)
sc_spR <- scores(rda_modelR, display="species", choices=c(1,2), scaling=1)
sc_bpR <- scores(rda_modelR, display="bp", choices=c(1, 2), scaling=1)

sc_si2R <- scores(rda_modelR, display="sites", choices=c(1,2), scaling=2)
sc_sp2R <- scores(rda_modelR, display="species", choices=c(1,2), scaling=2)
sc_bp2R <- scores(rda_modelR, display="bp", choices=c(1, 2), scaling=2)

orders_abrevR <- c("Other","Rubrobacterales","Cyanobacteriales","Phormidesmiles","Frankiales","Rhodobacterales",
                  "Euzebyales","Cytophagales","Gemmatimonadales","Microtrichales","Myxococcales","Tistrellales",
                  "Sphingomonadales","Rhodothermales","Nitriliruptorales","Oxyphotobacteria","Azospirillales",
                  "Longimicrobiales","Staphylococcales","SBR1031","Pseudomonadales","Corynebacteriales","Rhodospirillales")
env_abrev <- c("CO2","Soil Temp.","pH","Water act.","Cond.","Water gW/gDW","AFDW %DW")

## Custom triplot, step by step

# Set up a blank plot with scaling, axes, and labels
plot(rda_modelR,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = TRUE,
     # set axis limits
     xlim = c(-2,2), 
     ylim = c(-2,2),
     # label the plot (title, and axes)
     main = "RDA - scaling 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# # add points for site scores
# points(sc_si, 
#        pch = 22, # set shape (here, square with a fill colour)
#        col = "black", # outline colour
#        bg = "steelblue", # fill colour
#        cex = 1.5) # size
# # add text labels for species abbreviations
# text(sc_si + c(0.05, 0.09), # adjust text coordinates to avoid overlap with points 
#      labels = rownames(sc_si), 
#      col = "black", 
#      font = 1, # not bold
#      cex = 1)
# add points for species scores
points(sc_spR, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black",
       bg = "black", 
       cex = 1.5)
# add text labels for species abbreviations
text(sc_spR + c(0.05, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = orders_abrevR, 
     col = "black", 
     font = 1, # not bold
     cex = 1)
# add arrows for effects of the explanatory variables
arrows(0,0, # start them from (0,0)
       sc_bpR[,1]*7, sc_bpR[,2]*7, # end them at the score value
       col = "red", 
       lwd = 2)
# add text labels for arrows
text(x = sc_bpR[,1]*7 -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bpR[,2]*7 - 0.03, 
     labels = env_abrev, 
     col = "red", 
     cex = 1, 
     font = 1)


# Set up a blank plot with scaling, axes, and labels
plot(rda_modelR,
     scaling = 2, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = TRUE,
     # set axis limits
     xlim = c(-0.9,0.9), 
     ylim = c(-0.9,0.9),
     # label the plot (title, and axes)
     main = "RDA - scaling 2",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)

# add points for species scores
points(sc_sp2R, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black",
       bg = "black", 
       cex = 1.5)
# add text labels for species abbreviations
text(sc_sp2R + c(0.05, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = orders_abrevR, 
     col = "black", 
     font = 1, # not bold
     cex = 1)
# add arrows for effects of the explanatory variables
arrows(0,0, # start them from (0,0)
       sc_bp2R[,1]*1.2, sc_bp2R[,2]*1.2, # end them at the score value
       col = "red", 
       lwd = 2)
# add text labels for arrows
text(x = sc_bp2R[,1]*1.2 -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp2R[,2]*1.2 - 0.03, 
     labels = env_abrev, 
     col = "red", 
     cex = 1, 
     font = 1)

