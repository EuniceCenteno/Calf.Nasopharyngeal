
library(ggplot2)
library(tidyr) #separate function
library(reshape2) #melt function
library(dplyr)
library(naniar) # for replace_with_na_all function
library(data.table)
library(phyloseq)
library(qiime2R)
library(ggpubr)
library(tidyverse)

setwd("~/Desktop/eunice/PhD/CalfNasal/Qiime/OtherQiimeFiles/CorrectQiime/")
##

##### Taxonomy barplot
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
# Read the experimental design, and species classified documents
ASVs <- read_qza("grouped-tableNPRarified.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #3860 ASVs
#str(ASV_table)
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
colSums(ASV_table)

#Importing the metadata file
metadata <- read.table("CalfMetadata01.14.txt", header=TRUE, sep="\t")
str(metadata)
rownames(metadata) <- metadata$ID2
metadata$calf = as.factor(metadata$calf)
metadata$d <- as.numeric(metadata$d)
metadata$treatment <- as.factor(metadata$treatment)
metadata$Sickness <- as.factor(metadata$Sickness)
IDCorrect <- metadata$ID
IDCorrect2 <- metadata$ID2

#merging the abundance of each OTU with the metadata and the taxonomy file
##Adding taxonomy
#Taxonomy of each OTU
tax <- read_qza("CorrectTaxNoRarified.qza")
tax <- as.data.frame(tax$data)
tax2 = separate(tax, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#This warning means that some cells are empty and that R is replacing empty cells with NA. Now there are other cells that are unclassified that  say, for example `s__` or `g__`. All these need to be removed and replaced with `NA`. 
#All this is OK except that in future use of the taxonomy table, these ASVs will be ignored because they are not classified. Why are ASVs not classified? Its because there is not a close enough match in the database. Just because there is not a good match in the database does not mean they don’t exist, so I wanted to make sure this data was not lost. So in my new code, from lines 300 – 334 I make it so that ASVs that are unclassified at any level are classified as the lowest taxonomic level for which there is a classification.

#All the strings that need to be removed and replaced with NA
na_strings <- c(" s__", " g__", " f__", " o__", " c__")

tax3 = replace_with_na_all(tax2, condition = ~.x %in% na_strings)
#This code is great because an ASV that is unclassified at a certain level are all listed as `NA`.
#Unfortunately this command changed ou Feature.ID names

#Next, all these `NA` classifications with the last level that was classified
tax3[] <- t(apply(tax3, 1, zoo::na.locf))

tax3 <- as.data.frame(tax3)
row.names(tax3) <- tax3[,1]
tax3 = tax3[,-c(1)]
tax.clean <- as.data.frame(tax3)
tax.clean$OTUs <- rownames(tax.clean)
#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ASV_s),]

tax.final$Domain <- sub("D_0__*", "", tax.final[,1])
tax.final$Phylum <- sub("D_1__*", "", tax.final[,2])
tax.final$Phylum <- sub("D_2__*", "", tax.final[,2])
tax.final$Class <- sub("D_1__*", "", tax.final[,3])
tax.final$Class <- sub("D_2__*", "", tax.final[,3])
tax.final$Class <- sub("D_3__*", "", tax.final[,3])
tax.final$Order <- sub("D_1__*", "", tax.final[,4])
tax.final$Order <- sub("D_2__*", "", tax.final[,4])
tax.final$Order <- sub("D_3__*", "", tax.final[,4])
tax.final$Order <- sub("D_4__*", "", tax.final[,4])
tax.final$Family <- sub("D_1__*", "", tax.final[,5])
tax.final$Family <- sub("D_2__*", "", tax.final[,5])
tax.final$Family <- sub("D_3__*", "", tax.final[,5])
tax.final$Family <- sub("D_4__*", "", tax.final[,5])
tax.final$Family <- sub("D_5__*", "", tax.final[,5])
tax.final$Genus <- sub("D_1__*", "", tax.final[,6])
tax.final$Genus <- sub("D_2__*", "", tax.final[,6])
tax.final$Genus <- sub("D_3__*", "", tax.final[,6])
tax.final$Genus <- sub("D_4__*", "", tax.final[,6])
tax.final$Genus <- sub("D_5__*", "", tax.final[,6])
tax.final$Genus <- sub("D_6__*", "", tax.final[,6])
tax.final$Species <- sub("D_1__*", "", tax.final[,7])
tax.final$Species <- sub("D_2__*", "", tax.final[,7])
tax.final$Species <- sub("D_3__*", "", tax.final[,7])
tax.final$Species <- sub("D_4__*", "", tax.final[,7])
tax.final$Species <- sub("D_5__*", "", tax.final[,7])
tax.final$Species <- sub("D_6__*", "", tax.final[,7])
tax.final$Species <- sub("D_7__*", "", tax.final[,7])

#write.table(tax.final,"taxonomyNasal.txt",sep=",", row.names = FALSE) 
TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
row.names(TaxASV) <- TaxASV[,11]
TaxASV = TaxASV[,-c(1)]
#write.table(TaxASV,"TaxASV.txt",sep="\t", row.names = T, col.names = T)

#Preparing ASVtable
NewASVtable <- t(ASV_table)
NewASVtable <- merge(metadata, NewASVtable, by.x = "ID", by.y = 0)
str(NewASVtable)
NewASVtable$ID2 <- as.factor(NewASVtable$ID2)
levels(NewASVtable$ID2) <- list("d0_1"="d0_1",    "d7_1"="d7_1",    "ds0_1"="ds0_1",   "d0_2"="d0_2",    "d7_2"="d7_2",    "d14_2"="d14_2",   "d0_3"="d0_3",    "d7_3"="d7_3",    "d14_3"="d14_3",  
                                "d0_4"="d0_4",    "d7_4"="d7_4",    "d14_4"="d14_4",   "d0_5"="d0_5",    "d7_5"="d7_5",    "d14_5"="d14_5",   "d0_6"="d0_6",    "d7_6"="d7_6",    "d14_6"="d14_6", 
                                "d0_7"="d0_7",    "d7_7"="d7_7",    "d14_7"="d14_7",   "d0_8"="d0_8",    "d7_8"="d7_8",    "d14_8"="d14_8",   "d0_9"="d0_9",    "d7_9"="d7_9",    "ds0_9"="ds0_9",  
                                "ds1_9"="ds1_9",   "ds10_9"="ds10_9",  "d7_10"="d7_10",   "d14_10"="d14_10",  "d0_11"="d0_11",   "d7_11"="d7_11",   "d14_11"= "d14_11",  "d0_12"="d0_12",   "d7_12"="d7_12", 
                                "d14_12"="d14_12",  "d0_13"="d0_13",   "d7_13"="d7_13",   "d14_13"="d14_13",  "d0_14"="d0_14",   "d14_14"="d14_14",  "d0_15"="d0_15",   "d7_15"="d7_15",   "d14_15"="d14_15", 
                                "ds1_15"="ds1_15",  "ds5_15"= "ds5_15",  "ds10_15"="ds10_15", "d7_16"="d7_16",   "d14_16"="d14_16",  "d7_17"="d7_17",   "d14_17"="d14_17",  "ds0_18"="ds0_18",  "ds1_18"="ds1_18", 
                                "ds5_18"="ds5_18",  "ds10_18"="ds10_18", "d0_19"="d0_19",   "ds0_19"= "ds0_19",  "ds1_19"="ds1_19",  "ds5_19"="ds5_19",  "ds10_19"="ds10_19", "d0_20"="d0_20",   "d7_20"="d7_20",  
                                "d14_20"= "d14_20",  "ds0_21"="ds0_21",  "ds1_21"="ds1_21",  "ds5_21"="ds5_21",  "ds10_21"="ds10_21", "d0_22"="d0_22" ,   "d7_22"="d7_22",   "d14_22"="d14_22",  "ds1_23"="ds1_23", 
                                "ds5_23"="ds5_23",  "ds10_23"="ds10_23", "d7_24"="d7_24",   "d14_24"="d14_24",  "d7_25"="d7_25",   "d14_25"="d14_25",  "d7_26"="d7_26",   "d14_26"="d14_26",  "d7_27"="d7_27",  
                                "d7_28"= "d7_28",   "ds0_28"="ds0_28",  "ds1_28"="ds1_28",  "ds5_28"="ds5_28",  "ds10_28"="ds10_28", "ds0_29"="ds0_29",  "ds1_29"="ds1_29",  "ds5_29"="ds5_29",  "ds10_29"="ds10_29",
                                "d0_30"="d0_30",   "ds0_30"="ds0_30",  "ds1_30"="ds1_30",  "ds5_30"="ds5_30",  "ds10_30"="ds10_30")
NewASVtable <-NewASVtable[,-c(1:7,9)]
NewASVtable <- arrange(NewASVtable, ID2)
rownames(NewASVtable) <- NewASVtable$ID2
ASV_table <- NewASVtable[,-c(1)]
ASV_table <- as.matrix(ASV_table) #95 and 3860 
rowSums(ASV_table) #1889 ASVs

### now including rep seqs
seq <- read_qza("rep-seqs-NoPseudoNoRari.qza")
seq <- as.data.frame(seq$data)
seq1 <- merge(seq, TaxASV, by.x =0, by.y = "OTUs")
#write.csv(seq1, "seq1.csv")

#importing Control ASVs
ASVsC<- read_qza("tableControls.qza")
ASV_sC <- as.data.frame(ASVsC$data)
ASV_tableC <- as.data.frame(ASVsC$data) #274 ASVs
#str(ASV_table)
ASV_tableC$ASVnos <- paste0("ASV", 1:nrow(ASV_tableC))
ASV_tableC$ASVstring <- rownames(ASV_tableC)
rownames(ASV_tableC) <- ASV_tableC$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkeyC <- ASV_tableC[, (ncol(ASV_tableC)-1):ncol(ASV_tableC)] #the key withe the names
ASV_tableC <- ASV_tableC[,-(ncol(ASV_tableC)-1):-ncol(ASV_tableC)]
colSums(ASV_tableC)

#merging the abundance of each OTU with the metadata and the taxonomy file
##Adding taxonomy
#Taxonomy of each OTU
taxC <- read_qza("SilvaTaxControl.qza")
taxC <- as.data.frame(taxC$data)
tax2C = separate(taxC, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#This warning means that some cells are empty and that R is replacing empty cells with NA. Now there are other cells that are unclassified that  say, for example `s__` or `g__`. All these need to be removed and replaced with `NA`. 
#All this is OK except that in future use of the taxonomy table, these ASVs will be ignored because they are not classified. Why are ASVs not classified? Its because there is not a close enough match in the database. Just because there is not a good match in the database does not mean they don’t exist, so I wanted to make sure this data was not lost. So in my new code, from lines 300 – 334 I make it so that ASVs that are unclassified at any level are classified as the lowest taxonomic level for which there is a classification.

#All the strings that need to be removed and replaced with NA
na_strings <- c(" s__", " g__", " f__", " o__", " c__")

tax3C = replace_with_na_all(tax2C, condition = ~.x %in% na_strings)
#This code is great because an ASV that is unclassified at a certain level are all listed as `NA`.
#Unfortunately this command changed ou Feature.ID names

#Next, all these `NA` classifications with the last level that was classified
tax3C[] <- t(apply(tax3C, 1, zoo::na.locf))

tax3C <- as.data.frame(tax3C)
row.names(tax3C) <- tax3C[,1]
tax3C = tax3C[,-c(1:2)]
tax.cleanC <- as.data.frame(tax3C)
tax.cleanC$OTUs <- rownames(tax.cleanC)
#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.finalC = tax.cleanC[row.names(tax.cleanC) %in% row.names(ASV_sC),]

tax.finalC$Phylum <- sub("D_1__*", "", tax.finalC[,1])
tax.finalC$Phylum <- sub("D_2__*", "", tax.finalC[,1])
tax.finalC$Class <- sub("D_1__*", "", tax.finalC[,2])
tax.finalC$Class <- sub("D_2__*", "", tax.finalC[,2])
tax.finalC$Class <- sub("D_3__*", "", tax.finalC[,2])
tax.finalC$Order <- sub("D_1__*", "", tax.finalC[,3])
tax.finalC$Order <- sub("D_2__*", "", tax.finalC[,3])
tax.finalC$Order <- sub("D_3__*", "", tax.finalC[,3])
tax.finalC$Order <- sub("D_4__*", "", tax.finalC[,3])
tax.finalC$Family <- sub("D_1__*", "", tax.finalC[,4])
tax.finalC$Family <- sub("D_2__*", "", tax.finalC[,4])
tax.finalC$Family <- sub("D_3__*", "", tax.finalC[,4])
tax.finalC$Family <- sub("D_4__*", "", tax.finalC[,4])
tax.finalC$Family <- sub("D_5__*", "", tax.finalC[,4])
tax.finalC$Genus <- sub("D_1__*", "", tax.finalC[,5])
tax.finalC$Genus <- sub("D_2__*", "", tax.finalC[,5])
tax.finalC$Genus <- sub("D_3__*", "", tax.finalC[,5])
tax.finalC$Genus <- sub("D_4__*", "", tax.finalC[,5])
tax.finalC$Genus <- sub("D_5__*", "", tax.finalC[,5])
tax.finalC$Genus <- sub("D_6__*", "", tax.finalC[,5])
tax.finalC$Species <- sub("D_1__*", "", tax.finalC[,6])
tax.finalC$Species <- sub("D_2__*", "", tax.finalC[,6])
tax.finalC$Species <- sub("D_3__*", "", tax.finalC[,6])
tax.finalC$Species <- sub("D_4__*", "", tax.finalC[,6])
tax.finalC$Species <- sub("D_5__*", "", tax.finalC[,6])
tax.finalC$Species <- sub("D_6__*", "", tax.finalC[,6])
tax.finalC$Species <- sub("D_7__*", "", tax.finalC[,6])

#write.table(tax.final,"taxonomyNasal.txt",sep=",", row.names = FALSE) 
TaxASVC <- merge(tax.finalC, ASVkeyC, by.x = 0, by.y = "ASVstring")
row.names(TaxASVC) <- TaxASVC[,10]
TaxASVC = TaxASVC[,-c(1)]
# write.table(TaxASV,"TaxASV.txt",sep="\t", row.names = T, col.names = T)

#now rep seq
### now including rep seqs
seqC <- read_qza("rep-seqsControl.qza")
seqC <- as.data.frame(seqC$data)
seq1C <- merge(seqC, TaxASVC, by.x =0, by.y = "OTUs")
#write.csv(seq1C, "seq1C.csv")

##try to see if there is any sequence from the samples that is shared with the control
# we have 3860 sequences and 274 sequences in the control
cont <- merge(seq1, seq1C, by.x = "Row.names", by.y = "Row.names")
#0 shared based on the OTU names
#now lets try merging by sequence
cont <- merge(seq1, seq1C, by.x = "x", by.y = "x")

#no sequences based on similar sequences

### Calculate the relative abundance of the "contaminated ASVs"
contASV <- read.csv("Contaminaiton100percent.csv")

#---------------MIKING THE MODEL CON 
#generating a dataframe with all the response (days) from the samples
# Make one column for our outcome/response variable 
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq # [ 3860 taxa and 95 samples ]
colnames(tax_table(physeq_deseq))
## Filter any non-baxteria, chloroplast and mitochondria
physeq_deseq %>%
  subset_taxa(Family != "Mitochondria" & 
                Genus != "Mitochondria" &
                Species != "Mitochondria" &
                Order != "Chloroplast" &
                Family != "Chloroplast" &
                Genus != "Chloroplast" &
                Species != "Chloroplast") -> physeq_deseq
physeq_deseq ##[ 3838 taxa and 95 samples ]
## Random forest, we want to make the model so it classify samples based on the age
ntaxa(physeq_deseq) #total of 3838

physeq_deseq
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:11)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
#write.table(NewTax,"NewTax.txt",sep=",", row.names = TRUE)
NewASVtable <- prunetable
NewASVtable <- NewASVtable[,-c(2:11)]
row.names(NewASVtable) <- NewASVtable[,1]
#NewASVtable = NewASVtable[,-c(1)]
#write.csv(NewASVtable, "NewASVtable.csv")
NewASVtable = t(NewASVtable[,-c(1)]) #3838 and 8 

#### Separate by healthy and sick
healthy <- subset(metadata, Sickness=="Never")
str(healthy)
healthy$d <- as.factor(healthy$d)
levels(healthy$d) <- list("0"="0", "7"="1", "14"="2")

H_OTU <- NewASVtable[rownames(NewASVtable) %in% rownames(healthy),] #51, 3838

## -----------------------------------PHYLUM LEVEL, Healthy period
### CALCULATION OF THE ABUNDANCE OF EACH OTU  
## check the prop.table, it gives me NA in the sample Dairy 92
otu.summary <- prop.table((H_OTU), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:3838)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)


#merging the abundance of each OTU with the metadata and the taxonomy file
str(metadata)
str(melt_otu)
meta_otu2 <- merge(healthy, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu2)
meta_otu_tax <- merge(meta_otu2, NewTax, by.x = "ASV", by.y = 0)
str(meta_otu_tax)
order_groups <- healthy$ID2
meta_otu_tax$Row.names <- factor(meta_otu_tax$Row.names, levels = order_groups)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Sickness <- factor(meta_otu_tax$Sickness)
levels(meta_otu_tax$Sickness) <- list("Never"="Healthy", "Sick"="Got")
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$ASV <- factor(meta_otu_tax$ASV)
str(meta_otu_tax)

## relative abundance of all the ASVs
ASV <- meta_otu_tax %>% 
  group_by(ASVnos) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/51)*100) ## the total number of samples (51)
attach(ASV )
sum(ASV$Ave_Abundance)
ASV  <- ASV [order(-Ave_Abundance),]
#write.csv(ASV,"RabundASVs.csv")

ASV1 <- meta_otu_tax%>% 
  group_by(Row.names, d, ASVnos) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, ASVnos) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(ASV1)
sum(ASV1$taxa.average) #3
#write.csv(ASV1,"RabunASVsDay.csv")

#subset Clostridium sensu stricto, Escherichia Shigella, Lactobacillus and Mycoplasma
str(melt_otu)
str(contASV)

my_colors <- c(
  '#a6cee3','#1f78b4','#b3df8a','#33a03c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab3d6','#6a3d9a','#ffff99','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

ASV2 <- meta_otu_tax%>% 
  group_by(Row.names, d, Genus, ASV) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Genus, ASV) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(ASV2)
sum(ASV2$taxa.average) #3
x <-ASV2 %>% count(Genus) 
#write.csv(ASV2,"ContASVsDayGenus.csv")
asv <- c("Clostridium sensu stricto 1", "Lactobacillus", "Escherichia-Shigella")
filteredASV2 <- ASV2[ASV2$Genus %in% asv, ] # we should have 3664 ASVs
sum(filteredASV2 $taxa.average)
#write_tsv(filteredASV2, "filteredASVSpecies.tsv")
levels(filteredASV2$d)
levels(filteredASV2$d) <- list("0"="0", "7"="7", "14"="14")


#Plot the graph 
a <- ggplot(filteredASV2, aes(x = d, y = taxa.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  labs(fill = "All ASVs (Genus)") +
  ylim(c(0,0.5)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Day')

##now look only at the ASV
unique(contASV$ASV)
contaminants <- merge(contASV,filteredASV2, by.x = "ASV", by.y = "ASV")
unique(contaminants$ASV)
sum(contaminants$taxa.average) #< 1
levels(contaminants$d)
levels(contaminants$d) <- list("0"="0", "7"="7", "14"="14")
#write.csv(contaminants, "contaminantsH.csv")

b <- ggplot(contaminants, aes(x = d, y = taxa.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.5)) +
  labs(fill = "Only 100% similarity (Genus)") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Day')

c <- ggplot(contaminants, aes(x = d, y = taxa.average, fill =SharedASVs)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  labs(fill = "Only 100% similarity (ASV)") +
  ylim(c(0,0.5)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Day')

ggarrange(a,b,c, labels = c("a", "b","c"),vjust = 0.9,
          ncol=3)

### now check the sick period

# we have to divide them by antibiotic 
ant <- read.table("CalfMetadataSickAnimals.txt", header=TRUE, sep="\t")
rownames(ant) <- ant$ID2
ant$antibiotic <- as.factor(ant$antibiotic)
ant$d <- as.factor(ant$d)
levels(ant$d) <- list("0"="0", "1"="1", "5"="2", "10"="3")
Draxxin <- subset(ant, antibiotic=="Draxxin")
rownames(Draxxin) <- Draxxin$ID2
Nuflur <- subset(ant, antibiotic=="Nuflur")
rownames(Nuflur) <- Nuflur$ID2


#New to generate the OTU for the healthy period
D_OTU <- NewASVtable[rownames(NewASVtable) %in% rownames(Draxxin),] #22, 3838
str(D_OTU)
N_OTU <- NewASVtable[rownames(NewASVtable) %in% rownames(Nuflur),] #11, 3838
str(N_OTU)

## check the prop.table, it gives me NA in the sample Dairy 92
otu.summaryD <- prop.table((D_OTU), 1) 
str(otu.summaryD)
otu_abundD <- colSums(otu.summaryD)
otu_abund2D <- as.data.frame(otu_abundD)
otu.summaryD <- rbind(otu_abundD, otu.summaryD)
str(otu.summaryD)
otu.summary_sortedD <- otu.summaryD[,order(otu.summaryD[1,], decreasing = TRUE)]
str(otu.summary_sortedD)
melt_otuD <- reshape2::melt(otu.summary_sortedD[, c(1:3838)]) ###TOTAL NUMBER OF OTUS
str(melt_otuD)
colnames(melt_otuD) <- c("Sample", "ASV", "Abundance")
str(melt_otuD)
levels(melt_otuD$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
str(Draxxin)
str(melt_otuD)
meta_otuD <- merge(Draxxin, melt_otuD, by.x = 0, by.y = "Sample")
str(meta_otuD)
meta_otu_taxD <- merge(meta_otuD, NewTax, by.x = "ASV", by.y = 0)
str(meta_otu_taxD)
order_groupsD <- Draxxin$ID2
meta_otu_taxD$Row.names <- factor(meta_otu_taxD$Row.names, levels = order_groupsD)
summary(meta_otu_taxD$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu_taxD$Family <- factor(meta_otu_taxD$Family)
meta_otu_taxD$Sickness <- factor(meta_otu_taxD$Sickness)
levels(meta_otu_taxD$Sickness) <- list("Never"="Healthy", "Sick"="Got")
meta_otu_taxD$Genus <- factor(meta_otu_taxD$Genus)
meta_otu_taxD$Phylum <- factor(meta_otu_taxD$Phylum)
meta_otu_taxD$ASV <- factor(meta_otu_taxD$ASV)
str(meta_otu_taxD)

## relative abundance of all the ASVs
str(Draxxin)
ASVD <- meta_otu_taxD %>% 
  group_by(ASVnos) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/22)*100) ## the total number of samples (51)
attach(ASVD)
sum(ASVD$Ave_Abundance)
ASVD  <- ASVD[order(-Ave_Abundance),]
#write.csv(ASVD,"RabundASVsDraxxin.csv")

ASV1D <- meta_otu_taxD%>% 
  group_by(Row.names, d, ASVnos) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, ASVnos) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(ASV1D)
sum(ASV1D$taxa.average) #4
#write.csv(ASV1D,"RabunASVsDayDRAXXIN.csv")

#subset Clostridium sensu stricto, Escherichia Shigella, Lactobacillus and Mycoplasma
str(melt_otu)
str(contASV)

ASV2D <- meta_otu_taxD%>% 
  group_by(Row.names, d, Genus, ASV) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Genus, ASV) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(ASV2D)
sum(ASV2D$taxa.average) #4
#write.csv(ASV2D,"ContASVsDayGenusDRAXXIN.csv")
filteredASV2D <- ASV2D[ASV2D$Genus %in% asv, ] # we should have 3664 ASVs
sum(filteredASV2D$taxa.average)
#write_tsv(filteredASV2D, "filteredASVSpeciesDRAXXIN.tsv")
levels(filteredASV2D$d)
levels(filteredASV2D$d) <- list("0"="0", "1"="1", "5"="5", "10"="10")


#Plot the graph 
A <- ggplot(filteredASV2D, aes(x = d, y = taxa.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.5)) +
  labs(fill = "All ASVs (Genus)") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Day')

##now look only at the ASV
unique(contASV$ASV)
contaminantsD <- merge(contASV,filteredASV2D, by.x = "ASV", by.y = "ASV")
unique(contaminantsD$ASV)
sum(contaminantsD$taxa.average) #< 1.232868
levels(contaminantsD$d)
levels(contaminantsD$d) <- list("0"="0", "1"="1", "5"="5", "10"="10")

B <- ggplot(contaminantsD, aes(x = d, y = taxa.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.5)) +
  labs(fill = "Only 100% similarity (Genus)") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Day')

C <- ggplot(contaminantsD, aes(x = d, y = taxa.average, fill =SharedASVs)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.5)) +
  labs(fill = "Only 100% similarity (ASV)") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Day')

ggarrange(A,B,C, labels = c("a", "b","c"),vjust = 0.9,
          ncol=3)

#NOW NUFLUR
## check the prop.table, it gives me NA in the sample Dairy 92
otu.summaryN <- prop.table((N_OTU), 1) 
str(otu.summaryN)
otu_abundN <- colSums(otu.summaryN)
otu_abund2D <- as.data.frame(otu_abundN)
otu.summaryN <- rbind(otu_abundN, otu.summaryN)
str(otu.summaryN)
otu.summary_sortedN <- otu.summaryN[,order(otu.summaryN[1,], decreasing = TRUE)]
str(otu.summary_sortedN)
melt_otuN <- reshape2::melt(otu.summary_sortedN[, c(1:3838)]) ###TOTAL NUMBER OF OTUS
str(melt_otuN)
colnames(melt_otuN) <- c("Sample", "ASV", "Abundance")
str(melt_otuN)
levels(melt_otuN$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
str(Nuflur)
str(melt_otuN)
meta_otuN <- merge(Nuflur, melt_otuN, by.x = 0, by.y = "Sample")
str(meta_otuN)
meta_otu_taxN <- merge(meta_otuN, NewTax, by.x = "ASV", by.y = 0)
str(meta_otu_taxN)
order_groupsN <- Nuflur$ID2
meta_otu_taxN$Row.names <- factor(meta_otu_taxN$Row.names, levels = order_groupsN)
summary(meta_otu_taxN$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu_taxN$Family <- factor(meta_otu_taxN$Family)
meta_otu_taxN$Sickness <- factor(meta_otu_taxN$Sickness)
levels(meta_otu_taxN$Sickness) <- list("Never"="Healthy", "Sick"="Got")
meta_otu_taxN$Genus <- factor(meta_otu_taxN$Genus)
meta_otu_taxN$Phylum <- factor(meta_otu_taxN$Phylum)
meta_otu_taxN$ASV <- factor(meta_otu_taxN$ASV)
str(meta_otu_taxN)


## relative abundance of all the ASVs
str(Nuflur)
ASVN <- meta_otu_taxN %>% 
  group_by(ASVnos) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/11)*100) ## the total number of samples (51)
attach(ASVN)
sum(ASVN$Ave_Abundance)
ASVN  <- ASVN[order(-Ave_Abundance),]
#write.csv(ASVN,"RabundASVsNUFLUR.csv")

ASV1N <- meta_otu_taxN%>% 
  group_by(Row.names, d, ASVnos) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, ASVnos) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(ASV1N)
sum(ASV1N$taxa.average) #4
#write.csv(ASV1N,"RabunASVsDayNUFLUR.csv")

#subset Clostridium sensu stricto, Escherichia Shigella, Lactobacillus and Mycoplasma
str(melt_otu)
str(contASV)

ASV2N <- meta_otu_taxN%>% 
  group_by(Row.names, d, Genus, ASV) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Genus, ASV) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(ASV2N)
sum(ASV2N$taxa.average) #4
#write.csv(ASV2N,"ContASVsDayGenusNUFLUR.csv")
filteredASV2N <- ASV2N[ASV2N$Genus %in% asv, ] # we should have 3664 ASVs
sum(filteredASV2N$taxa.average)
#write_tsv(filteredASV2N, "filteredASVSpeciesNUFLUR.tsv")
levels(filteredASV2N$d)
levels(filteredASV2N$d) <- list("0"="0", "1"="1", "5"="5", "10"="10")


#Plot the graph 
X <- ggplot(filteredASV2N, aes(x = d, y = taxa.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.5)) +
  labs(fill = "All ASVs (Genus)") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Day')

##now look only at the ASV
unique(contASV$ASV)
contaminantsN <- merge(contASV,filteredASV2N, by.x = "ASV", by.y = "ASV")
unique(contaminantsN$ASV)
sum(contaminantsN$taxa.average) #< 1.188126
levels(contaminantsN$d)
levels(contaminantsN$d) <- list("0"="0", "1"="1", "5"="5", "10"="10")

Y <- ggplot(contaminantsN, aes(x = d, y = taxa.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.5)) +
  labs(fill = "Only 100% similarity (Genus)") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Day')

Z <- ggplot(contaminantsN, aes(x = d, y = taxa.average, fill =SharedASVs)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.5)) +
  labs(fill = "Only 100% similarity (ASV)") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Day')

ggarrange(a,b,c,A,B,C,X,Y,Z, labels = c("a", "d","c","d","e","f","g", "h", "i"),vjust = 0.9,
          nrow = 3,ncol=3)

##check if the ASV found in E.coli that is 100% similar is found in the mock or only in the water
#E. coli is matchin the ASV1 in the controls
### CALCULATION OF THE ABUNDANCE OF EACH OTU  
ASV_tableC1 <- t(ASV_tableC)
otu.summary <- prop.table((ASV_tableC1), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:274)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)


#merging the abundance of each OTU with the metadata and the taxonomy file
controls <- read.csv("water_metadata.csv")
str(controls)
str(melt_otu)
meta_otu2 <- merge(controls, melt_otu, by.x = "ID", by.y = "Sample")
str(meta_otu2)
meta_otu_tax <- merge(meta_otu2, TaxASVC, by.x = "ASV", by.y = 0)
str(meta_otu_tax)
order_groups <- controls$ID
meta_otu_tax$Row.names <- factor(meta_otu_tax$ID, levels = order_groups)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$ASV <- factor(meta_otu_tax$ASV)
str(meta_otu_tax)


ASVcontrols <- meta_otu_tax %>% 
  group_by(Genus) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/9)*100) ## the total number of samples (51)
attach(ASVcontrols)
sum(ASVcontrols$Ave_Abundance)
ASVcontrols  <- ASVcontrols[order(-Ave_Abundance),]
write.csv(ASVcontrols,"ASVcontrols.csv")

ASVcontrols2 <- meta_otu_tax%>% 
  group_by(ID,Genus, ASV) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(ID, ASV) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 


e.coli <- subset(meta_otu_tax, Genus =="Escherichia-Shigella")
unique(e.coli$Genus)
unique(e.coli$ASVnos)

lacto <- subset(meta_otu_tax, Genus =="Lactobacillus")
unique(lacto$Genus)
unique(lacto$ASVnos)

clos <- subset(meta_otu_tax, Genus =="Clostridium sensu stricto 1")
unique(clos$Genus)
unique(clos$ASVnos)

e.coli2 <- subset(ASVcontrols2, ASV =="ASV1")
unique(e.coli2$ASV)

## relative abundance of all the ASVs
ASV <- e.coli %>% 
  group_by(ASVnos) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/9)*100) ## the total number of samples (51)
attach(ASV )

ASV <- e.coli %>% 
  group_by(ASVnos, ID) %>% 
  summarise(Ave_Abundance = (sum(Abundance))) ## the total number of samples (51)
attach(ASV )


e <- ggplot(ASV, aes(x = ID, y = Ave_Abundance, fill =ASVnos)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  labs(title = "Escherichia-Shigella (Controls)") +
  labs(fill = "ASV") +
  ylim(c(0,0.5)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Sample') +
  theme(axis.text.x = element_text(angle = 90))


#Lactobacillus
ASVL <- lacto %>% 
  group_by(ASVnos) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/9)*100) ## the total number of samples (51)
attach(ASVL)

ASVL1 <- lacto %>% 
  group_by(ASVnos, ID) %>% 
  summarise(Ave_Abundance = (sum(Abundance))) ## the total number of samples (51)
attach(ASV )


i <- ggplot(ASVL1, aes(x = ID, y = Ave_Abundance, fill =ASVnos)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  labs(title = "Lactobacillus (Controls)") +
  labs(fill = "ASV") +
  ylim(c(0,0.5)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Sample') +
  theme(axis.text.x = element_text(angle = 90))

ggarrange(e,i,o, labels = c("a", "b","c"),vjust = 0.9,
          ncol=3)

#Lactobacillus
ASVC <- clos %>% 
  group_by(ASVnos) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/9)*100) ## the total number of samples (51)
attach(ASVC)

ASVC1 <- clos %>% 
  group_by(ASVnos, ID) %>% 
  summarise(Ave_Abundance = (sum(Abundance))) ## the total number of samples (51)
attach(ASVC1)


o <- ggplot(ASVC1, aes(x = ID, y = Ave_Abundance, fill =ASVnos)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  labs(title = "Clostridium sensu stricto 1 (Controls)") +
  labs(fill = "ASV") +
  ylim(c(0,0.5)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Sample') +
  theme(axis.text.x = element_text(angle = 90))

###subset only for the Mock samples
str(ASV_tableC)
ASV_mock <- ASV_tableC[,c(2:5)]
MockTax <- merge(ASV_mock, TaxASVC, by.x = 0, by.y = 0) #merging asv_table with taxa
Mecoli <- subset(MockTax, Genus =="Escherichia-Shigella") #subsetting for ecoli
SeqMecoli <- merge(Mecoli, seqC, by.x = "OTUs", by.y = 0)
write.csv(SeqMecoli, "SeqMecoli.csv")

## Remove ASV
conta <- read.csv("Contaminants.csv") #subsetting to only include Lactobacillus, ES, MY, CSS, MH, PA, MA
conta1 <- merge(conta, seq1, by.x = "ASV", by.y = "ASVnos")
#write.csv(conta1, "conta1.csv")
conta1 <- read.csv("conta1.csv") # we just added a column with the OTU name and ">"
names <- c("ASV412")

rowSums(ASV_table)
ASV_table1 <- t(ASV_table)
colSums(ASV_table1)
ASV_table1 <- as.data.frame(ASV_table1)
colSums(ASV_table1)
ASV_table1 <- tibble::rownames_to_column(ASV_table1, "ASV")
filteredASV_table <- ASV_table1[!ASV_table1$ASV %in% names, ] # we should have 3859 ASVs
filteredASV_table <- as.data.frame(filteredASV_table )
filteredASV_table <- merge(TaxASV, filteredASV_table, by.x = 0, by.y = "ASV")
filteredASV_table1 <- filteredASV_table 
filteredASV_table <- filteredASV_table[,-c(1:8,10)]
write.csv(filteredASV_table, "filteredASV_table.csv")


#filter taxonomy
str(filteredASV_table1)
filtered_taxonomy <- filteredASV_table1[!filteredASV_table1$ASVnos %in% names, ]
filtered_taxonomy <- filtered_taxonomy[,-c(11:106)]
write.table(filtered_taxonomy, "filtered_taxonomy.txt")



#Function
phyloseq_to_df <- function(physeq, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance"){
  
  # require(phyloseq)
  
  ## Data validation
  if(any(addtax == TRUE || sorting == "taxonomy")){
    if(is.null(phyloseq::tax_table(physeq, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }
  }
  
  ## Prepare data frame
  if(taxa_are_rows(physeq) == TRUE){
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), phyloseq::otu_table(physeq), stringsAsFactors = F)
  } else {
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), t(phyloseq::otu_table(physeq)), stringsAsFactors = F)
  }
  
  ## Check if the sample names were silently corrected in the data.frame
  if(any(!phyloseq::sample_names(physeq) %in% colnames(res)[-1])){
    if(addtax == FALSE){
      warning("Warning: Sample names were converted to the syntactically valid column names in data.frame. See 'make.names'.\n")
    }
    
    if(addtax == TRUE){
      stop("Error: Sample names in 'physeq' could not be automatically converted to the syntactically valid column names in data.frame (see 'make.names'). Consider renaming with 'sample_names'.\n")
    }
  }
  
  ## Add taxonomy
  if(addtax == TRUE){
    
    ## Extract taxonomy table
    taxx <- as.data.frame(phyloseq::tax_table(physeq), stringsAsFactors = F)
    
    ## Reorder taxonomy table
    taxx <- taxx[match(x = res$OTU, table = rownames(taxx)), ]
    
    ## Add taxonomy table to the data
    res <- cbind(res, taxx)
    
    ## Add max tax rank column
    if(addmaxrank == TRUE){
      
      ## Determine the lowest level of taxonomic classification
      res$LowestTaxRank <- get_max_taxonomic_rank(taxx, return_rank_only = TRUE)
      
      ## Reorder columns (OTU name - Taxonomy - Max Rank - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), "LowestTaxRank", phyloseq::sample_names(physeq))]
      
    } else {
      ## Reorder columns (OTU name - Taxonomy - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), phyloseq::sample_names(physeq))]
      
    } # end of addmaxrank
  }   # end of addtax
  
  ## Reorder OTUs
  if(!is.null(sorting)){
    
    ## Sort by OTU abundance
    if(sorting == "abundance"){
      otus <- res[, which(colnames(res) %in% phyloseq::sample_names(physeq))]
      res <- res[order(rowSums(otus, na.rm = T), decreasing = T), ]
    }
    
    ## Sort by OTU taxonomy
    if(sorting == "taxonomy"){
      taxtbl <- as.data.frame( phyloseq::tax_table(physeq), stringsAsFactors = F )
      
      ## Reorder by all columns
      taxtbl <- taxtbl[do.call(order, taxtbl), ]
      # taxtbl <- data.table::setorderv(taxtbl, cols = colnames(taxtbl), na.last = T)
      res <- res[match(x = rownames(taxtbl), table = res$OTU), ]
    }
  }
  
  ## Add OTU total abundance
  if(addtot == TRUE){
    res$Total <- rowSums(res[, which(colnames(res) %in% phyloseq::sample_names(physeq))])
  }
  
  rownames(res) <- NULL
  return(res)
}


