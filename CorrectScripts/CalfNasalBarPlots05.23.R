#             R-Script for Calf Fecal microbiota Taxonomy Barplots
#                     Author: Ruth Eunice Centeno Martinez, 2020
#          Written at Dr. Tim Johnson Lab, Dept. of Animal Sciences,       #
#                             Purdue University, 2020
#                               rcenteno@purdue.edu
library(ggplot2)
library(tidyr) #separate function
library(reshape2) #melt function
library(dplyr)
library(naniar) # for replace_with_na_all function
library(data.table)
library(phyloseq)
library(qiime2R)
library(ggpubr)


setwd("~/Desktop/eunice/PhD/CalfNasal/Qiime/CorrectFiltered/")
##
rm(ASV_table)
##### Taxonomy barplot
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
# Read the experimental design, and species classified documents
ASVs <- read_qza("filteredASV_table.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #3859 ASVs
#str(ASV_table)
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
colSums(ASV_table)

#Importing the metadata file
metadata <- read.csv("CalfMetadata08.30.csv")
str(metadata)
rownames(metadata) <- metadata$ID
metadata$calf = as.factor(metadata$calf)
metadata$d <- as.numeric(metadata$d)
metadata$Sickness <- as.factor(metadata$Sickness)
IDCorrect <- metadata$ID

#merging the abundance of each OTU with the metadata and the taxonomy file
##Adding taxonomy
#Taxonomy of each OTU
tax <- read_qza("filteredTax.qza")
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
tax3 = tax3[,-c(1:2,10:13)]
tax.clean <- as.data.frame(tax3)
tax.clean$OTUs <- rownames(tax.clean)
#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ASV_s),]

#write.table(tax.final,"taxonomyNasal.txt",sep=",", row.names = FALSE) 
TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
row.names(TaxASV) <- TaxASV[,10]
TaxASV = TaxASV[,-c(1,10)]
# write.table(TaxASV,"TaxASV.txt",sep="\t", row.names = T, col.names = T)

#Preparing ASVtable
NewASVtable <- t(ASV_table)
NewASVtable <- merge(metadata, NewASVtable, by.x = 0, by.y = 0)
str(NewASVtable)
NewASVtable$ID <- as.factor(NewASVtable$ID)
levels(NewASVtable$ID) <- list("d0_1"="d0_1",    "d7_1"="d7_1",    "ds0_1"="ds0_1",   "d0_2"="d0_2",    "d7_2"="d7_2",    "d14_2"="d14_2",   "d0_3"="d0_3",    "d7_3"="d7_3",    "d14_3"="d14_3",  
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
NewASVtable <-NewASVtable[,-c(2:9)]
NewASVtable <- arrange(NewASVtable, Row.names)
rownames(NewASVtable) <- NewASVtable$Row.names
ASV_table <- NewASVtable[,-c(1)]
ASV_table <- as.matrix(ASV_table) #95 and 3859 

#---------------MIKING THE MODEL CON 
#generating a dataframe with all the response (days) from the samples
# Make one column for our outcome/response variable 
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq # [ 3859 taxa and 95 samples ]
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
physeq_deseq ##[ 3837 taxa and 95 samples ]
## Random forest, we want to make the model so it classify samples based on the age
ntaxa(physeq_deseq) #total of 3837

physeq_deseq
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
#write.table(NewTax,"NewTax.txt",sep=",", row.names = TRUE)
NewASVtable <- prunetable
NewASVtable <- NewASVtable[,-c(2:9)]
row.names(NewASVtable) <- NewASVtable[,1]
NewASVtable = t(NewASVtable[,-c(1)]) #3837 and 8 

#### Separate by healthy and sick
healthy <- subset(metadata, Sickness=="Never")
str(healthy)
healthy$d <- as.factor(healthy$d)
levels(healthy$d) <- list("0"="0", "7"="1", "14"="2")

H_OTU <- NewASVtable[rownames(NewASVtable) %in% rownames(healthy),] #51, 3837

## -----------------------------------PHYLUM LEVEL, Healthy period
### CALCULATION OF THE ABUNDANCE OF EACH OTU  
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
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:3837)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
str(metadata)
str(melt_otu)
meta_otu <- merge(healthy, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu)
meta_otu_tax <- merge(meta_otu, NewTax, by.x = "ASV", by.y = 0)
str(meta_otu_tax)
order_groups <- healthy$ID
meta_otu_tax$Row.names <- factor(meta_otu_tax$Row.names, levels = order_groups)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (3837 total, same value from the taxonomy file) 
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Sickness <- factor(meta_otu_tax$Sickness)
levels(meta_otu_tax$Sickness) <- list("Never"="Healthy", "Sick"="Got")
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$ASV <- factor(meta_otu_tax$ASV)
str(meta_otu_tax)

### Abundance at a phlyum level
sum(meta_otu_tax$Abundance)
phylum <- meta_otu_tax %>% 
  group_by(Phylum) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/51)) ## the total number of samples (51)
attach(phylum)
sum(phylum$Ave_Abundance)
phylum <- phylum[order(-Ave_Abundance),]
write.csv(phylum,"phylumHealthyNasal.csv")

###Mycoplasma abundance in the healthy animals
myco <- subset(meta_otu_tax, Genus =="Mycoplasma")
str(myco)
Mycoplasmas <- as.data.frame(myco %>% count(Species))

##check how mycoplasma species are different across time points 
mycoH <- myco %>% 
  group_by(Row.names, d, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(mycoH)
sum(mycoH$taxa.average) #0.4454338
mycoH$Species <- factor(mycoH$Species)
write.csv(mycoH,"MycoplasmaHealthyNasal.csv")

my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)

g <- ggplot(mycoH, aes(x = d, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.3)) +
  guides(fill="none") +
  #guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  #guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Mycoplasma ~ Healthy ") +
  ylab(paste0("Mycoplasma species")) +  labs(x='Day')

#abundance of mycoplasma 
lacto <- subset(meta_otu_tax, Genus =="Lactobacillus")
str(lacto)
lactobacillus <- as.data.frame(lacto %>% count(Species))

##check how mycoplasma species are different across time points 
lactoH <- lacto %>% 
  group_by(Row.names, d, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(lactoH)
sum(lactoH$taxa.average) #0.2238292
lactoH$Species <- factor(lactoH$Species)
write.csv(lactoH,"LactoHealthyNasal.csv")

a <- ggplot(lactoH, aes(x = d, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.3)) +
  guides(fill="none") +
  #guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  #guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Lactobacillus ~ Healthy") +
  ylab(paste0("Lactobacillus species")) +  labs(x='Day')


####------------------------taxonomy sick period 
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
D_OTU <- NewASVtable[rownames(NewASVtable) %in% rownames(Draxxin),] #22, 3837
str(D_OTU)
N_OTU <- NewASVtable[rownames(NewASVtable) %in% rownames(Nuflur),] #11, 3837
str(N_OTU)

## check the prop.table,
otu.summaryD <- prop.table((D_OTU), 1) 
str(otu.summaryD)
otu_abundD <- colSums(otu.summaryD)
otu_abund2D <- as.data.frame(otu_abundD)
otu.summaryD <- rbind(otu_abundD, otu.summaryD)
str(otu.summaryD)
otu.summary_sortedD <- otu.summaryD[,order(otu.summaryD[1,], decreasing = TRUE)]
str(otu.summary_sortedD)
melt_otuD <- reshape2::melt(otu.summary_sortedD[, c(1:3837)]) ###TOTAL NUMBER OF OTUS
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

#abundance of mycoplasma 
mycoD <- subset(meta_otu_taxD, Genus =="Mycoplasma")
str(mycoD)
MycoplasmasD <- as.data.frame(mycoD %>% count(Species))

##check how mycoplasma species are different across time points 
mycoAD <- mycoD %>% 
  group_by(Row.names, d, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(mycoAD)
sum(mycoAD$taxa.average) #0.2694382
mycoAD$Species <- factor(mycoAD$Species)
write.csv(mycoAD,"MycoplasmaDraxxinNasal.csv")

h <- ggplot(mycoAD, aes(x = d, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.3)) +
  guides(fill="none") +
  #guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  #guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ggtitle("Mycoplasma ~ Tulathromycin") +
  ylab(paste0("Mycoplasma species")) +  labs(x='Day')

ggarrange(f,g, labels = c("a", "b"),vjust = 0.9,
          ncol = 2)
#abundance of mycoplasma 
lactoD <- subset(meta_otu_taxD, Genus =="Lactobacillus")
str(lactoD)
lactobacillusD <- as.data.frame(lactoD %>% count(Species))

##check how mycoplasma species are different across time points 
lactoAD <- lactoD %>% 
  group_by(Row.names, d, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(lactoAD)
sum(lactoAD$taxa.average) #0.5831663
lactoAD$Species <- factor(lactoAD$Species)
write.csv(lactoAD,"LactoDraxxinNasalTreatment.csv")

b <- ggplot(lactoAD, aes(x = d, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.3)) +
  guides(fill="none") +
  #guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  #guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ggtitle("Lactobacillus ~ Tulathromycin") +
  ylab(paste0("Lactobacillus species")) +  labs(x='Day')


#### Draxxin 
str(meta_otu_taxD)
meta_otu_taxD$calf <- as.factor(meta_otu_taxD$calf)
levels(meta_otu_taxD$calf)
phylumD <- meta_otu_taxD %>% 
  group_by(Phylum) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/22)*100) ## the total number of samples (21)
attach(phylumD)
sum(phylumD$Ave_Abundance)
phylumD <- phylumD[order(-Ave_Abundance),]
write.csv(phylumD,"phylumDraxxinNasal.csv")


##3--------------------Nuflur
## check the prop.table, it gives me NA in the sample Dairy 92
otu.summaryN <- prop.table((N_OTU), 1) 
str(otu.summaryN)
otu_abundN <- colSums(otu.summaryN)
otu_abund2D <- as.data.frame(otu_abundN)
otu.summaryN <- rbind(otu_abundN, otu.summaryN)
str(otu.summaryN)
otu.summary_sortedN <- otu.summaryN[,order(otu.summaryN[1,], decreasing = TRUE)]
str(otu.summary_sortedN)
melt_otuN <- reshape2::melt(otu.summary_sortedN[, c(1:3837)]) ###TOTAL NUMBER OF OTUS
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

#abundance of mycoplasma 
mycoN <- subset(meta_otu_taxN, Genus =="Mycoplasma")
str(mycoN)
MycoplasmasN <- as.data.frame(mycoN %>% count(Species))

##check how mycoplasma species are different across time points 
mycoAN <- mycoN %>% 
  group_by(Row.names, d, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(mycoAN)
sum(mycoAN$taxa.average) #0.05900168
mycoAN$Species <- factor(mycoAN$Species)
write.csv(mycoAN,"MycoplasmaNuflurNasal.csv")

i <- ggplot(mycoAN, aes(x = d, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.3)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Mycoplasma ~ Florfenicol") +
  ylab(paste0("Mycoplasma species")) +  labs(x='Day')

#abundance of mycoplasma 
lactoN <- subset(meta_otu_taxN, Genus =="Lactobacillus")
str(lactoN)
write.csv()
lactobacillusN <- as.data.frame(lactoN %>% count(Species))

##check how mycoplasma species are different across time points 
lactoAN <- lactoN %>% 
  group_by(Row.names,  d, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(lactoAN)
sum(lactoAN$taxa.average) #0.6399811
lactoAN$Species <- factor(lactoAN$Species)
write.csv(lactoAN,"LactoNuflurNasal.csv")

c <- ggplot(lactoAN, aes(x = d, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.3)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Lactobacillus ~ Florfenicol") +
  ylab(paste0("Lactobacillus species")) +  labs(x='Day')

ggarrange(i,c, labels = c("c", "f"),
          nrow=2, ncol = 1)

###now its time to identify how the lactobacillus is distributed based on the new taxonomy
#we need to remove the ASV identified just as Lactobacillus to include the new taxonomy
str(lacto)
lacto$Species = as.factor(lacto$Species)
levels(lacto$Species)
lacto<-lacto[!(lacto$Species=="Lactobacillus"),]
levels(lacto$Species)

lactoD$Species = as.factor(lactoD$Species)
levels(lactoD$Species)
lactoD<-lactoD[!(lactoD$Species=="Lactobacillus"),] 
levels(lactoD$Species)
lactoD = lactoD[,-c(4,10,12:17)]

lactoN$Species = as.factor(lactoN$Species)
levels(lactoN$Species)
lactoN<-lactoN[!(lactoN$Species=="Lactobacillus"),] 
levels(lactoN$Species)
lactoN = lactoN[,-c(4,10,12:17)]

## now we include the new taxonomy
lacTaxNew <- read.csv("CorrectLactoTaxonomyUnclassified.csv")
rownames(lacTaxNew) = lacTaxNew$Feature.ID

lacTaxH = merge(meta_otu_tax[,-c(12:18)], lacTaxNew, by.x = "OTUs", by.y = 0)
lacTaxH = lacTaxH[,-c(1,13,14)] #to make sure the column number is the same
str(lacTaxH)
str(lacto)
lacto = lacto[,-c(19)]
lacTaxH = rbind(lacto, lacTaxH)
levels(lacTaxH$Species)


lacTaxD = merge(meta_otu_taxD[,-c(4,10,12:16,18:24)], lacTaxNew, by.x = "OTUs", by.y = 0)
str(lacTaxD)
str(lactoD)
lacTaxD = lacTaxD[,-c(1,12,13)] #to make sure the column number is the same
lacTaxD = rbind(lactoD, lacTaxD)
levels(lacTaxD$Species)

lacTaxN = merge(meta_otu_taxN[,-c(4,10,12:16,18:24)], lacTaxNew, by.x = "OTUs", by.y = 0)
lacTaxN = lacTaxN[,-c(1,12,13)]
lacTaxN = rbind(lactoN, lacTaxN)
str(lacTaxN)
str(lactoN)


##now calculate abundance
lactoH1 <-lacTaxH %>% 
  group_by(ID,  d, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(lactoH1)
sum(lactoH1$taxa.average) # 0.2238292
lactoH1$Species <- factor(lactoH1$Species)
write.csv(lactoH1,"LactoHealthySpeciesNasal.csv")

x <- ggplot(lactoH1, aes(x = d, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.3)) +
  guides(fill="none")+
  #guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  #guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Lactobacillus sp ~ Healthy") +
  ylab(paste0("Lactobacillus species")) +  labs(x='Day')

lactoD1 <-lacTaxD %>% 
  group_by(ID,  d, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(lactoD1)
sum(lactoD1$taxa.average) # 0.5831663
lactoD1$Species <- factor(lactoD1$Species)
write.csv(lactoD1,"LactoDraxxinSpeciesNasal.csv")

y <- ggplot(lactoD1, aes(x = d, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.3)) +
  guides(fill="none")+
  #guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  #guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Lactobacillus sp ~ Tulathromycin") +
  ylab(paste0("Lactobacillus species")) +  labs(x='Day')

lactoN1 <-lacTaxN %>% 
  group_by(ID,  d, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(lactoN1)
sum(lactoN1$taxa.average) # 0.6399811
lactoN1$Species <- factor(lactoN1$Species)
write.csv(lactoN1,"LactoNuflurSpeciesNasal.csv")

z <- ggplot(lactoN1, aes(x = d, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,0.3)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Lactobacillus sp ~ Florfenicol") +
  ylab(paste0("Lactobacillus species")) +  labs(x='Day')


ggarrange(x,y, labels = c("a", "b"),
          nrow=1, ncol = 2)

#### Nuflur
str(meta_otu_taxN)
phylumN <- meta_otu_taxN %>% 
  group_by(Phylum) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/11)*100) ## the total number of samples (51)
attach(phylumN)
sum(phylumN$Ave_Abundance)
phylumN <- phylumN[order(-Ave_Abundance),]
write.csv(phylumN,"phylumNuflurNasal.csv")

top10phylumN <- phylumN[c(1:10),] #to select the top 10 most abundant phylum
sum(top10phylumN$Ave_Abundance) #the total of the community composed of the top 10
top10PN <- merge(meta_otu_taxN, top10phylumN, by.x = "Phylum", by.y="Phylum")
top10PN$Phylum <- as.factor(top10PN$Phylum)
levels(top10PN$Phylum)
str(top10PN)

phylumAN <- top10PN %>% 
  group_by(Row.names, d, Phylum) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Phylum) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(phylumAN)
sum(phylumAN$taxa.average)
phylumAN$Phylum <- factor(phylumAN$Phylum)
write.csv(phylumAN,"phylumNuflurNasalTop10.csv")

##### We then check the correlation between Psychobacter and Mycoplasma species using a repeated measures correlation
library(rmcorr)

str(meta_otu_tax)
corrH <- subset(meta_otu_tax, subset = Genus %in% c("Psychrobacter", "Mycoplasma")) #healthy
corrD <- subset(meta_otu_taxD, subset = Genus %in% c("Psychrobacter", "Mycoplasma")) #Draxxin
corrN <- subset(meta_otu_taxN, subset = Genus %in% c("Psychrobacter", "Mycoplasma")) #Nuflur
write.csv(corrH, "corrH.csv")

#663 + 633
str(corrH)

corrH <- read.csv("corrD.csv")
corrH$calf <- as.factor(corrH$calf)
corrH$d <- as.factor(corrH$d)
CH <- rmcorr(calf, Psychrobacter, Mycoplasma, corrH) #no effect including all the days
#subset by day
CH0 <- subset(corrH, d=="10")
library("ggpubr")
ggscatter(CH0, x = "Psychrobacter", y = "Mycoplasma", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Psychrobacter rel abun", ylab = "Mycoplasma Rel Abun") 

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
