library(qiime2R)
library(phyloseq)
library(zoo)
library(ggsci)
library(naniar) 
library(tidyr)
library("lmtest")
library(dplyr)
library(broom)
library(ggpubr)
library(tidyverse)

#install.packages("R2admb")
#install.packages("glmmADMB", 
library(glmmADMB)
library(broom.mixed)
library(lme4)
library("DESeq2")

setwd("~/Desktop/eunice/PhD/CalfNasal/Qiime/CorrectFiltered/")

#Importing the metadata file
metadata <- read.csv("CalfMetadata08.30.csv")
str(metadata)
rownames(metadata) <- metadata$ID
metadata$calf = as.factor(metadata$calf)
metadata$d <- as.numeric(metadata$d)
metadata$d <- as.numeric(metadata$d)
metadata$Sickness <- as.factor(metadata$Sickness)
IDCorrect <- metadata$ID

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

# Phylum level beta regression
tax<-read_qza("filteredTax.qza") #phylum level
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
NewASVtable = NewASVtable[,-c(1)] #3837 and 8 
#write.csv(NewASVtable,"NewASVtable.csv")

#----------------------BETA REGRESSION---------------------------------

##importing the new taxonomy
taxonomy<-read_qza("table-level1.qza") #phylum level
head(taxonomy$data)

taxa <- data.frame(t(taxonomy$data))
str(taxa)

## we need to filter the data based on sick and healthy animals
rm(healthy)
#I'm using the d as a continuous variable
healthy <- subset(metadata, Sickness=="Never") ### 51 samples
healthy$sample <- rownames(healthy)
str(healthy)
healthy$calf <- as.factor(healthy$calf)
rm(calf)
calf <- healthy$calf

###combining the healthy metadata with the taxa (ASV_table)
H_OTUP <- as.data.frame(taxa[rownames(taxa) %in% rownames(healthy),]) #51, 37 phylum
H_OTU <- t(H_OTUP)
table <- as.data.frame(H_OTU)
sum(table$d0_2) #1661 sequences
H_OTU <- prop.table((H_OTU), 2) #to calculate relative abundance for each sample but not grouped
colSums(H_OTU)
H_OTU <- H_OTU[apply(H_OTU, 1, function(x) mean(x) > 0.0001),] #from 37 to 27
colSums(H_OTU)
names <- healthy$sample
taxa2 <- cbind(names, H_OTU)
taxa2 <- taxa2[,-c(1)]
taxa3<- reshape2::melt(taxa2[, c(1:51)])
str(taxa3)
taxa3$value <- as.numeric(taxa3$value)
sum(taxa3$value)
colnames(taxa3) <- c("Phylum", "sample", "RA")
sample <- subset(taxa3, sample == "d0_2")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 0.9957857

##now mergin the summarized phylum table with the metadata
phylumH <- healthy %>% 
  inner_join(taxa3, by = "sample")
phylumH$Phylum <- as.character(phylumH$Phylum)
phylumH$RA <- as.numeric(phylumH$RA)
#phylumH$d <- factor(phylumH$d, levels = c("0","7","14"))
summary(phylumH)

healthyC <- phylumH[,-c(3:7)] # keep only the variables you mentioned 1377
healthyC <- filter(healthyC, RA != 0) #removing 0's 490

##now running the loop
str(healthyC)
healthyC$Phylum <- as.factor(healthyC$Phylum)
safe_betareg <- possibly(glmmadmb, NA_real_)
regH2 <-healthyC %>% 
  group_by(Phylum) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

select <- regH2
select[which(select$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
select[which(select$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
select <- select[order(select$term), ]

select_betareg = select[-c(1:32,49:91),]
select_betareg <- select_betareg[order(select_betareg$estimate), ]


A <- ggplot(select_betareg, aes(x=reorder(Phylum, estimate), y=estimate)) +
  geom_segment(aes(x = reorder(Phylum, estimate), xend =Phylum, y = 0, yend = estimate)) +
  geom_point(aes(shape = sig), size = 3) +
  labs(title = "Healthy Samples") +
  scale_shape_manual(values = c(16, 1)) +
  coord_flip() +
  labs(x='', y='Regression Coefficient')+
  ylim(-0.55, 0.5) +
  geom_hline(yintercept=0) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.title.x = element_text(color="black", size=10, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 8), axis.text.y = element_text(color = "black", size = 10, face ="italic")) 

### genus

##adding the taxonomy file ~genus
taxonomy2<-read_qza("table-level5.qza")
head(taxonomy2$data)
taxaG <- data.frame(taxonomy2$data)
taxaG$ASVnos <- paste0("ASV", 1:nrow(taxaG))
taxaG$ASVstring <- rownames(taxaG)
rownames(taxaG) <- taxaG$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- taxaG[, (ncol(taxaG)-1):ncol(taxaG)] #the key withe the names
rownames(ASVkey) <- ASVkey$ASVstring
ASVkey  <- cbind(ASVkey , read.table(text=row.names(ASVkey), sep=";", 
                                     header=FALSE, col.names = paste0("col", 1:5), stringsAsFactors=FALSE))
taxaG <- taxaG[,-(ncol(taxaG)-1):-ncol(taxaG)]

taxaG <- t(taxaG)
taxaG <- data.frame(taxaG)

###combining the healthy metadata with the taxa (ASV_table)
H_OTUGP <- as.data.frame(taxaG[rownames(taxaG) %in% rownames(healthy),]) #916
H_OTUG <- t(H_OTUGP)
table <- as.data.frame(H_OTUG)
sum(table$d0_2) #1661 sequences
H_OTUG <- prop.table((H_OTUG), 2) #to calculate relative abundance for each sample but not grouped
colSums(H_OTUG)
H_OTUG <- H_OTUG[apply(H_OTUG, 1, function(x) mean(x) > 0.001),] #filter the ASVs with an abundance < than 0.001, from 916 to 110
colSums(H_OTUG)
names <- healthy$sample
taxaG2 <- cbind(names, H_OTUG)
taxaG2 <- taxaG2[,-c(1)]
taxaG3<- reshape2::melt(taxaG2[, c(1:51)])
str(taxaG3)
taxaG3$value <- as.numeric(taxaG3$value)
sum(taxaG3$value)
colnames(taxaG3) <- c("Genus", "sample", "RA")
sample <- subset(taxaG3, sample == "d0_2")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 0.8531005


##now mergin the summarized phylum table with the metadata
genus.C <- healthy %>% 
  inner_join(taxaG3, by = "sample")
summary(genus.C)
genus.C$Genus <- as.character(genus.C$Genus)
genus.C$RA <- as.numeric(genus.C$RA)

####separate by treatment CON and SCFP
healthyCG <- genus.C[,-c(3:7)] # keep only the variables you mentioned 5610
healthyCG <- filter(healthyCG, RA != 0) #removing 0's 2178

##now running the loop
str(healthyCG)
healthyCG$Genus <- as.factor(healthyCG$Genus)
regHG <-healthyCG %>% 
  group_by(Genus) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectG <- regHG
selectG[which(selectG$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectG[which(selectG$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectG <- selectG[order(selectG$term), ]

select_betaregG = selectG[-c(1:160,241:430),]
select_betaregG <- select_betaregG[order(select_betaregG$estimate), ]

select_betaregSig <- filter(select_betaregG, sig=="P-Value < 0.05")
select_betaregSig <- select_betaregSig[-2]
#merge taxonomy with the significant Gens
SigGHealthyNT <- merge(select_betaregSig, ASVkey, by.x = "Genus", by.y = "ASVnos")
write.csv(SigGHealthyNT, "SigGHealthyNT.csv")
SigGHealthyNT <- read.csv("SigGHealthyNT.csv")
SigGHealthyNT <- SigGHealthyNT[order(SigGHealthyNT$estimate), ]

my_colors <- c(
 '#1f78b4','#b3df8a',"coral","#653936","lightpink",
  '#fdbf6f','#ff7f00','#cab3d6','#6a3d9a','#ffff99','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

a <- ggplot(SigGHealthyNT, aes(x=reorder(genus, estimate), y=estimate, color=col1)) +
  geom_segment(aes(x = reorder(genus, estimate), xend =genus, y = 0, yend = estimate)) +
  geom_point(aes(shape = sig), size = 3) +
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = c(16, 1)) +
  labs(title = "Genus") +
  labs(x='', y='Regression Coefficient')+
  theme_bw() +
  labs(color= "Phylum") +
  ylim(-0.55, 0.5) +
  geom_hline(yintercept=0) +
  coord_flip() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.title.x = element_text(color="black", size=10, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 9), axis.text.y = element_text(color = "black", size = 10, face ="italic")) 
a

##adding the taxonomy file ~species
taxonomy3<-read_qza("table-level6.qza")
head(taxonomy3$data)
taxaS <- data.frame(taxonomy3$data)
taxaS$ASVnos <- paste0("ASV", 1:nrow(taxaS))
taxaS$ASVstring <- rownames(taxaS)
rownames(taxaS) <- taxaS$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey2 <- taxaS[, (ncol(taxaS)-1):ncol(taxaS)] #the key withe the names
rownames(ASVkey2) <- ASVkey2$ASVstring
ASVkey2  <- cbind(ASVkey2 , read.table(text=row.names(ASVkey2), sep=";", 
                                     header=FALSE, col.names = paste0("col", 1:6), stringsAsFactors=FALSE))
taxaS <- taxaS[,-(ncol(taxaS)-1):-ncol(taxaS)]

taxaS <- t(taxaS)
taxaS <- data.frame(taxaS)

###combining the healthy metadata with the taxa (ASV_table)
H_OTUS <- as.data.frame(taxaS[rownames(taxaS) %in% rownames(healthy),]) #1507
H_OTUS <- t(H_OTUS)
table <- as.data.frame(H_OTUS)
sum(table$d0_2) #1661 sequences
H_OTUS <- prop.table((H_OTUS), 2) #to calculate relative abundance for each sample but not grouped
colSums(H_OTUS)
H_OTUSpecies = H_OTUS
H_OTUS <- H_OTUS[apply(H_OTUS, 1, function(x) mean(x) > 0.001),] #filter the ASVs with an abundance < than 0.001, from 1507 to 125
colSums(H_OTUS)
names <- healthy$sample
taxaS2 <- cbind(names, H_OTUS)
taxaS2 <- taxaS2[,-c(1)]
taxaS3<- reshape2::melt(taxaS2[, c(1:51)])
str(taxaS3)
taxaS3$value <- as.numeric(taxaS3$value)
sum(taxaS3$value)
colnames(taxaS3) <- c("Specie", "sample", "RA")
sample <- subset(taxaS3, sample == "d0_2")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 0.8007225


##now mergin the summarized phylum table with the metadata
species.C <- healthy %>% 
  inner_join(taxaS3, by = "sample")
summary(species.C)
species.C$Specie <- as.character(species.C$Specie)
species.C$RA <- as.numeric(species.C$RA)

healthyCS <- species.C[,-c(3:7)] # keep only the variables you mentioned 6375
healthyCS <- filter(healthyCS, RA != 0) #removing 0's 2161

##now running the loop
str(healthyCS)
healthyCS$Specie <- as.factor(healthyCS$Specie)
regHS <-healthyCS %>% 
  group_by(Specie) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectS <- regHS
selectS[which(selectS$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectS[which(selectS$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectS <- selectS[order(selectS$term), ]

select_betaregS = selectS[-c(1:196,295:517),]
select_betaregS <- select_betaregS[order(select_betaregS$estimate), ]

select_betaregSigS <- filter(select_betaregS, sig=="P-Value < 0.05")
#merge taxonomy with the significant Gens
SigGHealthyNS <- merge(select_betaregSigS, ASVkey2, by.x = "Specie", by.y = "ASVnos") %>% 
  unnest(data)
write.csv(SigGHealthyNS, "SigGHealthyNS.csv")
SigGHealthyNS <- read.csv("SigGHealthyNS.csv")
SigGHealthyNS <- SigGHealthyNS[order(SigGHealthyNS$estimate), ]


b <- ggplot(SigGHealthyNS, aes(x=reorder(species, estimate), y=estimate, color=col1)) +
  geom_segment(aes(x = reorder(species, estimate), xend =species, y = 0, yend = estimate)) +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values = my_colors) +
  labs(title = "Species") +
  labs(x='', y='Regression Coefficient')+
  ylim(-0.5, 0.4) +
  geom_hline(yintercept=0) +
  theme_bw() +
  labs(color= "Phylum") +
  coord_flip() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.title.x = element_text(color="black", size=10, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 9), axis.text.y = element_text(color = "black", size = 10, face ="italic")) 
b

ggarrange(A,a,b, labels = c("a", "b", "c"),
          ncol=2, nrow=2)
A
a
b

##identify if the BRD pathogens significantly change through time
#only include the pathogens species, work with ASVKey2 
Hs = subset(ASVkey2, col6 == "Histophilus")
Pm = subset(ASVkey2, col6 == "Pasteurella multocida")
Mh = subset(ASVkey2, col6 == "Mannheimia")
Mb = subset(ASVkey2, col6 == "Mycoplasma bovis")
patho = rbind(Hs, Pm, Mh, Mb)
rownames(patho) <- patho$ASVnos

taxBRD <- H_OTUSpecies[rownames(H_OTUSpecies) %in% rownames(patho),]
names <- healthy$sample
taxaS2.1 <- cbind(names, taxBRD)
taxaS2.1 <- taxaS2.1[,-c(1)]
taxaS3.1<- reshape2::melt(taxaS2.1[, c(1:51)])
str(taxaS3.1)
taxaS3.1$value <- as.numeric(taxaS3.1$value)
sum(taxaS3.1$value)
colnames(taxaS3.1) <- c("Specie", "sample", "RA")
sample <- subset(taxaS3.1, sample == "d0_11")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 0.0006157635

##now mergin the summarized phylum table with the metadata
species.C1 <- healthy %>% 
  inner_join(taxaS3.1, by = "sample")
summary(species.C1)
species.C1$Specie <- as.character(species.C1$Specie)
species.C1$RA <- as.numeric(species.C1$RA)

####separate by treatment CON and SCFP
healthyCS1 <- species.C1[,-c(3:7)] # keep only the variables you mentioned 6375
healthyCS1 <- filter(healthyCS1, RA != 0) #removing 0's 18

##now running the loop
str(healthyCS1)
healthyCS1$Specie <- as.factor(healthyCS1$Specie)
regHS <-healthyCS1 %>% 
  group_by(Specie) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectS <- regHS
selectS[which(selectS$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectS[which(selectS$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectS <- selectS[order(selectS$term), ]

####----------------------------------------now we go with the sick period
ant <- read.table("CalfMetadataSickAnimals.txt", header=TRUE, sep="\t")
ant <- ant[,-c(11:14)]
#also using d as continuous factoe
rownames(ant) <- ant$ID2
ant$antibiotic <- as.factor(ant$antibiotic)
ant$sample <- ant$ID2

Draxxin <- subset(ant, antibiotic=="Draxxin") #22 samples
Draxxin$sample <- Draxxin$ID2
rownames(Draxxin) <- Draxxin$ID2
Draxxin$calf <- as.factor(Draxxin$calf)
calf <- Draxxin$calf #6 animals

Nuflur <- subset(ant, antibiotic=="Nuflur") #11 samples
rownames(Nuflur) <- Nuflur$ID2
Nuflur$calf <- as.factor(Nuflur$calf)
calf <- Nuflur$calf #3 animals

ant$calf <- as.factor(ant$calf)
calf <- ant$calf #10 animals

###creating the asv_table 
D_OTUP <- as.data.frame(taxa[rownames(taxa) %in% rownames(Draxxin),]) #51, 37 phylum
D_OTU <- t(D_OTUP)
table <- as.data.frame(D_OTU)
sum(table$ds0_9) #1889 sequences
D_OTU <- prop.table((D_OTU), 2) #to calculate relative abundance for each sample but not grouped
colSums(D_OTU)
D_OTU <- D_OTU[apply(D_OTU, 1, function(x) mean(x) > 0.0001),] #from 37 to 22
colSums(D_OTU)
names <- Draxxin$sample
taxaD2 <- cbind(names, D_OTU)
taxaD2 <- taxaD2[,-c(1)]
taxaD3<- reshape2::melt(taxaD2[, c(1:22)])
str(taxaD3)
taxaD3$value <- as.numeric(taxaD3$value)
sum(taxaD3$value)
colnames(taxaD3) <- c("Phylum", "sample", "RA")
sample <- subset(taxaD3, sample == "ds0_9")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts=1

##now mergin the summarized phylum table with the metadata
phylumD <- Draxxin %>% 
  inner_join(taxaD3, by = "sample")
summary(phylumD)
phylumD$Phylum <- as.character(phylumD$Phylum)
phylumD$RA <- as.numeric(phylumD$RA)

phylumD <- phylumD[,-c(2:8,10)] # keep only the variables you mentioned 550
phylumD  <- filter(phylumD , RA != 0) #removing 0's 234
##now running the loop
str(phylumD)
phylumD$Phylum <- as.factor(phylumD$Phylum)
regD <-phylumD %>% 
  group_by(Phylum) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectD <- regD
selectD[which(selectD$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectD[which(selectD$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectD <- selectD[order(selectD$term), ]

select_betaregD = selectD[-c(1:36,55:97),]
select_betaregD <- select_betaregD[order(select_betaregD$estimate), ]
select_betaregD <- select_betaregD[-2]

X <- ggplot(select_betaregD, aes(x=reorder(Phylum, estimate), y=estimate)) +
  geom_segment(aes(x = reorder(Phylum, estimate), xend =Phylum, y = 0, yend = estimate)) +
  geom_point(aes(shape = sig), size = 3) +
  labs(title = "Tulathromycin - Phylum") +
  scale_shape_manual(values = c(16, 1)) +
  coord_flip() +
  ylim(-0.5, 0.5) +
  geom_hline(yintercept=0) +
  labs(x='', y='Regression Coefficient')+
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.title.x = element_text(color="black", size=10, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 8), axis.text.y = element_text(color = "black", size = 10, face ="italic")) 
X

###creating the asv_table for Nuflur
N_OTUP <- as.data.frame(taxa[rownames(taxa) %in% rownames(Nuflur),]) #51, 37 phylum
N_OTU<- t(N_OTUP)
table <- as.data.frame(N_OTU)
sum(table$ds1_15) #1736 sequences
N_OTU<- prop.table((N_OTU), 2) #to calculate relative abundance for each sample but not grouped
colSums(N_OTU)
N_OTU<- N_OTU[apply(N_OTU, 1, function(x) mean(x) > 0.0001),] #from 37 to 19
colSums(N_OTU)
names <- Nuflur$sample
taxaN2 <- cbind(names, N_OTU)
taxaN2 <- taxaN2[,-c(1)]
taxaN3<- reshape2::melt(taxaN2[, c(1:10)])
str(taxaN3)
taxaN3$value <- as.numeric(taxaN3$value)
sum(taxaN3$value)
colnames(taxaN3) <- c("Phylum", "sample", "RA")
sample <- subset(taxaN3, sample == "ds1_15")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 1

##now mergin the summarized phylum table with the metadata
phylumN <- Nuflur %>% 
  inner_join(taxaN3, by = "sample")
summary(phylumN)
phylumN$Phylum <- as.character(phylumN$Phylum)
phylumN$RA <- as.numeric(phylumN$RA)

phylumN <- phylumN[,-c(2:8,10)] # keep only the variables you mentioned 190
phylumN  <- filter(phylumN , RA != 0) #removing 0's 91

##now running the loop
str(phylumN)
phylumN$Phylum <- as.factor(phylumN$Phylum)
regN <-phylumN %>% 
  group_by(Phylum) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectN <- regN
selectN[which(selectN$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectN[which(selectN$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectN <- selectN[order(selectN$term), ]

select_betaregN = selectN[-c(1:24,37:67),]
select_betaregN <- select_betaregN[order(select_betaregN$estimate), ]
select_betaregN <- select_betaregN[-2]

Y <- ggplot(select_betaregN, aes(x=reorder(Phylum, estimate), y=estimate)) +
  geom_segment(aes(x = reorder(Phylum, estimate), xend =Phylum, y = 0, yend = estimate)) +
  geom_point(aes(shape = sig), size = 3) +
  labs(title = "Florfenicol - Phylum") +
  scale_shape_manual(values = c(16, 1)) +
  coord_flip() +
  labs(x='', y='Regression Coefficient')+
  ylim(-0.5, 0.5) +
  geom_hline(yintercept=0) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.title.x = element_text(color="black", size=10, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 8), axis.text.y = element_text(color = "black", size = 10, face ="italic")) 
Y

###### Draxxin and Nuflur genus
D_OTUGP <- as.data.frame(taxaG[rownames(taxaG) %in% rownames(Draxxin),]) #916 genera
D_OTUG <- t(D_OTUGP)
table <- as.data.frame(D_OTUG)
sum(table$ds0_9) #1739 sequences
D_OTUG <- prop.table((D_OTUG), 2) #to calculate relative abundance for each sample but not grouped
colSums(D_OTUG)
D_OTUG <- D_OTUG[apply(D_OTUG, 1, function(x) mean(x) > 0.001),] #filter the ASVs with an abundance < than 0.001, from 915 to 97
colSums(D_OTUG)
names <- Draxxin$sample
taxaDG2 <- cbind(names, D_OTUG)
taxaDG2 <- taxaDG2[,-c(1)]
taxaDG3<- reshape2::melt(taxaDG2[, c(1:22)])
str(taxaDG3)
taxaDG3$value <- as.numeric(taxaDG3$value)
sum(taxaDG3$value)
colnames(taxaDG3) <- c("Genus", "sample", "RA")
sample <- subset(taxaDG3, sample == "ds0_9")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 0.6992524

##now mergin the summarized phylum table with the metadata
genus.DC <- Draxxin %>% 
  inner_join(taxaDG3, by = "sample")
summary(genus.DC)
genus.DC$Genus <- as.character(genus.DC$Genus)
genus.DC$RA <- as.numeric(genus.DC$RA)

####separate by treatment CON and SCFP
genus.DC <- genus.DC[,-c(2:8, 10)] # keep only the variables you mentioned 2134
genus.DC <- filter(genus.DC, RA != 0) #removing 0's 874

##now running the loop
str(genus.DC)
genus.DC$Genus <- as.factor(genus.DC$Genus)
regDG <-genus.DC %>% 
  group_by(Genus) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectDG <- regDG
selectDG[which(selectDG$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectDG[which(selectDG$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectDG <- selectDG[order(selectDG$term), ]

select_betaregDG = selectDG[-c(1:144,217:385),]
select_betaregDG <- select_betaregDG[order(select_betaregDG$estimate), ]
select_betaregDG <- select_betaregDG[-2]
select_betaregSigDG <- filter(select_betaregDG, sig=="P-Value < 0.05")

#merge taxonomy with the significant Gens
SigGDraxxin <- merge(select_betaregSigDG, ASVkey, by.x = "Genus", by.y = "ASVnos") 
write.csv(SigGDraxxin, "SigGDraxxin.csv")
SigGDraxxin <- read.csv("SigGDraxxin.csv")
SigGDraxxin <- SigGDraxxin[order(SigGDraxxin$estimate), ]

my_colors1 <- c(
  '#1f78b4','#b3df8a',"#8A7C64","coral","lightpink",
  '#fdbf6f','#a6cee3','#cab3d6','#6a3d9a','#ffff99','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

x <- ggplot(SigGDraxxin, aes(x=reorder(genus, estimate), y=estimate, color=col1)) +
  geom_segment(aes(x = reorder(genus, estimate), xend =genus, y = 0, yend = estimate)) +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values = my_colors1) +
  labs(title = "Genus") +
  ylim(-0.5, 0.7) +
  geom_hline(yintercept=0) +
  labs(x='', y='Regression Coefficient')+
  theme_bw() +
  labs(color= "Phylum") +
  coord_flip() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.title.x = element_text(color="black", size=10, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 9), axis.text.y = element_text(color = "black", size = 10, face ="italic")) 
x

##Species Draxxin
D_OTUSP <- as.data.frame(taxaS[rownames(taxaS) %in% rownames(Draxxin),]) #1507 genera
D_OTUS <- t(D_OTUSP)
table <- as.data.frame(D_OTUS)
sum(table$ds0_9) #1739 sequences
D_OTUS <- prop.table((D_OTUS), 2) #to calculate relative abundance for each sample but not grouped
colSums(D_OTUS)
D_OTUSpecies = D_OTUS
D_OTUS <- D_OTUS[apply(D_OTUS, 1, function(x) mean(x) > 0.01),] #filter the ASVs with an abundance < than 0.001, from 1507 to 21
colSums(D_OTUS)
names <- Draxxin$sample
taxaDS2 <- cbind(names, D_OTUS)
taxaDS2 <- taxaDS2[,-c(1)]
taxaDS3<- reshape2::melt(taxaDS2[, c(1:22)])
str(taxaDS3)
taxaDS3$value <- as.numeric(taxaDS3$value)
sum(taxaDS3$value)
colnames(taxaDS3) <- c("Species", "sample", "RA")
sample <- subset(taxaDS3, sample == "ds0_9")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 0.2455434

##now mergin the summarized phylum table with the metadata
species.DC <- Draxxin %>% 
  inner_join(taxaDS3, by = "sample")
summary(species.DC)
species.DC$Species <- as.character(species.DC$Species)
species.DC$RA <- as.numeric(species.DC$RA)

species.DC <- species.DC[,-c(2:8, 10)] # keep only the variables you mentioned 462
species.DC <- filter(species.DC, RA != 0) #removing 0's 299

##now running the loop
species.DC$Species <- as.factor(species.DC$Species)
regDS <-species.DC %>% 
  group_by(Species) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectDS <- regDS
selectDS[which(selectDS$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectDS[which(selectDS$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectDS <- selectDS[order(selectDS$term), ]

select_betaregDS = selectDS[-c(1:40,61:101),]
select_betaregDS <- select_betaregDS[order(select_betaregDS$estimate), ]
select_betaregDS <- select_betaregDS[-2]
select_betaregSigDS <- filter(select_betaregDS, sig=="P-Value < 0.05")

#merge taxonomy with the significant Gens
SigSDraxxin <- merge(select_betaregSigDS, ASVkey2, by.x = "Species", by.y = "ASVnos") 
#only lactobacillus was significant

#Nuflur
N_OTUGP <- as.data.frame(taxaG[rownames(taxaG) %in% rownames(Nuflur),]) #915 genera
N_OTUG <- t(N_OTUGP)
table <- as.data.frame(N_OTUG)
sum(table$ds1_15) #1889 sequences
N_OTUG <- prop.table((N_OTUG), 2) #to calculate relative abundance for each sample but not grouped
colSums(N_OTUG)
N_OTUG <- N_OTUG[apply(N_OTUG, 1, function(x) mean(x) > 0.001),] #filter the ASVs with an abundance < than 0.001, from 96 to 10
colSums(N_OTUG)
names <- Nuflur$sample
taxaNG2 <- cbind(names, N_OTUG)
taxaNG2 <- taxaNG2[,-c(1)]
taxaNG3<- reshape2::melt(taxaNG2[, c(1:10)])
str(taxaNG3)
taxaNG3$value <- as.numeric(taxaNG3$value)
sum(taxaNG3$value)
colnames(taxaNG3) <- c("Genus", "sample", "RA")
sample <- subset(taxaNG3, sample == "ds1_15")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 0.9562212

##now mergin the summarized phylum table with the metadata
genus.NC <- Nuflur %>% 
  inner_join(taxaNG3, by = "sample")
summary(genus.NC)
genus.NC$Genus <- as.character(genus.NC$Genus)
genus.NC$RA <- as.numeric(genus.NC$RA)

genus.NC <- genus.NC[,-c(2:8, 10)] # keep only the variables you mentioned 960
genus.NC <- filter(genus.NC, RA != 0) #removing 0's 386

##now running the loop
str(genus.NC)
genus.NC$Genus <- as.factor(genus.NC$Genus)
regNG <-genus.NC %>% 
  group_by(Genus) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectNG <- regNG
selectNG[which(selectNG$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectNG[which(selectNG$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectNG <- selectNG[order(selectNG$term), ]

select_betaregNG = selectNG[-c(1:132,199:360),]
select_betaregNG <- select_betaregNG[order(select_betaregNG$estimate), ]
select_betaregNG <- select_betaregNG[-2]
select_betaregSigNG <- filter(select_betaregNG, sig=="P-Value < 0.05")
#merge taxonomy with the significant Gens
SigGNuflur <- merge(select_betaregSigNG, ASVkey, by.x = "Genus", by.y = "ASVnos") 
write.csv(SigGNuflur, "SigGNuflur.csv")
SigGNuflur <- read.csv("SigGNuflur.csv")
SigGNuflur <- SigGNuflur[order(SigGNuflur$estimate), ]

my_colors2 <- c(
    '#b3df8a',"#5E738F","coral","lightpink",
    '#a6cee3','#a6cee3','#cab3d6','#6a3d9a','#ffff99','#b15938', 
    "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
    "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
    "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
  )

y <- ggplot(SigGNuflur, aes(x=reorder(genus, estimate), y=estimate, color=col1)) +
  geom_segment(aes(x = reorder(genus, estimate), xend =genus, y = 0, yend = estimate)) +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values = my_colors2) +
  labs(title = "Genus") +
  labs(x='', y='Regression Coefficient')+
  ylim(-1.3, 0.7) +
  geom_hline(yintercept=0) +
  theme_bw() +
  labs(color= "Phylum") +
  coord_flip() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.title.x = element_text(color="black", size=10, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 9), axis.text.y = element_text(color = "black", size = 10, face ="italic")) 
y

ggarrange(X,Y,x,y,xx,yy, labels = c("a", "b", "c","d", "e","f"),
          nrow=3, ncol=2)
ggarrange(X,Y,x,y, labels = c("a", "b", "c","d"),
          nrow=2, ncol=2)

X
x
xx
Y
y
yy

##Species Nuflur
N_OTUSP <- as.data.frame(taxaS[rownames(taxaS) %in% rownames(Nuflur),]) #1507 genera
N_OTUS <- t(N_OTUSP)
table <- as.data.frame(N_OTUS)
sum(table$ds0_19) #1847 sequences
N_OTUS <- prop.table((N_OTUS), 2) #to calculate relative abundance for each sample but not grouped
colSums(N_OTUS)
N_OTUSpecies = N_OTUS
N_OTUS <- N_OTUS[apply(N_OTUS, 1, function(x) mean(x) > 0.01),] #filter the ASVs with an abundance < than 0.001, from 1507 to 21
colSums(N_OTUS)
names <- Nuflur$sample
taxaNS2 <- cbind(names, N_OTUS)
taxaNS2 <- taxaNS2[,-c(1)]
taxaNS3<- reshape2::melt(taxaNS2[, c(1:10)])
str(taxaNS3)
taxaNS3$value <- as.numeric(taxaNS3$value)
sum(taxaNS3$value)
colnames(taxaNS3) <- c("Species", "sample", "RA")
sample <- subset(taxaNS3, sample == "ds0_19")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 0.8971305

##now mergin the summarized phylum table with the metadata
species.NC <- Nuflur %>% 
  inner_join(taxaNS3, by = "sample")
summary(species.NC)
species.NC$Species <- as.character(species.NC$Species)
species.NC$RA <- as.numeric(species.NC$RA)

####separate by treatment CON and SCFP
species.NC <- species.NC[,-c(2:8, 10)] # keep only the variables you mentioned 210
species.NC <- filter(species.NC, RA != 0) #removing 0's 151

##now running the loop
species.NC$Species <- as.factor(species.NC$Species)
regNS <-species.NC %>% 
  group_by(Species) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectNS <- regNS
selectNS[which(selectNS$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectNS[which(selectNS$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectNS <- selectNS[order(selectNS$term), ]

select_betaregNS = selectNS[-c(1:30,46:81),]
select_betaregNS <- select_betaregNS[order(select_betaregNS$estimate), ]
select_betaregNS <- select_betaregNS[-2]
select_betaregSigNS <- filter(select_betaregNS, sig=="P-Value < 0.05")

#merge taxonomy with the significant Gens
SigSNuflur <- merge(select_betaregSigNS, ASVkey2, by.x = "Species", by.y = "ASVnos") 

##now look for the BRD pathogens on the sick samples
#draxxin
taxBRD.D <- D_OTUSpecies[rownames(D_OTUSpecies) %in% rownames(patho),]
names <- Draxxin$sample
taxaS2.1 <- cbind(names, taxBRD.D)
taxaS2.1 <- taxaS2.1[,-c(1)]
taxaS3.1<- reshape2::melt(taxaS2.1[, c(1:22)])
str(taxaS3.1)
taxaS3.1$value <- as.numeric(taxaS3.1$value)
sum(taxaS3.1$value) #0.7960569
colnames(taxaS3.1) <- c("Specie", "sample", "RA")
sample <- subset(taxaS3.1, sample == "ds0_18")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts=0.005221932

##now mergin the summarized phylum table with the metadata
species.D1 <- Draxxin %>% 
  inner_join(taxaS3.1, by = "sample")
summary(species.D1)
species.D1$Specie <- as.character(species.D1$Specie)
species.D1$RA <- as.numeric(species.D1$RA)

####separate by treatment CON and SCFP
DraxxinCS1 <- species.D1[,-c(2:8,10,11)] # keep only the variables you mentioned 6375
DraxxinCS1 <- filter(DraxxinCS1, RA != 0) #removing 0's 15

### test analysis
test <- subset(DraxxinCS1, Specie=="ASV1374")
test$calf <- as.factor(test$calf)
test <- filter(test, RA != 0)
x <- glmmadmb(RA~d+(d|calf),
              data=test, family="beta")
summary(x)

##now running the loop
str(DraxxinCS1)
DraxxinCS1$Specie <- as.factor(DraxxinCS1$Specie)
regHS <-DraxxinCS1 %>% 
  group_by(Specie) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectS <- regHS
selectS[which(selectS$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectS[which(selectS$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectS <- selectS[order(selectS$term), ]

## nuflur
taxBRD.N <- N_OTUSpecies[rownames(N_OTUSpecies) %in% rownames(patho),]
names <- Nuflur$sample
taxaS2.1 <- cbind(names, taxBRD.N)
taxaS2.1 <- taxaS2.1[,-c(1)]
taxaS3.1<- reshape2::melt(taxaS2.1[, c(1:10)])
str(taxaS3.1)
taxaS3.1$value <- as.numeric(taxaS3.1$value)
sum(taxaS3.1$value) #0.08496721
colnames(taxaS3.1) <- c("Specie", "sample", "RA")
sample <- subset(taxaS3.1, sample == "ds0_19")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 0.002165674

##now mergin the summarized phylum table with the metadata
species.N1 <- Nuflur %>% 
  inner_join(taxaS3.1, by = "sample")
summary(species.N1)
species.N1$Specie <- as.character(species.N1$Specie)
species.N1$RA <- as.numeric(species.N1$RA)

####separate by treatment CON and SCFP
NuflurCS1 <- species.N1[,-c(2:8,10,11)] # keep only the variables you mentioned 6375
NuflurCS1 <- filter(NuflurCS1, RA != 0) #removing 0's 8

##now running the loop
str(NuflurCS1)
NuflurCS1$Specie <- as.factor(NuflurCS1$Specie)
regHS <-NuflurCS1 %>% 
  group_by(Specie) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectS <- regHS
selectS[which(selectS$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectS[which(selectS$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectS <- selectS[order(selectS$term), ]


## -----------------------------------DESeq for the sickness effect-----------
## Comparing healthy and BRD animals
## We will use the day h7 and day s0
str(metadata)
comp7_0 <- subset(metadata, subset=day %in% c("h7", "s0"))
comp7_0 <- subset(comp7_0, subset=Antibiotic %in% c("No"))
str(comp7_0)
comp7_0$calf = as.factor(comp7_0$calf)
comp7_0 <- comp7_0[-c(1,7,12,20),]

NewASVtable2 = t(NewASVtable)
S_OTU <- as.data.frame(NewASVtable2[rownames(NewASVtable2) %in% rownames(comp7_0),]) #51, 37 phylum

#DeSeq

S_OTU  <- S_OTU  + 1

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(S_OTU), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
comp7_0$Sickness <- as.factor(comp7_0$Sickness)
meta.physeq = sample_data(comp7_0)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq #[ 3837 taxa and 19 samples ]

levels(comp7_0$Sickness)

# establishing the model
diagdds = phyloseq_to_deseq2(physeq_deseq, ~Sickness)
#PenCode needs to be factor
diagdds = DESeq(diagdds, test="Wald",fitType = "local")
head(diagdds)
resultsNames(diagdds)

par <- estimateDispersions(diagdds, fitType = "parametric")
loc <- estimateDispersions(diagdds, fitType = "local")
mea <- estimateDispersions(diagdds, fitType = "mean")
plotDispEsts(par, main= "dispEst: parametric")
plotDispEsts(loc, main= "dispEst: local")
plotDispEsts(mea, main= "dispEst")

my_contrast = c("Sickness", "Got", "Never")
res = results(diagdds, contrast = my_contrast, cooksCutoff = FALSE, alpha=0.05, pAdjustMethod = "BH")
summary(res)
res
res <- as.data.frame(res)
res2 <- merge(res, NewTax, by.x = 0, by.y = 0)
alpha = 0.05
#sigtab = res ### No significant results
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))
sigtab

sigtab$High_low <- ifelse(
  sigtab$log2FoldChange < -2.00, 'High in Healthy',
  ifelse(sigtab$log2FoldChange > 2.00, 'High in Sick',
         'Mid Change'))
write.table(sigtab,"sigtabSicknessCalfNasal.txt",sep=",", row.names = TRUE)

###--------------let's compare day 0 to day 1 to identify any antibiotic effect
str(Draxxin)
Draxxin$day = as.factor(Draxxin$day)
levels(Draxxin$day)
d0_1 = Draxxin[Draxxin$day %in% c('s0', 's1'), ]

Nuflur$day = as.factor(Nuflur$day)
levels(Nuflur$day)
n0_1 = Nuflur[Nuflur$day %in% c('s0', 's1'), ]

###Draxxin
D_01 <- as.data.frame(taxaG[rownames(taxaG) %in% rownames(d0_1),]) #916 genera
D_01G <- t(D_01)
table <- as.data.frame(D_01G)
sum(table$ds0_9) #1739 sequences
D_01G <- prop.table((D_01G), 2) #to calculate relative abundance for each sample but not grouped
colSums(D_01G)
D_01G <- D_01G[apply(D_01G, 1, function(x) mean(x) > 0.001),] #filter the ASVs with an abundance < than 0.001, from 916 to 100
colSums(D_01G)
names <- d0_1$sample
taxDG2 <- cbind(names, D_01G)
taxDG2 <- taxDG2[,-c(1)]
taxDG3<- reshape2::melt(taxDG2[, c(1:11)])
str(taxDG3)
taxDG3$value <- as.numeric(taxDG3$value)
sum(taxDG3$value)
colnames(taxDG3) <- c("Genus", "sample", "RA")
sample <- subset(taxDG3, sample == "ds0_9")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 0.7096032

##now mergin the summarized phylum table with the metadata
genus.D01 <- d0_1 %>% 
  inner_join(taxDG3, by = "sample")
summary(genus.D01)
genus.D01$Genus <- as.character(genus.D01$Genus)
genus.D01$RA <- as.numeric(genus.D01$RA)

genus.D01 <- genus.D01[,-c(2:8, 10)] # keep only the variables you mentioned 1100
genus.D01 <- filter(genus.D01, RA != 0) #removing 0's 420

##now running the loop
genus.D01$d <- as.factor(genus.D01$d)
genus.D01$Genus <- as.factor(genus.D01$Genus)
regD01 <-genus.D01 %>% 
  group_by(Genus) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectD01 <- regD01
selectD01[which(selectD01$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectD01[which(selectD01$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectD01 <- selectD01[order(selectD01$term), ]

select_betaregD01 = selectD01[-c(1:106,160:312),]
select_betaregD01 <- select_betaregD01[order(select_betaregD01$estimate), ]
select_betaregD01 <- select_betaregD01[-2]
select_betaregSigD01 <- filter(select_betaregD01, sig=="P-Value < 0.05")
SigDG01 <- merge(select_betaregSigD01, ASVkey, by.x = "Genus", by.y = "ASVnos")
write.csv(SigDG01, "SigDraxxinG01.csv")

##Species Draxxin
D_01S <- as.data.frame(taxaS[rownames(taxaS) %in% rownames(d0_1),]) #1507 genera
D01S <- t(D_01S)
table <- as.data.frame(D01S)
sum(table$ds0_9) #1739 sequences
D01S <- prop.table((D01S), 2) #to calculate relative abundance for each sample but not grouped
colSums(D01S)
D01SSpecies = D01S
D01S <- D01S[apply(D01S, 1, function(x) mean(x) > 0.01),] #filter the ASVs with an abundance < than 0.001, from 1507 to 11
colSums(D01S)
names <- d0_1$sample
taxDS2 <- cbind(names, D01S)
taxDS2 <- taxDS2[,-c(1)]
taxDS3<- reshape2::melt(taxDS2[, c(1:11)])
str(taxDS3)
taxDS3$value <- as.numeric(taxDS3$value)
sum(taxDS3$value)
colnames(taxDS3) <- c("Species", "sample", "RA")
sample <- subset(taxDS3, sample == "ds0_9")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 0.2455434


##now mergin the summarized phylum table with the metadata
species.D01 <-d0_1 %>% 
  inner_join(taxDS3, by = "sample")
summary(species.D01)
species.D01$Species <- as.character(species.D01$Species)
species.D01$RA <- as.numeric(species.D01$RA)

species.D01 <- species.D01[,-c(2:8, 10)] # keep only the variables you mentioned 253
species.D01 <- filter(species.D01, RA != 0) #removing 0's 135

##now running the loop
species.D01$d <- as.factor(species.D01$d)
species.D01$Species <- as.factor(species.D01$Species)
regDS <-species.D01 %>% 
  group_by(Species) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectDS <- regDS
selectDS[which(selectDS$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectDS[which(selectDS$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectDS <- selectDS[order(selectDS$term), ]

select_betaregDS = selectDS[-c(1:28,43:79),]
select_betaregDS <- select_betaregDS[order(select_betaregDS$estimate), ]
select_betaregDS <- select_betaregDS[-2]

select_betaregSigDS <- filter(select_betaregDS, sig=="P-Value < 0.05")
SigDSp <- merge(select_betaregSigDS, ASVkey2, by.x = "Species", by.y = "ASVnos") 
write.csv(SigDSp, "SigDraxxinSp01.csv")

AntDraxxin01 <- read.csv("AntDraxxin01.csv")
DrGenus = subset(AntDraxxin01, test == "Genus")

a <- ggplot(DrGenus, aes(x=reorder(sigASV, estimate), y=estimate)) +
  geom_segment(aes(x = reorder(sigASV, estimate), xend =sigASV, y = 0, yend = estimate)) +
  geom_point(aes(shape = sig), size = 3) +
  facet_free(test~.) +
  theme(strip.text.y = element_text(size = 15, color = "black", face = "bold.italic")) +
  scale_shape_manual(values = c(16, 1)) +
  labs(x='', y='Regression Coefficient')+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Tulathromycin: Comparison d0 vs d1") +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 10, face ="italic")) 

#----------Nuflur samples
N_01 <- as.data.frame(taxaG[rownames(taxaG) %in% rownames(n0_1),]) #916 genera
N_01G <- t(N_01)
table <- as.data.frame(N_01G)
sum(table$ds0_19) #1847 sequences
N_01G <- prop.table((N_01G), 2) #to calculate relative abundance for each sample but not grouped
colSums(N_01G)
N_01G <- N_01G[apply(N_01G, 1, function(x) mean(x) > 0.001),] #filter the ASVs with an abundance < than 0.001, from 915 to 81
colSums(N_01G)
names <- n0_1$sample
taxNG2 <- cbind(names, N_01G)
taxNG2 <- taxNG2[,-c(1)]
taxNG3<- reshape2::melt(taxNG2[, c(1:5)])
str(taxNG3)
taxNG3$value <- as.numeric(taxNG3$value)
sum(taxNG3$value)
colnames(taxNG3) <- c("Genus", "sample", "RA")
sample <- subset(taxNG3, sample == "ds0_19")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts= 0.983216


##now mergin the summarized phylum table with the metadata
genus.N01 <- n0_1 %>% 
  inner_join(taxNG3, by = "sample")
summary(genus.N01)
genus.N01$Genus <- as.character(genus.N01$Genus)
genus.N01$RA <- as.numeric(genus.N01$RA)

genus.N01 <- genus.N01[,-c(2:8, 10)] # keep only the variables you mentioned 405
genus.N01 <- filter(genus.N01, RA != 0) #removing 0's 194

##now running the loop
genus.N01$d <- as.factor(genus.N01$d)
genus.N01$Genus <- as.factor(genus.N01$Genus)
regN01 <-genus.N01 %>% 
  group_by(Genus) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectN01 <- regN01
selectN01[which(selectN01$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectN01[which(selectN01$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectN01 <- selectN01[order(selectN01$term), ]

select_betaregN01 = selectN01[-c(1:58,88:197),]
select_betaregN01 <- select_betaregN01[order(select_betaregN01$estimate), ]
select_betaregN01 <- select_betaregN01[-2]

select_betaregSigN01 <- filter(select_betaregN01, sig=="P-Value < 0.05")
SigNG01 <- merge(select_betaregSigN01, ASVkey, by.x = "Genus", by.y = "ASVnos")
write.csv(SigNG01, "SigNuflurG01.csv")

##Species Draxxin
N_01S <- as.data.frame(taxaS[rownames(taxaS) %in% rownames(n0_1),]) #1507 genera
N01S <- t(N_01S)
table <- as.data.frame(N01S)
sum(table$ds0_19) #1847 sequences
N01S <- prop.table((N01S), 2) #to calculate relative abundance for each sample but not grouped
colSums(N01S)
N01SSpecies = N01S
N01S <- N01S[apply(N01S, 1, function(x) mean(x) > 0.01),] #filter the ASVs with an abundance < than 0.001, from 1507 to 21
colSums(N01S)
names <- n0_1$sample
taxNS2 <- cbind(names, N01S)
taxNS2 <- taxNS2[,-c(1)]
taxNS3<- reshape2::melt(taxNS2[, c(1:5)])
str(taxNS3)
taxNS3$value <- as.numeric(taxNS3$value)
sum(taxNS3$value)
colnames(taxNS3) <- c("Species", "sample", "RA")
sample <- subset(taxNS3, sample == "ds0_19")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ##

##now mergin the summarized phylum table with the metadata
species.N01 <-n0_1 %>% 
  inner_join(taxNS3, by = "sample")
summary(species.N01)
species.N01$Species <- as.character(species.N01$Species)
species.N01$RA <- as.numeric(species.N01$RA)

species.N01 <- species.N01[,-c(2:8, 10)] # keep only the variables you mentioned 105
species.N01 <- filter(species.N01, RA != 0) #removing 0's 68

##now running the loop
species.N01$d <- as.factor(species.N01$d)
species.N01$Species <- as.factor(species.N01$Species)
regNS <-species.N01 %>% 
  group_by(Species) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectNS <- regNS
selectNS[which(selectNS$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectNS[which(selectNS$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectNS <- selectNS[order(selectNS$term), ]

select_betaregNS = selectNS[-c(1:12,19:45),]
select_betaregNS <- select_betaregNS[order(select_betaregNS$estimate), ]
select_betaregNS <- select_betaregNS[-2]

select_betaregSigNS <- filter(select_betaregNS, sig=="P-Value < 0.05")
SigNSp <- merge(select_betaregSigNS, ASVkey2, by.x = "Species", by.y = "ASVnos") 
write.csv(SigNSp, "SigNuflurSp01.csv")

AntNuflur01 <- read.csv("AntNuflur01.csv")
NuGenus = subset(AntNuflur01, test =="Genus")

b = ggplot(NuGenus, aes(x=reorder(sigASV, estimate), y=estimate)) +
  geom_segment(aes(x = reorder(sigASV, estimate), xend =sigASV, y = 0, yend = estimate)) +
  geom_point(aes(shape = sig), size = 3) +
  facet_free(test~.) +
  theme(strip.text.y = element_text(size = 15, color = "black", face = "bold.italic")) +
  scale_shape_manual(values = c(16, 1)) +
  labs(x='', y='Regression Coefficient')+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Florfenicol: Comparison d0 vs d1") +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 10, face ="italic")) 

ggarrange(a,b, labels = c("a", "b"),
          nrow=2)

########### Now we compare the entire community to make sure how the antibiotic affect before treat to after treatment
d0_5 = Draxxin[Draxxin$day %in% c('s0', 's5'), ]
d0_10 = Draxxin[Draxxin$day %in% c('s0', 's10'), ]

n0_5 = Nuflur[Nuflur$day %in% c('s0', 's5'), ]
n0_10 = Nuflur[Nuflur$day %in% c('s0', 's10'), ]

## Let's start with the draxxin samples, you change in the code the metadata
###creating the asv_table 
#D_day <- as.data.frame(taxa[rownames(taxa) %in% rownames(d0_10),]) #for phylum
D_day <- as.data.frame(taxaG[rownames(taxaG) %in% rownames(d0_10),]) #for genus
D_dayG <- t(D_day)
table <- as.data.frame(D_dayG)
sum(table$ds0_9) #1739 sequences
D_dayG <- prop.table((D_dayG), 2) #not subset
names <- d0_1$sample
#names <- d0_5$sample
#names <- d0_10$sample
taxaD2 <- cbind(names, D_dayG)
taxaD2 <- taxaD2[,-c(1)]
taxaD3<- reshape2::melt(taxaD2[, c(1:11)]) # change this value based on the # of samples your metadata has
str(taxaD3)
taxaD3$value <- as.numeric(taxaD3$value)
sum(taxaD3$value)
#colnames(taxaD3) <- c("Phylum", "sample", "RA")
colnames(taxaD3) <- c("Genus", "sample", "RA")
sample <- subset(taxaD3, sample == "ds0_9")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts=1

##now mergin the summarized phylum table with the metadata
phylum.D <- d0_10 %>% #change metadata here
  inner_join(taxaD3, by = "sample")
summary(phylum.D)
phylum.D$Phylum <- as.character(phylum.D$Phylum)
phylum.D$RA <- as.numeric(phylum.D$RA)

phylum.D <- phylum.D[,-c(2:8, 10)] # keep only the variables you mentioned 307 (d0_1 and d01_5), 407 (d0_10)
phylum.D <- filter(phylum.D, RA != 0) #removing 0's 117 (d0_1) 106 (d0_5), 126 (d0_10)

##now running the loop
phylum.D$d <- as.factor(phylum.D$d)
phylum.D$Phylum <- as.factor(phylum.D$Phylum)
regD <-phylum.D %>% 
  group_by(Phylum) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectD <- regD
selectD[which(selectD$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectD[which(selectD$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectD <- selectD[order(selectD$term), ]

select_betaregD = selectD[-c(1:22,34:69),]
select_betaregD <- select_betaregD[order(select_betaregD$estimate), ]
select_betaregD <- select_betaregD[-2]
select_betaregSigD <- filter(select_betaregD, sig=="P-Value < 0.05") 
#difference d0_1 Euryarchaeota -2.34805
#no difference d0_5
# d0_10 only chloroflexi

#### Genus
##now mergin the summarized phylum table with the metadata
genus.d <- d0_10 %>% ##change the metadata gere
  inner_join(taxaD3, by = "sample")
summary(genus.d)
genus.d$Genus <- as.character(genus.d$Genus)
genus.d$RA <- as.numeric(genus.d$RA)
genus.d <- genus.d[,-c(2:8, 10)] # keep only the variables you mentioned 10076 (d0_1 and d0_10), 9160 (d0_5)
genus.d <- filter(genus.d, RA != 0) #removing 0's 717(d0_1) 676 (d0_5), 834 (d0_10) 

##now running the loop
genus.d$d <- as.factor(genus.d$d)
genus.d$Genus <- as.factor(genus.d$Genus)
regD <-genus.d %>% 
  group_by(Genus) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectD <- regD
selectD[which(selectD$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectD[which(selectD$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectD <- selectD[order(selectD$term), ]
select_betaregD = selectD[-c(1:238,358:810),]
select_betaregD <- select_betaregD[order(select_betaregD$estimate), ]
select_betaregD <- select_betaregD[-2]
select_betaregSigD <- filter(select_betaregD, sig=="P-Value < 0.05") #no difference
SigDG010 <- merge(select_betaregSigD, ASVkey, by.x = "Genus", by.y = "ASVnos") 
write.csv(SigDG010, "SigDraxxinG010.csv")

## now plot it
AntDraxxin <- read.csv("BetaRegressionDraxxin0.5.10.csv")
AntDraxxin$tax = as.factor(AntDraxxin$tax)
levels(AntDraxxin$tax) <- list("d0 vs d1"= "d0 vs d1","d0 vs d5"= "d0 vs d5", "d0 vs d10"= "d0 vs d10")

a = ggplot(AntDraxxin, aes(x=reorder(genus, estimate), y=estimate)) +
  geom_segment(aes(x = reorder(genus, estimate), xend =genus, y = 0, yend = estimate)) +
  geom_point(aes(shape = sig), size = 3) +
  facet_free(tax~.) +
  theme(strip.text.y = element_text(size = 15, color = "black", face = "bold.italic")) +
  scale_shape_manual(values = c(16, 1)) +
  labs(x='', y='Regression Coefficient')+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Tulathromycin") +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 10, face ="italic")) 
a
#### Now for nuflur samples
#N_day <- as.data.frame(taxa[rownames(taxa) %in% rownames(n0_10),]) #for phylum
N_day <- as.data.frame(taxaG[rownames(taxaG) %in% rownames(n0_10),]) #for genus
N_dayG <- t(N_day)
table <- as.data.frame(N_dayG)
sum(table$ds0_19) #1847 sequences
N_dayG <- prop.table((N_dayG), 2) #not subset
names <- n0_10$sample #change metadata
taxaN2 <- cbind(names, N_dayG)
taxaN2 <- taxaN2[,-c(1)]
taxaN3<- reshape2::melt(taxaN2[, c(1:5)]) #missing ds5_30
str(taxaN3)
taxaN3$value <- as.numeric(taxaN3$value)
sum(taxaN3$value)
#colnames(taxaN3) <- c("Phylum", "sample", "RA")
colnames(taxaN3) <- c("Genus", "sample", "RA")
sample <- subset(taxaN3, sample == "ds0_19")
sample$RA <- as.numeric(sample$RA)
sum(sample$RA) ## RA=1 counts=1

##now mergin the summarized phylum table with the metadata
phylum.N <- n0_10 %>%  #change metadata
  inner_join(taxaN3, by = "sample")
summary(phylum.N)
phylum.N$Phylum <- as.character(phylum.N$Phylum)
phylum.N$RA <- as.numeric(phylum.N$RA)

phylum.N <- phylum.N[,-c(2:8, 10)] # keep only the variables you mentioned 185(d0_1), 148 (d0_5). 185 (d0_10)
phylum.N <- filter(phylum.N, RA != 0) #removing 0's 44 (d0_1) 33 (d0_5),45 (d0_10)

##now running the loop
phylum.N$d <- as.factor(phylum.N$d)
phylum.N$Phylum <- as.factor(phylum.N$Phylum)
regN <-phylum.N %>% 
  group_by(Phylum) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectN <- regN
selectN[which(selectN$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectN[which(selectN$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectN <- selectN[order(selectN$term), ]

select_betaregN = selectN[-c(1:12,19:42),]
select_betaregN <- select_betaregN[order(select_betaregN$estimate), ]
select_betaregN <- select_betaregN[-2]
select_betaregSigN <- filter(select_betaregN, sig=="P-Value < 0.05") 
# no difference d0_1
#no difference d0_5
# d0_10 only chloroflexi

#### Genus
##now mergin the summarized phylum table with the metadata
genus.n <- n0_10 %>% 
  inner_join(taxaN3, by = "sample")
summary(genus.n)
genus.n$Genus <- as.character(genus.n$Genus)
genus.n$RA <- as.numeric(genus.n$RA)
genus.n <- genus.n[,-c(2:8, 10)] # keep only the variables you mentioned 4580 (d0_1), 3664 (d0_5), 4580
genus.n <- filter(genus.n, RA != 0) #removing 0's 276 (d0_1) 226 (d0_5), 309 (d0_10) 

##now running the loop
genus.n$d <- as.factor(genus.n$d)
genus.n$Genus <- as.factor(genus.n$Genus)
regN <-genus.n %>% 
  group_by(Genus) %>% 
  nest() %>% 
  mutate(model_results = map(data, ~tidy(safe_betareg(RA ~ d + (d|calf),family="beta", .)))) %>% 
  unnest(model_results)  %>% 
  mutate(p.adj = p.adjust(p.value, "BH"))###give the regression value for every phylum

selectN <- regN
selectN[which(selectN$'p.adj' < 0.05),'sig'] <- 'P-Value < 0.05'
selectN[which(selectN$'p.adj' >= 0.05),'sig'] <- 'P-Value > 0.05'
selectN <- selectN[order(selectN$term), ]
select_betaregN = selectN[-c(1:67,100:274),]
select_betaregN <- select_betaregN[order(select_betaregN$estimate), ]
select_betaregN <- select_betaregN[-2]
select_betaregSigN <- filter(select_betaregN, sig=="P-Value < 0.05") #no difference
SigNG010 <- merge(select_betaregSigN, ASVkey, by.x = "Genus", by.y = "ASVnos") 
write.csv(SigNG010, "SigNuflurG010.csv")

## now plot it
AntNuflur <- read.csv("BetaRegressionNuf0.5.10.csv")
AntNuflur$tax = as.factor(AntNuflur$tax)
levels(AntNuflur$tax) <- list("d0 vs d1"= "d0 vs d1","d0 vs d5"= "d0 vs d5", "d0 vs d10"= "d0 vs d10")

b = ggplot(AntNuflur, aes(x=reorder(genus, estimate), y=estimate)) +
  geom_segment(aes(x = reorder(genus, estimate), xend =genus, y = 0, yend = estimate)) +
  geom_point(aes(shape = sig), size = 3) +
  facet_free(tax~.) +
  theme(strip.text.y = element_text(size = 15, color = "black", face = "bold.italic")) +
  scale_shape_manual(values = c(16, 1)) +
  labs(x='', y='Regression Coefficient')+
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Florfenicol") +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 10, face ="italic")) 
b
ggarrange(a,b, labels = c("a", "b"),
          ncol=2)

### Now check the abundance of BRD pathobionts on day 0 compared to the other days
#obtain the relative abundance of the data calculated with all the other samples
#change the metadata 
D_010S <- as.data.frame(taxaS[rownames(taxaS) %in% rownames(n0_10),]) #for genus
D_010S <- t(D_010S)
table <- as.data.frame(D_010S)
D_010S <- prop.table((D_010S), 2) #not subset
#names <- d0_1$sample
#names <- d0_5$sample
#names <- d0_10$sample
#names <- n0_1$sample
#names <- n0_5$sample
names <- n0_10$sample
taxaD010 <- cbind(names, D_010S)
taxaD010 <- taxaD010[,-c(1)]
taxaD010 <- reshape2::melt(taxaD010[, c(1:5)]) # change this value based on the # of samples your metadata has
str(taxaD010)
taxaD010$value <- as.numeric(taxaD010$value)
sum(taxaD010$value)
colnames(taxaD010) <- c("Species", "sample", "RA")

##relative abundance for each day
taxaD01
taxaD05
taxaD0110

#now select the ASVs that are classified as BRD pathobionts
pathos = rownames(patho)
taxaD01 = taxaD01[taxaD01$Species %in% pathos, ]
taxaD05 = taxaD05[taxaD05$Species %in% pathos, ]
taxaD010 = taxaD010[taxaD010$Species %in% pathos, ]

#now calculate the beta regression
pathoTax <- n0_10 %>%  #change metadata
  inner_join(taxaD010, by = "sample")
summary(pathoTax)
pathoTax$Species <- as.character(pathoTax$Species)
pathoTax$RA <- as.numeric(pathoTax$RA)
pathoTax <- pathoTax[,-c(2:8, 10)] # samples day0_1 (44), d0_5 (11), d0_10 (11)
pathoTax$RelAbund = pathoTax$RA * 100
pathoTax$comparison = c("d0 vs d10")
FlorComp =pathoTax
FlorComp = rbind(FlorComp, pathoTax)
write.csv(FlorComp, "FlorComp.csv")

##------Let's check the abundance of the BRD-pathobionts on the healthy samples
str(healthy)
healhty2 = healthy
healhty2$d = as.factor(healhty2$d)

H0 = subset(healhty2, d == 0)
H7 = subset(healhty2, d == 7)
H14 = subset(healhty2, d == 14)

H_14 <- as.data.frame(taxaS[rownames(taxaS) %in% rownames(H14),]) #for genus
H_14S <- t(H_14)
table <- as.data.frame(H_14S)
H_14S <- prop.table((H_14S), 2) #not subset
names <- H14$sample
taxa14H <- cbind(names, H_14S)
taxa14H <- taxa14H[,-c(1)]
taxa14H <- reshape2::melt(taxa14H[, c(1:19)]) # change this value based on the # of samples your metadata has
str(taxa14H)
taxa14H$value <- as.numeric(taxa14H$value)
sum(taxa14H$value)
colnames(taxa14H) <- c("Species", "sample", "RA")

##relative abundance for each day
taxaH0
taxaH7
taxa14H

#now select the ASVs that are classified as BRD pathobionts
pathos = rownames(patho)
taxaH0 = taxaH0[taxaH0$Species %in% pathos, ]
taxaH7 = taxaH7[taxaH7$Species %in% pathos, ]
taxa14H = taxa14H[taxa14H$Species %in% pathos, ]

#now calculate the beta regression
pathoTax <- H14 %>%  #change metadata
  inner_join(taxa14H, by = "sample")
summary(pathoTax)
pathoTax$Species <- as.character(pathoTax$Species)
pathoTax$RA <- as.numeric(pathoTax$RA) # samples day0 (52), d0_5 (11), d0_10 (11)
pathoTax$RelAbund = pathoTax$RA * 100
pathoTax$comparison = c("d14")
#FlorComp =pathoTax
FlorComp = rbind(FlorComp, pathoTax)
write.csv(FlorComp, "HealthyComp.csv")

## functions
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

