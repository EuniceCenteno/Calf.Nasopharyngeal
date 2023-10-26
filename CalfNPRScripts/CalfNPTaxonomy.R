## This script goes through making a taxa bar plot and then running
## DESeq2 to find differentially abundant ASVs


# for help installing phyloseq, see this website
# https://bioconductor.org/packages/release/bioc/html/phyloseq.html

# to install phyloseq:
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("phyloseq")

library(qiime2R)
library(phyloseq)
library(zoo)
library(tidyverse)
library(tidyr) #for separate function
library(naniar)# Ffor replace all function
library("ape") #to create the tree
library(ggpubr)

##############################################
setwd("~/Desktop/eunice/PhD/CalfNasal/Qiime/CorrectFiltered/")
##
rm(ASV_table)
##### Taxonomy barplot
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
# Read the experimental design, and species classified documents
ASVs <- read_qza("filteredASV_table.qza")
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

# creating a phylosep object to change the taxa names
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
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

##now we separate the data of the healthy and sicka animal
#### Separate by healthy and sick
healthy <- subset(metadata, Sickness=="Never")
str(healthy)
healthy$d <- as.factor(healthy$d)
levels(healthy$d) <- list("0"="0", "7"="1", "14"="2")

H_OTU <- NewASVtable[rownames(NewASVtable) %in% rownames(healthy),] #51, 3837
write.csv(H_OTU, "HealthyASVTable.csv")
H_OTU2 = t(H_OTU)
H_NewTax <- NewTax[rownames(NewTax) %in% rownames(H_OTU2),]
write.csv(H_NewTax, "HealthyASVTaxonomy.csv")

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

#################################################################
##Taxa barplot
#################################################################
#Healthy animals
OTU.physeqH = otu_table(as.matrix(H_OTU), taxa_are_rows=FALSE)
tax.physeqH = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqH = sample_data(healthy)
physeq_H = phyloseq(OTU.physeqH, tax.physeqH, meta.physeqH)

##Making the tree, we need the phyloseq object 
random_tree = rtree(ntaxa(physeq_H), rooted=TRUE, tip.label=taxa_names(physeq_H))
##Merging the tree with the phyloseq object
physeqH1 = merge_phyloseq(physeq_H, meta.physeq, random_tree)
physeqH1 #This command should show the otu table, sample data, tax table and the phy tree
physeq_bar_plot = physeqH1

my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)

#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Genus")
my_column <- "d"  #this is the metadata column that we will use in the taxa barplot

rm(taxa.summary)

abund_filter <- 0.02  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$d = factor(physeq_meta_filtered$d, c("0", "7", "14"))
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = get(my_column), y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~LitterTreatment) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    theme_bw()+
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=8, face = "italic", vjust = 1.5)) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample")) 
  ggsave(paste0("", ml, "BarPlot_", my_column, ".png"), height = 5, width = 4)
}
write.csv(physeq_meta_filtered, "physeq_meta_filteredHGenus.csv")

### we need to know the abundance of the BRD associated pathogens probably on day 0
#only include the pathogens species
Hs = subset(NewTax, Genus == "Histophilus")
Pm = subset(NewTax, Genus == "Pasteurella")
Mh = subset(NewTax, Genus == "Mannheimia")
Mb = subset(NewTax, Genus == "Mycoplasma")
patho = rbind(Hs, Pm, Mh, Mb)

#calculating the relative abundance of the pathogens in the healthy group
otu.summary <- prop.table((H_OTU), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
#a <- as.data.frame(otu_abund)
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

str(melt_otu)
meta_otu <- merge(healthy, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu)
meta_otu_tax <- merge(meta_otu, patho, by.x = "ASV", by.y = 0)
str(meta_otu_tax)

SpeciesH <- meta_otu_tax %>% 
  group_by(Row.names, d, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 

my_colors2 <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#92C5DE","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
) ### for dairy samples

a <- ggplot(Species, aes(x = d, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors2) +
  ylim(c(0,0.3)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Healthy samples") +
  ylab(paste0("Relative Abundance")) +  labs(x='Day')

## now we do the same thing for the draxxin nuflur samples
#draxxin animals
str(Draxxin)

OTU.physeqD = otu_table(as.matrix(D_OTU), taxa_are_rows=FALSE)
tax.physeqD = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqD = sample_data(Draxxin)
physeq_D = phyloseq(OTU.physeqD, tax.physeqD, meta.physeqD)

##Making the tree, we need the phyloseq object 
random_treeD = rtree(ntaxa(physeq_D), rooted=TRUE, tip.label=taxa_names(physeq_D))
##Merging the tree with the phyloseq object
physeqD1 = merge_phyloseq(physeq_D, meta.physeqD, random_treeD)
physeqD1 #This command should show the otu table, sample data, tax table and the phy tree
physeq_bar_plotD = physeqD1

# Set colors for plotting
my_colorsD <- c(
  "#D6604D", '#a6cee3',"#F4A582", '#1f78b4','#b2df8a','#33a02c',
  '#5E738F',"#D1A33D","violet",'turquoise1','#fdbf6f',"#009999",
  '#ff7f00','#cab2d6', "olivedrab1","lightslateblue","slategray2",'#6a3d9a','#ffff99','#b15928', 
  "#CBD588")

my_colorsD1 <- c(
  "lightpink",'#a6cee3','#1f78b4','#b2df8a',"plum1",'powderblue','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00',"#599861", "cadetblue1",'#cab2d6',"#D1A33D",'#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black") #for famuly
my_colorsD2 <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)#for phylum

#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Phylum")
my_column <- "d"  #this is the metadata column that we will use in the taxa barplot

rm(taxa.summary)

abund_filter <- 0.02  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plotD %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$d = factor(physeq_meta_filtered$d, c("0", "1", "5", "10"))
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = get(my_column), y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~LitterTreatment) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colorsD2) +
    theme_bw()+
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=8, face = "italic", vjust = 1.5)) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample")) 
  ggsave(paste0("", ml, "BarPlotDraxxin_", my_column, ".png"), height = 5, width = 4)
}
write.csv(physeq_meta_filtered, "physeq_meta_filteredPhylumDraxxin.csv")

###Now the bacteria species for the Draxxin group
otu.summary <- prop.table((D_OTU), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
#a <- as.data.frame(otu_abund)
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

str(melt_otu)
meta_otu <- merge(Draxxin, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu)
meta_otu_taxD <- merge(meta_otu, patho, by.x = "ASV", by.y = 0)
str(meta_otu_taxD)

SpeciesD <- meta_otu_taxD %>% 
  group_by(Row.names, d, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 

my_colors2 <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#92C5DE","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
) ### for dairy samples

b <- ggplot(Species, aes(x = d, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors2) +
  ylim(c(0,0.3)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Tulathromycin samples") +
  ylab(paste0("Relative Abundance")) +  labs(x='Day')

#now nuflur
str(Nuflur)

OTU.physeqN = otu_table(as.matrix(N_OTU), taxa_are_rows=FALSE)
tax.physeqN = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqN = sample_data(Nuflur)
physeq_N = phyloseq(OTU.physeqN, tax.physeqN, meta.physeqN)

##Making the tree, we need the phyloseq object 
random_treeN = rtree(ntaxa(physeq_N), rooted=TRUE, tip.label=taxa_names(physeq_N))
##Merging the tree with the phyloseq object
physeqN1 = merge_phyloseq(physeq_N, meta.physeqN, random_treeN)
physeqN1 #This command should show the otu table, sample data, tax table and the phy tree
physeq_bar_plotN = physeqN1

# Set colors for plotting
# Set colors for plotting
my_colorsN <- c(
  "#D6604D", '#a6cee3', '#1f78b4','#b2df8a','#33a02c',"lightcoral",
  '#5E738F',"#D1A33D","#599861", "gray","#e31a1c", "#fdbf6f",
  '#ff7f00','#cab2d6', "lightseagreen","rosybrown2","lightcyan1","tomato","slategray2",'#6a3d9a',"orange", 
  "#b15928","#CBD588")

my_colorsN1 <- c(
  '#a6cee3','#1f78b4','#b2df8a',"plum1",'powderblue',"gray",'#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00',"cadetblue1",'#cab2d6',"#D1A33D","palegreen3",'#6a3d9a',"#b15928", 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black") #for famuly
my_colorsN2 <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)#for phylum

#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Phylum")
my_column <- "d"  #this is the metadata column that we will use in the taxa barplot

rm(taxa.summary)

abund_filter <- 0.02  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plotN %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$d = factor(physeq_meta_filtered$d, c("0", "1", "5", "10"))
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = get(my_column), y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~LitterTreatment) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colorsN2) +
    theme_bw()+
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=8, face = "italic", vjust = 1.5)) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample")) 
  ggsave(paste0("", ml, "BarPlotNuflur_", my_column, ".png"), height = 5, width = 4)
}
#write.csv(physeq_meta_filtered, "physeq_meta_filteredPhylumNuflur.csv")

###Now the bacteria species for the Draxxin group
otu.summary <- prop.table((N_OTU), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
#a <- as.data.frame(otu_abund)
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

str(melt_otu)
meta_otu <- merge(Nuflur, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu)
meta_otu_taxN <- merge(meta_otu, patho, by.x = "ASV", by.y = 0)
str(meta_otu_taxN)

SpeciesN <- meta_otu_taxN %>% 
  group_by(Row.names, d, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(d, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 

my_colors2 <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#92C5DE","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
) ### for dairy samples

c <- ggplot(Species, aes(x = d, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors2) +
  ylim(c(0,0.3)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Florfenicol samples") +
  ylab(paste0("Relative Abundance")) +  labs(x='Day')
write.csv(Species, "SpeciesNuflur.csv")

ggarrange(a,b,c, labels = c("a", "b", "c"),
          nrow=1, ncol = 3)
a
b
c

###we need to determine if the sequences of ASV assigned to Mycoplasma are the same or have the same sequence
# first we need to filter the ASV assigned to Mycoplasma
str(tax.clean)
myco = as.data.frame(subset(tax.clean, Genus == "Mycoplasma"))

#now we import the rep.seqs obtained from Qiime.
repseq <- read_qza("rep-seqs-filtered.qza")
rep_seq  <- as.data.frame(repseq$data)
rep_seq$ID = rownames(rep_seq)

#now we combine the mycoplasma data with the rep_seq
MycoSeq = merge(myco, rep_seq, by.x=0, by.y=0)
write.csv(MycoSeq, "MycoSeq.csv")

###we need to determine if the sequences of ASV assigned to Lactobacillus are the same or have the same sequence
# first we need to filter the ASV assigned to Mycoplasma
str(tax.clean)
Lacto = as.data.frame(subset(tax.clean, Genus == "Lactobacillus"))

#now we import the rep.seqs obtained from Qiime.
repseq <- read_qza("rep-seqs-filtered.qza")
rep_seq  <- as.data.frame(repseq$data)
rep_seq$ID = rownames(rep_seq)

#now we combine the mycoplasma data with the rep_seq
LactoSeq = merge(Lacto, rep_seq, by.x=0, by.y=0)
write.csv(LactoSeq, "LactoSeq2.csv")


##abundance of Mannheimia and Pasteurella on day s0 and ds1
day <- c("0", "1")
str(meta_otu_taxD)

MannD = subset(meta_otu_taxD, Genus =="Mannheimia")
MannD<-MannD[which(MannD$d==day),]
MannN = subset(meta_otu_taxN, Genus =="Mannheimia")
MannN<-MannN[which(MannN$d==day),]
Mann = rbind(MannD, MannN)

#remove 0
Mann0 = subset(Mann, Abundance > 0)

PauD = subset(meta_otu_taxD, Genus =="Pasteurella")
PauD<-PauD[which(PauD$d==day),]
PauN = subset(meta_otu_taxN, Genus =="Pasteurella")
PauN<-PauN[which(PauN$d==day),]
Pau = rbind(PauD, PauN)

#remove 0
Pau0 = subset(Pau, Abundance > 0)

#Create plot 
a = ggplot(Mann, aes(d, Abundance)) +
  geom_point(aes(colour = factor(calf), shape = factor(antibiotic)), size=3) +
  geom_line(aes(d, Abundance, group=factor(calf))) +
  theme_bw() +
  ggtitle("Mannheimia d0 and d1") # a lot of 0 

b = ggplot(Mann0, aes(d, Abundance)) +
  geom_point(aes(colour = factor(calf), shape = factor(antibiotic)), size=3) +
  geom_line(aes(d, Abundance, group=factor(calf))) +
  theme_bw() + # a lot of 0 
  ggtitle("Mannheimia - without relative abundance of 0") # a lot of 0 

c = ggplot(Pau, aes(d, Abundance, colour = factor(calf), shape = factor(antibiotic), size=3)) +
  geom_point(aes(colour = factor(calf), shape = factor(antibiotic)), size=3) +
  geom_line(aes(d, Abundance, group=factor(calf))) +
  theme_bw() +
  ggtitle("Pasteurella d0 and d1") # a lot of 0 

d = ggplot(Pau0, aes(d, Abundance)) +
  geom_point(aes(colour = factor(calf), shape = factor(antibiotic)), size=3) +
  geom_line(aes(d, Abundance, group=factor(calf))) +
  theme_bw() + # a lot of 0 
  ggtitle("Pasteurella - without relative abundance of 0") # a lot of 0

ggarrange(a,c,b,d, labels = c("a", "b", "c", ""),
          nrow=2, ncol = 2)

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


