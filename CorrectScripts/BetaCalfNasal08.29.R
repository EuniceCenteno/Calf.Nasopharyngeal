BETA diversity Swabs

library(vegan) 
library(ggplot2)
library(ggpubr)
library(data.table)
library(phyloseq)
library(qiime2R)
library(tidyr)
library(naniar)
library(raster)
library("ape")
library(dplyr)
#transpose function

setwd("~/Desktop/eunice/PhD/CalfNasal/Qiime/CorrectFiltered/") #sets new working directory for Windows systems (remember to replace … with your filepath)
metadata <-read.csv("CalfMetadata08.30.csv")
#OTU table (shared file)
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
str(metadata)
metadata$calf <- as.factor(metadata$calf)
metadata$Sickness <- as.factor(metadata$Sickness)
metadata$diagnostic <- factor(metadata$diagnostic)
levels(metadata$Sickness) <- list("Sick"="Got", "Healthy"="Never")
metadata$Antibiotic <- factor(metadata$Antibiotic)
rownames(metadata) <- metadata$ID

#OTU table, we use the rarified table
ASVs <- read_qza("rarefied_tableFiltered.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #3528 ASVs and 89
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
ASV_table <- t(ASV_table)
ASV_table <- merge(metadata, ASV_table, by.x = "ID", by.y = 0)
order_groups <- ASV_table$ID
row.names(ASV_table) = ASV_table[,1]
ASV_table = ASV_table[,-(1:8)]

#Taxonomy of each OTU
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


### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq
# otu_table()   OTU Table:         [ 3528 taxa and 89 samples ]
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
physeq_deseq
#otu_table()   OTU Table:         [ 3513 taxa and 89 samples ]
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
#write.table(NewTax,"NewTax.txt",sep=",", row.names = TRUE)
NewASVTable <- prunetable
NewASVTable <- NewASVTable[,-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]
NewASVTable = t(NewASVTable)

##Subsetting healthy and sick animals
healthy1 <- subset(metadata, Sickness=="Healthy")
str(healthy1)

H_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(healthy1),]
healthy1 <- healthy1[rownames(healthy1) %in% rownames(H_OTU),]

#Creating the Physeq object for the healthy animals
### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(H_OTU), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(healthy1)

physeq_H = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_H
# otu_table()   OTU Table:         [ 3513 taxa and 47 samples ]
colnames(tax_table(physeq_H))

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
dist.brayH <- phyloseq::distance(physeq_H, method = "bray")
x <- as.matrix(dist.brayH)
dist.brayH <- as.dist(dist.brayH)
write.csv(x, "dist.brayH.csv")

## PERMANOVA
#Bray curtis
BC_H <- adonis(dist.brayH ~ healthy1$d, strata=healthy1$calf, permutations = 999)
BC_H


## Weighted Unifrac
#we will use the phylos

##Making the tree, we need the phyloseq object 
random_tree = rtree(ntaxa(physeq_H), rooted=TRUE, tip.label=taxa_names(physeq_H))
##Merging the tree with the phyloseq object
physeqH1 = merge_phyloseq(physeq_H, meta.physeq, random_tree)
physeqH1 #This command should show the otu table, sample data, tax table and the phy tree

OTU_UnifracH <- UniFrac(physeqH1, weighted = TRUE) ## Calculating the weighted unifrac distances, this is a distance matrix
y <- as.matrix(OTU_UnifracH)
write.csv(y, "WU.H.csv")
## Permanova using the distance matrix
M3 <- adonis(OTU_UnifracH~ healthy1$d, strata=healthy1$calf, permutations = 999)
M3

#### Sick Animals
ant <- read.table("CalfMetadataSickAnimals.txt", header=TRUE, sep="\t")
ant$antibiotic <- factor(ant$antibiotic)
ant$treatment <- factor(ant$treatment)
ant$calf <- factor(ant$calf)
rownames(ant) <- ant$ID2

S_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(ant),] #33
ant <- ant[rownames(ant) %in% rownames(S_OTU),] #33

#Creating the Physeq object for the healthy animals
### Creating the Phyloseq Object
OTU.physeqS = otu_table(as.matrix(S_OTU), taxa_are_rows=FALSE)
tax.physeqS = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqS = sample_data(ant)

physeq_S = phyloseq(OTU.physeqS, tax.physeqS, meta.physeqS)
physeq_S
# otu_table()   OTU Table:         [ 3513 taxa and 33 samples ]
colnames(tax_table(physeq_S))

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
dist.brayS <- phyloseq::distance(physeq_S, method = "bray")
dist.brayS <- as.dist(dist.brayS)
x <- as.matrix(dist.brayS)
write.csv(x, "dist.brayS.csv")

## PERMANOVA
#Bray curtis
BC_S <- adonis(dist.brayS ~ ant$d, strata=ant$calf, permutations = 999)
BC_S #No significant

#Weighted Unifrac
#we will use the phylos

##Making the tree, we need the phyloseq object 
random_treeS = rtree(ntaxa(physeq_S), rooted=TRUE, tip.label=taxa_names(physeq_S))
##Merging the tree with the phyloseq object
physeqS1 = merge_phyloseq(physeq_S, meta.physeqS, random_treeS)
physeqS1 #This command should show the otu table, sample data, tax table and the phy tree

OTU_UnifracS <- UniFrac(physeqS1, weighted = TRUE) ## Calculating the weighted unifrac distances, this is a distance matrix
y <- as.matrix(OTU_UnifracS)
write.csv(y, "WU.S.csv")
## Permanova using the distance matrix
S3 <- adonis(OTU_UnifracS~ ant$d, strata=ant$calf, permutations = 999)
S3

####treatment and antibiotic effect each dat
## ---------------------------------------ds1, antibiotic and treatment
ant$day <- as.factor(ant$d)
days1 <- subset(ant, day == "1")
days1$calf <- as.factor(days1$calf)
ds1 <- NewASVTable[rownames(NewASVTable) %in% rownames(days1),]

OTU.physeq = otu_table(as.matrix(ds1), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(days1)

physeq_H = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_H
# otu_table()   OTU Table:         [ 3838 taxa and 9 samples ]

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
dist.brayH <- phyloseq::distance(physeq_H, method = "bray")
dist.brayH <- as.dist(dist.brayH)

## PERMANOVA ------ Day s1
#Bray curtis
BC_H1 <- adonis(dist.brayH ~ days1$antibiotic, permutations = 999) #Significant antibiotic
BC_H1 #ant

#Weighted UniFrac
random_tree = rtree(ntaxa(physeq_H), rooted=TRUE, tip.label=taxa_names(physeq_H))
##Merging the tree with the phyloseq object
physeqH1 = merge_phyloseq(physeq_H, meta.physeq, random_tree)
physeqH1 #This command should show the otu table, sample data, tax table and the phy tree

OTU_UnifracH <- UniFrac(physeqH1, weighted = TRUE) ## Calculating the weighted unifrac distances, this is a distance matrix

## Permanova using the distance matrix
M3 <- adonis(OTU_UnifracH~  days1$antibiotic,  permutations = 999) #significant antibiotic
M3 #ant


## Weighted unifrac 
ordu.wt.uni <- ordinate(physeqH1 , "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(physeqH1, 
                              ordu.wt.uni, color="antibiotic") 
wuaxis1 <- wt.unifrac[["data"]][["Axis.1"]]
wuaxis1 <- as.data.frame(wuaxis1)
wuaxis1$number <- rownames(wuaxis1)
wuaxis2 <- wt.unifrac[["data"]][["Axis.2"]]
wuaxis2 <- as.data.frame(wuaxis2)
wuaxis2$number <- rownames(wuaxis2)
wuaxis3 <- wt.unifrac[["data"]][["calf"]]
wuaxis3 <- as.data.frame(wuaxis3)
wuaxis3$number <- rownames(wuaxis3)
wuaxis4 <- wt.unifrac[["data"]][["antibiotic"]]
wuaxis4 <- as.data.frame(wuaxis4)
wuaxis4$number <- rownames(wuaxis4)

wuaxis <- merge(wuaxis3, wuaxis4, by.x = "number", by.y = "number")
wuaxis <- merge(wuaxis, wuaxis1, by.x = "number", by.y = "number")
wuaxis <- merge(wuaxis, wuaxis2, by.x = "number", by.y = "number")


wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() 
print(wt.unifrac)

##other plot
str(wuaxis)

centroids <- aggregate(wuaxis[,4:5], list(Group=wuaxis$wuaxis4), mean)
colnames(centroids) <- c('antibiotic','groupX', 'groupY')

wuaxis <- merge(wuaxis, centroids, by.x = "wuaxis4", by.y = "antibiotic")

str(wuaxis)
b <- ggplot(wuaxis, aes(x=wuaxis1, y=wuaxis2, color=wuaxis4)) + 
  geom_point(size=2) + 
  theme_classic() + #stat_ellipse() +
  guides(size=FALSE) +
  #geom_point(data= wuaxis, aes(x=groupX, y=groupX, shape=wuaxis4, size=5, fill=wuaxis4)) +
  scale_shape_manual(values=c(25, 22))+
  labs(color= "Antibiotic") +
  labs(fill= "Health Status") +
  labs(shape= "Centroids") +
  #geom_segment(data= wuaxis, aes(x=wuaxis1, y=wuaxis2, xend=groupX, yend=groupY, color= wuaxis4), size = .05) +
  labs(x='Axis 1 (37%)', y= 'Axis 2  (23.4%)', caption = paste('Distance between centroids: 0.0.269')) +
  ggtitle("Weighted UniFrac") +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 


#distance betweeen centroids
pointDistance(c(-0.08268166, -0.05782184), c(0.13780277, 0.09636973), lonlat=FALSE)
# 0.269051

## Bray_curtis
ordu.bc <- ordinate(physeq_H, "PCoA", "bray")
Bray <- plot_ordination(physeq_H, 
                        ordu.bc, color="antibiotic") 
Bray1 <- Bray[["data"]][["Axis.1"]]
Bray1 <- as.data.frame(Bray1)
Bray1$number <- rownames(Bray1)
Bray2 <- Bray[["data"]][["Axis.2"]]
Bray2 <- as.data.frame(Bray2)
Bray2$number <- rownames(Bray2)
Bray3 <- Bray[["data"]][["calf"]]
Bray3 <- as.data.frame(Bray3)
Bray3$number <- rownames(Bray3)
Bray4 <- Bray[["data"]][["antibiotic"]]
Bray4 <- as.data.frame(Bray4)
Bray4$number <- rownames(Bray4)

bray <- merge(Bray3, Bray4, by.x = "number", by.y = "number")
bray <- merge(bray, Bray1, by.x = "number", by.y = "number")
bray <- merge(bray, Bray2, by.x = "number", by.y = "number")

Bray <- Bray + ggtitle("Bray Curtis") + geom_point(size = 2)
Bray <- Bray + theme_classic()
print(Bray)

str(bray)

centroids2 <- aggregate(bray[,4:5], list(Group=bray$Bray4), mean)
colnames(centroids2) <- c('antibiotic','groupX', 'groupY')

bray <- merge(bray, centroids2, by.x = "Bray4", by.y = "antibiotic")
## Distance between centroids
pointDistance(c(-0.08733063, -0.05533778), c(0.14555105, 0.09222964), lonlat=FALSE)
#0.2756991

a <- ggplot(bray, aes(x=Bray1, y=Bray2, color=Bray4)) + 
  geom_point(size=2) + 
  theme_classic() + #stat_ellipse() +
  guides(size=FALSE) +
  #geom_point(data= bray, aes(x=groupX, y=groupX, shape=Bray4, size=5, fill=Bray4)) +
  scale_shape_manual(values=c(25, 22))+
  labs(color= "Antibiotic") +
  labs(shape= "Centroids") +
  labs(x='Axis 1 (32.5%)', y= 'Axis 2 (21.3%)', caption = paste('Distance between centroids: 0.276')) +
  ggtitle("Bray-Curtis dissimilarity") +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 

ggarrange(a,b,labels = c("a", "b"),
          ncol = 2, font.label = list(size = 20))

# Dispersion
D1 <- betadisper(dist.brayH, type = c("centroid"), group = days1$antibiotic)
D1
boxplot(D1)
permutest(D1, permutations = 999)#, pairwise = TRUE) #no difference 

D1 <- betadisper(OTU_UnifracH, type = c("centroid"), group = days1$antibiotic)
D1
boxplot(D1)
permutest(D1, permutations = 999)#, pairwise = TRUE) #no difference 

##------------------------------------- ds5, antibiotic and treatment
days5 <- subset(ant, day == "5")
days5$calf <- as.factor(days5$calf)
ds5 <- NewASVTable[rownames(NewASVTable) %in% rownames(days5),]

OTU.physeq = otu_table(as.matrix(ds5), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(days5)

physeq_H = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_H
# otu_table()   OTU Table:         [ 3513 taxa and 8 samples ]

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
dist.brayH <- phyloseq::distance(physeq_H, method = "bray")
dist.brayH <- as.dist(dist.brayH)

## PERMANOVA ------ Day s5
#Bray curtis
BC_H2 <- adonis(dist.brayH ~ days5$antibiotic, permutations = 999) #No significant
BC_H2 #NS

#Weighted UniFrac
random_tree = rtree(ntaxa(physeq_H), rooted=TRUE, tip.label=taxa_names(physeq_H))
##Merging the tree with the phyloseq object
physeqH1 = merge_phyloseq(physeq_H, meta.physeq, random_tree)
physeqH1 #This command should show the otu table, sample data, tax table and the phy tree

OTU_UnifracH <- UniFrac(physeqH1, weighted = TRUE) ## Calculating the weighted unifrac distances, this is a distance matrix

## Permanova using the distance matrix
M3 <- adonis(OTU_UnifracH~ days5$antibiotic,  permutations = 999) #no significant
M3 #NS

# Dispersion
D3 <- betadisper(dist.brayH, type = c("centroid"), group = days5$antibiotic)
D3
boxplot(D3)
permutest(D3, permutations = 999)#, pairwise = TRUE)  

D5 <- betadisper(OTU_UnifracH, type = c("centroid"), group = days5$antibiotic)
D5
boxplot(D5)
permutest(D5, permutations = 999)#, pairwise = TRUE) #no difference  


##-------------------------------- ds10, antibiotic and treatment
days10 <- subset(ant, day == "10")
days10$calf <- as.factor(days10$calf)
ds10 <- NewASVTable[rownames(NewASVTable) %in% rownames(days10),]

OTU.physeq = otu_table(as.matrix(ds10), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(days10)

physeq_H = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_H
# otu_table()   OTU Table:         [ 3513 taxa and 9 samples ]

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
dist.brayH <- phyloseq::distance(physeq_H, method = "bray")
dist.brayH <- as.dist(dist.brayH)

## PERMANOVA ------ Day s5
#Bray curtis
BC_H2 <- adonis(dist.brayH ~  days10$antibiotic, permutations = 999) #No significant
BC_H2 #NS

#Weighted UniFrac
random_tree = rtree(ntaxa(physeq_H), rooted=TRUE, tip.label=taxa_names(physeq_H))
##Merging the tree with the phyloseq object
physeqH1 = merge_phyloseq(physeq_H, meta.physeq, random_tree)
physeqH1 #This command should show the otu table, sample data, tax table and the phy tree

OTU_UnifracH <- UniFrac(physeqH1, weighted = TRUE) ## Calculating the weighted unifrac distances, this is a distance matrix

## Permanova using the distance matrix
M3 <- adonis(OTU_UnifracH~days10$antibiotic,  permutations = 999) #no significant
M3 #NS

# Dispersion
D5 <- betadisper(dist.brayH, type = c("centroid"), group = days10$antibiotic)
D5
boxplot(D5)
permutest(D5, permutations = 999)#, pairwise = TRUE) #No difference

# Dispersion
D5 <- betadisper(OTU_UnifracH, type = c("centroid"), group = days10$antibiotic)
D5
boxplot(D5)
permutest(D5, permutations = 999)#, pairwise = TRUE) #No difference

#### Sickness effect 
comp7_0 <- subset(metadata, subset=day %in% c("h7", "s0"))
comp7_0 <- comp7_0[-c(1,10,28 ),]
comp7_0 <- subset(comp7_0, subset=Antibiotic %in% c("No"))
comp7_0 %>% count(Sickness)
comp7_0 %>% count(day)

#### calculating beta diversity based on sickness
Si_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(comp7_0),]

OTU.physeq = otu_table(as.matrix(Si_OTU), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(comp7_0)

physeq_H = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_H
# otu_table()   OTU Table:         [ 3513 taxa and 20 samples ]

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
dist.brayH <- phyloseq::distance(physeq_H, method = "bray")
dist.brayH <- as.dist(dist.brayH)

## PERMANOVA ------ Day sickness effect
#Bray curtis
BC_H2 <- adonis(dist.brayH ~ comp7_0$Sickness, permutations = 999) #No significant
BC_H2 #NS

#Weighted UniFrac
random_tree = rtree(ntaxa(physeq_H), rooted=TRUE, tip.label=taxa_names(physeq_H))
##Merging the tree with the phyloseq object
physeqH1 = merge_phyloseq(physeq_H, meta.physeq, random_tree)
physeqH1 #This command should show the otu table, sample data, tax table and the phy tree

OTU_UnifracH <- UniFrac(physeqH1, weighted = TRUE) ## Calculating the weighted unifrac distances, this is a distance matrix

## Permanova using the distance matrix
M3 <- adonis(OTU_UnifracH~ comp7_0$Sickness,  permutations = 999) #no significant
M3 #NS

# Dispersion
D5 <- betadisper(dist.brayH, type = c("centroid"), group = comp7_0$Sickness)
D5
boxplot(D5)
permutest(D5, permutations = 999)#, pairwise = TRUE) #No difference

# Dispersion
D5 <- betadisper(OTU_UnifracH, type = c("centroid"), group = comp7_0$Sickness)
D5
boxplot(D5)
permutest(D5, permutations = 999)#, pairwise = TRUE) #No difference

#### test the antibiotic effect between d0 and d1, separate the sample based on antibiotic treatment
str(ant)
ant$day = as.factor(ant$day)
levels(ant$day)
antd0_1 = ant[ant$day %in% c('s0', 's1'), ]
antd0_1 =antd0_1[-c(1),] #removing sample 1

## prepare the samples
ds01 <- NewASVTable[rownames(NewASVTable) %in% rownames(antd0_1),]
antd0_1 <- antd0_1[rownames(antd0_1) %in% rownames(ds01),] #to make sure we have the same samples
str(antd0_1)
antd0_1$d = as.factor(antd0_1$d)

OTU.physeq = otu_table(as.matrix(ds01), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(antd0_1)

physeq_d01 = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_d01
# otu_table()   OTU Table:         [ 3513 taxa and 15 samples ]

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
dist.brayH <- phyloseq::distance(physeq_d01, method = "bray")
dist.brayH <- as.dist(dist.brayH)

## PERMANOVA ------ 
#Bray curtis
#BC_H2 <- adonis2(dist.brayH ~  antd0_1$d + antd0_1$antibiotic, strata=antd0_1$calf, permutations=999) #No significant
BC_H2 #NS



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

##Pairwiese adonis function
#we can also performe a pairwise comparison with the function 
# Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
#https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

