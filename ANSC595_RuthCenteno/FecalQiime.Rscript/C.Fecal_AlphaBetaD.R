#####            Qiime to R
#### Analysis of calf fecal samples alpha and beta diversity
###  Measuring the effect of sickness (animals diagnosed with BRD and not)
###  Measuring the effect of antibiotic 
###         Author: Ruth Eunice Centeno Martinez


####  Setting the working directory and the packges needed
setwd("~/Desktop/eu/BRD/qiime/ANSC595_RuthCenteno/FecalQiime.Rscript/")
library(ggplot2)
library(ggfortify) #for autoplot
library(summarytools) #for "descr" function
library(car) ##for type II SS 
library(dplyr) # for tally function
library(ggpubr) #for ggarrange function
library(vegan) #for the ordiellipse function


# 1. ---------- ALPHA DIVERSITY ANALYSIS -----
# All the data used in the script was generated using Qiime2 

#Metadata. 
metadata <- read.table("CalfFecal_metadata.new.txt", header=TRUE, row.names=1, sep="\t") # contains all the information from the samples sent to sequenced
### the analysis of sickness in the nasopharyngeal samples was performed only analyzing 16 calves (n=8 healthy animals vs n=8 sick animals or diagnosed with BRD)
animals <- read.table("Sickness.effect.fecal/animals.fecal.txt", header=TRUE, sep="\t")# indicating day in the study and day block
meta <-  merge(metadata, animals, by.x = 0, by.y = "ID")

faith_pd <- read.table("Sickness.effect.fecal/faith_pd.tsv", header=TRUE, row.names=1, sep="\t")
shannon <-read.table("Sickness.effect.fecal/shannon.tsv", header=TRUE, row.names=1, sep="\t")
alpha_diversity <- merge(faith_pd, shannon, by.x = 0, by.y = 0)
meta <- merge(meta, alpha_diversity, by.x = "Row.names", by.y = "Row.names")
row.names(meta) <- meta$Row.names
meta <- meta[,-1]

###Alpha Diversity tables
##merge table with healthy animals and sick animals
animals <- meta
str(animals)
animals$day.block <- factor(animals$day.block)
animals$day.study <- factor(animals$day.study)
animals$sickness <- factor(animals$sickness)
animals$antibiotic <- factor(animals$antibiotic)
animals$diagnostic <- factor(animals$diagnostic)
animals$calf <- factor(animals$calf)
animals$day <- factor(animals$day)
animals$treatment <- factor(animals$treatment)
str(animals)
animals$day.block <- factor(animals$day.block)
animals$day.study <- factor(animals$day.study)

#To test for normalcy statistically, we can run the Shapiro-Wilk test of normality.
#if the value is less than 0.05, then the data is not normally distributed
shapiro.test(animals$shannon)
#data:  animals$shannon
#W = 0.95177, p-value = 0.5183
shapiro.test(animals$faith_pd) #not normal distributed
#data:  animals$faith_pd
#W = 0.86437, p-value = 0.02233

#Plots
hist(animals$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
hist(animals$shannon, main="Shannon", xlab="", breaks=15)

##We check the homogeneity of variance and normalitiy assumption 
#we first create the glm model
#SHANNON
resp_shan <- lm(shannon ~ sickness + day.block, data = animals) 

# we get for each continuous variable: mean, SD, CV, skewness, kurtosis, histogram, and normal probability plot. 
# this is just for continuous variables - please, check the mean, SD, CV, Skewness, Kurtosis, N valid
# Continuous dependent variables: Chao and Shannon
descr(animals$shannon, style="rmarkdown")
#**Skewness** |    0.16 |
#|    **Kurtosis** |   -1.09 |

descr(animals$faith_pd, style="rmarkdown")
#|    **Skewness** |     1.35 |
#|    **Kurtosis** |     1.63 |

# normal probability plots (to establish visually deviations from normality (run the following two lines of code together - and for each continuous factor)
qqnorm(animals$shannon) ##run at the same time the two codes
qqline(animals$shannon)

qqnorm(animals$faith_pd) ##run at the same time the two codes
qqline(animals$faith_pd)

## Checking assumptions
##checking for assumptions: Homogeneity of variance and normality of the residuals 
autoplot(resp_shan, smooth.color = NA) 


##---------------- ANOVA ANALYSIS------------- # 
str(animals)
#chao and shannon are number
animals %>% tally()
animals %>% count(sickness)
#Got         8
#Never       8
animals %>% count(day.block)
#0             8
#7             6
#14            2
#the design is balance, not including the effect of the blocking factor only including the effect of Sickness

## ANOVA TYPE II SUM OF SQUARE
#you used the lm model
options(contrasts=c("contr.sum", "contr.poly"))

###----------------------------CHAO ---------------------# 
# DEFAULT ANOVA
#Create the model
#--------------- CHAO -------------
#1. OTUS COUNTS

#Anova type II SS--no interaction is assumed (only one effect and blocking factor)
respshan <- Anova(resp_shan, type = 2) ###significant
respshan

#### calculating effect size, partial-eta square 
#It’s the sums of squares for each effect divided by the error SS
respshan
respshan$pes <- c(respshan$'Sum Sq'[-nrow(respshan)], NA)/(respshan$'Sum Sq' + respshan$'Sum Sq'[nrow(respshan)])
respshan

## Stats for non-parametric data
sickant<-kruskal.test(faith_pd ~ sickness, data=animals)
sickant

#Visualization
#Sickness factor, this calculates the error bars
SHA_summarySick <- animals %>% # the names of the new data frame and the data frame to be summarised
  group_by(sickness) %>%   # the grouping variable
  summarise(mean_sha = mean(shannon),  # calculates the mean of each group
            sd_sha = sd(shannon), # calculates the standard deviation of each group
            n_sha = n(),  # calculates the sample size per group
            se_sha = sd(shannon)/sqrt(n())) # calculates the standard error of each group

Faith_summarySick <- animals %>% # the names of the new data frame and the data frame to be summarised
  group_by(sickness) %>%   # the grouping variable
  summarise(mean_faith = mean(faith_pd),  # calculates the mean of each group
            sd_faith = sd(faith_pd), # calculates the standard deviation of each group
            n_faith = n(),  # calculates the sample size per group
            se_faith = sd(faith_pd)/sqrt(n())) # calculates the standard error of each group

#PLOTS
#Sicknes
my_colors2 <- c('tan3',"turquoise4", 'skyblue1') 
a <- ggplot(data = SHA_summarySick, aes(x = sickness, y = mean_sha, fill = sickness)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  scale_fill_manual(values = my_colors2) +
  geom_errorbar(aes(ymin = mean_sha - se_sha, ymax = mean_sha + se_sha), width=0.2) +
  guides(fill=FALSE) +
  xlab("Sickness") + ylab("Shannon") + 
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  #theme(strip.background = element_text(color="black", size=10,face= "bold")) +
  theme(axis.title.y = element_text(color="black", size=10, face="bold")) + 
  #theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_text(color = "black", size = 10, face="bold"), axis.ticks = element_blank())

b <- ggplot(data = Faith_summarySick, aes(x = sickness, y = mean_faith, fill = sickness)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  scale_fill_manual(values = my_colors2) +
  geom_errorbar(aes(ymin = mean_faith - se_faith, ymax = mean_faith + se_faith), width=0.2) +
  guides(fill=FALSE) +
  xlab("Sickness") + ylab("Faith") + 
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  #theme(strip.background = element_text(color="black", size=10,face= "bold")) +
  theme(axis.title.y = element_text(color="black", size=10, face="bold")) + 
  #theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_text(color = "black", size = 10, face="bold"), axis.ticks = element_blank())

ggarrange(a, b, labels = c("A", "B"),
          ncol = 2)

#####-----ANTIBIOTIC EFFECT
#animals d7 vs ds5
Anti <- read.table("Antibiotic.effect.fecal/AntibioticFecal.txt", header=TRUE, sep="\t")
#this effect was only analyzed in 8 calves (n=8 healthy calves vs n=8 sick calves)
A.animals2 <-  merge(metadata, Anti, by.x = 0, by.y = "ID")
faith_pdA <- read.table("Antibiotic.effect.fecal/faith_pd.tsv", header=TRUE, row.names=1, sep="\t")
shannonA <-read.table("Antibiotic.effect.fecal/shannon.tsv", header=TRUE, row.names=1, sep="\t")
alpha_diversityA <- merge(faith_pdA, shannonA, by.x = 0, by.y = 0)
A.animals2 <- merge(A.animals2, alpha_diversityA, by.x = "Row.names", by.y = "Row.names")
row.names(A.animals2) <- A.animals2$Row.names
A.animals2 <- A.animals2[,-1]
str(A.animals2)

A.animals2$day.block <- factor(A.animals2$day.block)
A.animals2$day.study <- factor(A.animals2$day.study)
A.animals2$sickness <- factor(A.animals2$sickness)
A.animals2$antibiotic <- factor(A.animals2$antibiotic)
A.animals2$diagnostic <- factor(A.animals2$diagnostic)
A.animals2$calf <- factor(A.animals2$calf)
A.animals2$day <- factor(A.animals2$day)
A.animals2$treatment <- factor(A.animals2$treatment)

#To test for normalcy statistically, 
#we can run the Shapiro-Wilk test of normality.
#value less than <0.05, not normal
shapiro.test(A.animals2$shannon)
#data:  A.animals2$shannon
#W = 0.92947, p-value = 0.5113

shapiro.test(A.animals2$faith_pd)
#data:  A.animals2$faith_pd
#W = 0.89667, p-value = 0.2696

#histograms for continuous factors
ggplot(A.animals2, aes(x= shannon)) + geom_histogram(binwidth=0.4) ###no following a normal distribution
ggplot(A.animals2, aes(x= faith_pd)) + geom_histogram(binwidth=0.4) ###no following a normal distribution

# normal probability plots (to establish visually deviations from normality (run the following two lines of code together - and for each continuous factor)
qqnorm(A.animals2$shannon) ##run at the same time the two codes
qqline(A.animals2$shannon)

qqnorm(A.animals2$faith_pd) 
qqline(A.animals2$faith_pd)

ant_sha <- lm(shannon ~ antibiotic, data = A.animals2)
ant_faith <- lm(faith_pd ~ antibiotic, data = A.animals2)

## ANOVA TYPE II SUM OF SQUARE
A.animals2 %>% tally()
A.animals2 %>% count(antibiotic)
#No             4
#Yes            4
#the model is balance

#you used the lm model
options(contrasts=c("contr.sum", "contr.poly"))

###----------------------------CHAO ---------------------# 
# DEFAULT ANOVA
#Create the model
#--------------- CHAO -------------
#1. OTUS COUNTS

#Anova type II SS-- no interaction due to blocking facter
antsha <- Anova(ant_sha, type = 2)
antsha

#### calculating effect size, partial-eta square 
#It’s the sums of squares for each effect divided by the error SS
antsha
antsha$pes <- c(antsha$'Sum Sq'[-nrow(antsha)], NA)/(antsha$'Sum Sq' + antsha$'Sum Sq'[nrow(antsha)])
antsha


antfaith <- Anova(ant_faith, type = 2)
antfaith
#### calculating effect size, partial-eta square 

antfaith
antfaith$pes <- c(antfaith$'Sum Sq'[-nrow(antfaith)], NA)/(antfaith$'Sum Sq' + antfaith$'Sum Sq'[nrow(antfaith)])
antfaith

#Antibiotic factor, this calculates the error bars
SHA_summaryA <- A.animals2 %>% # the names of the new data frame and the data frame to be summarised
  group_by(antibiotic) %>%   # the grouping variable
  summarise(mean_sha = mean(shannon),  # calculates the mean of each group
            sd_sha = sd(shannon), # calculates the standard deviation of each group
            n_sha = n(),  # calculates the sample size per group
            se_sha = sd(shannon)/sqrt(n())) # calculates the standard error of each group

Faith_summaryA <- A.animals2 %>% # the names of the new data frame and the data frame to be summarised
  group_by(antibiotic) %>%   # the grouping variable
  summarise(mean_faith = mean(faith_pd),  # calculates the mean of each group
            sd_faith = sd(faith_pd), # calculates the standard deviation of each group
            n_faith = n(),  # calculates the sample size per group
            se_faith = sd(faith_pd)/sqrt(n())) # calculates the standard error of each group

#PLOTS
#Sicknes
my_colors <- c("lightpink", "mediumorchid3") 
x <- ggplot(data = SHA_summaryA, aes(x = antibiotic, y = mean_sha, fill = antibiotic)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  scale_fill_manual(values = my_colors) +
  geom_errorbar(aes(ymin = mean_sha - se_sha, ymax = mean_sha + se_sha), width=0.2) +
  guides(fill=FALSE) +
  xlab("Antibiotic") + ylab("Shannon") + 
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  #theme(strip.background = element_text(color="black", size=10,face= "bold")) +
  theme(axis.title.y = element_text(color="black", size=10, face="bold")) + 
  #theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_text(color = "black", size = 10, face="bold"), axis.ticks = element_blank())

y <- ggplot(data = Faith_summaryA, aes(x = antibiotic, y = mean_faith, fill = antibiotic)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  scale_fill_manual(values = my_colors) +
  geom_errorbar(aes(ymin = mean_faith - se_faith, ymax = mean_faith + se_faith), width=0.2) +
  guides(fill=FALSE) +
  xlab("Antibiotic") + ylab("Faith") + 
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  #theme(strip.background = element_text(color="black", size=10,face= "bold")) +
  theme(axis.title.y = element_text(color="black", size=10, face="bold")) + 
  #theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_text(color = "black", size = 10, face="bold"), axis.ticks = element_blank())

ggarrange(x, y, labels = c("C", "D"),
          ncol = 2)

###----------------2. BETA DIVERSITY---------------------------
## All the information related to PCoA coordinates were obtained from Qiim2

##### SICKNESS EFFECT-------------------------------

#Metadata. 
meta #The same 16 animals analyzed in the alpha diversity measuring the effect of sickness

## ---- UNWEIGHTED UNIFRAC PCoA Sickness
unifracUN <- read.table("Sickness.effect.fecal//UnwUniFecal.txt", header=TRUE, sep="\t", stringsAsFactors = F) ##Unweighted UniFrac Coordinates
meta.pcoaSUN<-merge(meta,unifracUN, by.x = 0, by.y="ID")
str(meta.pcoaSUN)
meta.pcoaSUN$Row.names <- factor(meta.pcoaSUN$Row.names)
meta.pcoaSUN$day <- factor(meta.pcoaSUN$day)
meta.pcoaSUN$calf <- factor(meta.pcoaSUN$calf)
meta.pcoaSUN$day.block <- factor(meta.pcoaSUN$day.block)
meta.pcoaSUN$day.study <- factor(meta.pcoaSUN$day.study)
meta.pcoaSUN$sickness <- factor(meta.pcoaSUN$sickness)
meta.pcoaSUN$antibiotic <- factor(meta.pcoaSUN$antibiotic)
meta.pcoaSUN$diagnostic <- factor(meta.pcoaSUN$diagnostic)
str(meta.pcoaSUN)

row.names(meta.pcoaSUN) <- meta.pcoaSUN$Row.names
meta.pcoaSUN <-meta.pcoaSUN[,-1]
pcoa.sickUN <- meta.pcoaSUN[,-c(2:4,6,7,9,10)]
str(pcoa.sickUN)

#creating a variable for the colors
colorU <- c('tan3',"turquoise4", 'skyblue1')
colsU<-colorU[pcoa.sickUN$sickness]

par(mar = c(4,4,1,2))
plot(pcoa.sickUN[,4:5], pch=17,cex=0.9, col=colsU, xlab="PCoA1 (16.78%)", ylab="PCoA2 (13.11%)", xlim=c(-0.5, 0.7), ylim=c(-0.6, 0.6))
legend(0.45, 0.64,legend=levels(pcoa.sickUN$sickness), col=colorU, pch=17, cex=0.75)
ordiellipse(pcoa.sickUN[,4:5], pcoa.sickUN$sickness, col=colorU, conf = 0.95)
text(pcoa.sickUN[,4:5],labels=pcoa.sickUN$day.block, cex=0.7,pos=3)


## ---- WEIGHTED UNIFRAC pca Sickness
meta #The same 16 animals analyzed in the alpha diversity measuring the effect of sickness

unifracwe <- read.table("Sickness.effect.fecal/pcoaw.fecal.txt", header=TRUE, sep="\t", stringsAsFactors = F) ##Weighted UniFrac Coordinates
meta.pcoaS<-merge(meta,unifracwe, by.x = 0, by.y="ID")
str(meta.pcoaS)
meta.pcoaS$Row.names <- factor(meta.pcoaS$Row.names)
meta.pcoaS$day <- factor(meta.pcoaS$day)
meta.pcoaS$calf <- factor(meta.pcoaS$calf)
meta.pcoaS$day.block <- factor(meta.pcoaS$day.block)
meta.pcoaS$day.study <- factor(meta.pcoaS$day.study)
meta.pcoaS$sickness <- factor(meta.pcoaS$sickness)
meta.pcoaS$antibiotic <- factor(meta.pcoaS$antibiotic)
meta.pcoaS$diagnostic <- factor(meta.pcoaS$diagnostic)
str(meta.pcoaS)

row.names(meta.pcoaS) <- meta.pcoaS$Row.names
meta.pcoaS <-meta.pcoaS[,-1]
pcoa.sick <- meta.pcoaS[,-c(2:4,6,7,9,10)]
str(pcoa.sick)

#creating a variable for the colors
color <- c('tan3',"turquoise4", 'skyblue1')
cols<-color[pcoa.sick$sickness]

par(mar = c(4,4,1,2))
plot(pcoa.sick[,4:5], pch=17,cex=0.9, col=cols, xlab="PCoA1 (22.24%)", ylab="PCoA2 (12.24%)", xlim=c(-0.4, 0.5), ylim=c(-0.4, 0.5))
legend(0.32, 0.52,legend=levels(pcoa.sick$sickness), col=color, pch=17, cex=0.75)
ordiellipse(pcoa.sick[,4:5], pcoa.sick$sickness, col=color, conf = 0.95)
text(pcoa.sick[,4:5],labels=pcoa.sick$day.block, cex=0.7,pos=3)


##### UNWEIGHTED ANTIBIOTIC EFFECT--------------------------------------

###metadata and OTU table
#These tables will be merged for convenience and added to the metadata table as the original tutorial was organized.
A.animals2 ##the same 8 analyzed in the alpha diversity
unifracA.UN <- read.table("Antibiotic.effect.fecal/UNwF.AntFecal.txt", header=TRUE, sep="\t", stringsAsFactors = F)
meta.pcoaAUN<-merge(A.animals2,unifracA.UN, by.x = 0, by.y="ID")
str(meta.pcoaAUN)
meta.pcoaAUN$Row.names <- factor(meta.pcoaAUN$Row.names)
meta.pcoaAUN$day <- factor(meta.pcoaAUN$day)
meta.pcoaAUN$calf <- factor(meta.pcoaAUN$calf)
meta.pcoaAUN$day.block <- factor(meta.pcoaAUN$day.block)
meta.pcoaAUN$day.study <- factor(meta.pcoaAUN$day.study)
meta.pcoaAUN$sickness <- factor(meta.pcoaAUN$sickness)
meta.pcoaAUN$antibiotic <- factor(meta.pcoaAUN$antibiotic)
meta.pcoaAUN$diagnostic <- factor(meta.pcoaAUN$diagnostic)
str(meta.pcoaAUN)

row.names(meta.pcoaAUN) <- meta.pcoaAUN$Row.names
meta.pcoaAUN <-meta.pcoaAUN[,-1]
pcoa.A.UN <- meta.pcoaAUN[,-c(2:5,7,9,10)]
str(pcoa.A.UN)

#creating a variable for the colors
colorAU <- c("lightpink", "mediumorchid3")
colsAU<-colorAU[pcoa.A.UN$antibiotic]

par(mar = c(4,4,1,2))
plot(pcoa.A.UN[,4:5], pch=17,cex=0.9, col=colsAU,xlab="PCoA1 (18.99%)", ylab="PCoA2 (11.66%)", xlim=c(-0.65, 0.72), ylim=c(-0.6, 0.7))
legend(0.48, 0.75,legend=levels(pcoa.A.UN$antibiotic), col=colorAU, pch=17, cex=0.75)
ordiellipse(pcoa.A.UN[,4:5], pcoa.A.UN$antibiotic, col=colorAU, conf = 0.95)
text(pcoa.A.UN[,4:5],labels=pcoa.A.UN$day.block, cex=0.7,pos=1)

##-----   WIEGHTED UNIFRAC pca Antibiotic
A.animals2 ##the same 8 analyzed in the alpha diversity
unifracwe <- read.table("Antibiotic.effect.fecal/WUF.AntFecal.txt", header=TRUE, sep="\t", stringsAsFactors = F)
meta.pcoaA<-merge(A.animals2,unifracwe, by.x = 0, by.y="ID")
str(meta.pcoaA)
meta.pcoaA$Row.names <- factor(meta.pcoaA$Row.names)
meta.pcoaA$day <- factor(meta.pcoaA$day)
meta.pcoaA$calf <- factor(meta.pcoaA$calf)
meta.pcoaA$day.block <- factor(meta.pcoaA$day.block)
meta.pcoaA$day.study <- factor(meta.pcoaA$day.study)
meta.pcoaA$sickness <- factor(meta.pcoaA$sickness)
meta.pcoaA$antibiotic <- factor(meta.pcoaA$antibiotic)
meta.pcoaA$diagnostic <- factor(meta.pcoaA$diagnostic)
str(meta.pcoaA)

row.names(meta.pcoaA) <- meta.pcoaA$Row.names
meta.pcoaA <-meta.pcoaA[,-1]
pcoa.ant <- meta.pcoaA[,-c(2:5,7,9,10)]
str(pcoa.ant)

#creating a variable for the colors
color2 <- c("lightpink", "mediumorchid3")
cols2<-color2[pcoa.ant$antibiotic]

par(mar = c(4,4,1,2))
plot(pcoa.ant[,4:5], pch=17, col=cols2, xlab="PCoA1 (24.98%)", ylab="PCoA2 (18.66%)", xlim=c(-0.6, 0.7), ylim=c(-0.5, 0.6))
legend(0.48, 0.64,legend=levels(pcoa.ant$antibiotic), col=color2, pch=17, cex=0.75)
ordiellipse(pcoa.ant[,4:5], pcoa.ant$antibiotic, col=color2, conf = 0.95)
text(pcoa.ant[,4:5],labels=pcoa.ant$day.block, cex=0.7,pos=2)








