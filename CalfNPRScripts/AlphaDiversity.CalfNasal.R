## Repeated measures design 
# https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/
yes
#install.packages("ggfortify")
library(stringi)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(emmeans)
library(lme4)
library(lattice)
library("latticeExtra")
library(dplyr)
library(summarytools)
library(ggfortify)

setwd("~/Desktop/eunice/PhD/CalfNasal/Qiime/CorrectFiltered/") #sets new working directory for Windows systems (remember to replace … with your filepath)
meta <- read.csv("CalfMetadata08.30.csv")
rownames(meta) <- meta$ID

otu_table <- read.table("observed_otus.tsv", header=TRUE, row.names=1, sep="\t")
shannon <-read.table("shannon.tsv", header=TRUE, row.names=1, sep="\t")
faith <- read.table("faith_pd.tsv", header=TRUE, row.names=1, sep="\t")
alpha_diversity <- merge(otu_table, shannon, by.x = 0, by.y = 0)
alpha_diversity <- merge(alpha_diversity, faith, by.x = "Row.names", by.y = 0)
meta <- merge(meta, alpha_diversity, by.x = 0, by.y = "Row.names")
row.names(meta) <- meta$ID2
str(meta)
meta$calf <- as.factor(meta$calf)
meta$Sickness <- as.factor(meta$Sickness)
levels(meta$Sickness)
meta %>% count(Sickness)
meta$Antibiotic <- factor(meta$Antibiotic)
#meta$d <- factor(meta$d)
meta$diagnostic <- factor(meta$diagnostic)
levels(meta$diagnostic)
levels(meta$Sickness) <- list("Sick"="Got", "Healthy"="Never")

### healthy animals repeated measures
str(meta$Sickness)
healthy <- subset(meta, Sickness=="Healthy")
str(healthy) #47 samples
healthy$Sickness <- factor(healthy$Sickness)
healthy$Antibiotic <- factor(healthy$Antibiotic)
healthy$diagnostic <- factor(healthy$diagnostic)
healthy$calf <- factor(healthy$calf)
str(healthy)
str(healthy)
#write.table(healthy,"Healthy.txt",sep=",", row.names = TRUE) 

#### Statistical OTUs and Evenness
set_sum_contrasts() # important for afex
levels(healthy$calf) #20
# we need to add the withing factor in the error term
m1 <- mixed(observed_otus ~ d+ (d|calf), data = healthy,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)))
m1
summary(m1)$varcor
m2 <- mixed(observed_otus ~ d+ (d||calf), data = healthy,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
m2
summary(m2)$varcor # without the correlations the SD for day goes down
 
#compare both models
left_join(nice(m1), nice(m2), by = "Effect", 
          suffix = c("_full", "_final")) 

plot(m2$full_model)
qqnorm(residuals(m2$full_model))
qqline(residuals(m2$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
with(healthy, shapiro.test(observed_otus))

ggplot(data = healthy, aes(x = d, y = observed_otus)) + 
  geom_smooth(method = "lm", se = FALSE) +  
  geom_point() + geom_smooth(method = "lm", se=TRUE) + 
  scale_color_brewer(palette = "Dark2") + theme_classic() + xlab ("Day") + 
  ylab ("Observed ASVs") +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) 
  

ggsummarystats(
  healthy, x = "d", y = "observed_otus", 
  ggfunc = ggboxplot, add = "jitter"
)

m3 <- mixed(faith_pd ~ d+ (d|calf), data = healthy,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)))
summary(m3)$varcor # better model
m4 <- mixed(faith_pd ~ d+ (d||calf), data = healthy,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
summary(m4)$varcor
plot(m3$full_model)
qqnorm(residuals(m3$full_model))
with(healthy, shapiro.test(faith_pd))

#compare both models
left_join(nice(m3), nice(m4), by = "Effect", 
          suffix = c("_full", "_final")) 

ggplot(data = healthy, aes(x = d, y = faith_pd)) + 
  geom_smooth(method = "lm", se = FALSE) +  
  geom_point() + geom_smooth(method = "lm", se=TRUE) + 
  scale_color_brewer(palette = "Dark2") + theme_classic() + xlab ("Day") + 
  ylab ("Faith's PD Phylotegetic relationship")

### data visualziation sick animals
# Data filtering 15 Ge feature show
ant <- read.table("CalfMetadataSickAnimals.txt", header=TRUE, sep="\t")
ant <- na.omit(ant)
ant$antibiotic <- factor(ant$antibiotic)
ant$calf <- factor(ant$calf)
ant$d 

#### Statistical OTUs and Evenness
set_sum_contrasts() # important for afex
levels(ant$calf)

y1 <- mixed(observed_otus ~ d + antibiotic+ d*antibiotic+ (d|calf), data = ant,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)))
summary(y1)$varcor 
y2 <- mixed(observed_otus ~ d + antibiotic+ d*antibiotic+ (d||calf), data = ant,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
summary(y2)$varcor # better model

plot(y2$full_model)
qqnorm(residuals(y2$full_model))
qqline(residuals(y2$full_model))
with(ant, shapiro.test(observed_otus))

#compare both models
left_join(nice(y1), nice(y2), by = "Effect", 
          suffix = c("_full", "_final")) 

ggplot(data = ant, aes(x = d, y = observed_otus)) + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(antibiotic~.) +
  geom_point() + geom_smooth(method = "lm", se=TRUE) + 
  scale_color_brewer(palette = "Dark2") + theme_bw() + xlab ("Day") + 
  ylab ("Observed ASVs")


y3 <- mixed(faith_pd ~ d+ antibiotic + d*antibiotic+(d|calf), data = ant,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)))
summary(y3)$varcor #better model
y4 <- mixed(faith_pd ~ d+ antibiotic + d*antibiotic+ (d||calf), data = ant,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
summary(y4)$varcor 

plot(y4$full_model)
qqnorm(residuals(y4$full_model))
qqline(residuals(y4$full_model))
with(ant, shapiro.test(faith_pd))

#compare both models
left_join(nice(y3), nice(y4), by = "Effect", 
          suffix = c("_full", "_final")) 


## Comparing healthy and BRD animals
## We will use the day h7 and day s0
comp7_0 <- subset(meta, subset=day %in% c("h7", "s0"))
comp7_0 <- comp7_0[-c(1,16,23 ),]
comp7_0 <- subset(comp7_0, subset=Antibiotic %in% c("No"))
comp7_0 %>% count(Sickness)
comp7_0 %>% count(d)
#Running normal anova because there are no repeated measures
compASV7 <- lm(observed_otus ~ Sickness, data = comp7_0) 

#t test 
t.test(observed_otus ~ Sickness, data = comp7_0)
t.test(faith_pd ~ Sickness, data = comp7_0)

