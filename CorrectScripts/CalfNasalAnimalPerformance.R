#### Calf Nasal animal performance

library(dplyr)
library(broom)
library(ggplot2)
library(ggpubr)

setwd("~/Desktop/eunice/PhD/CalfNasal/Qiime/CorrectFiltered/")
#Animal performance data
data2 <- read.csv("CalfNasalAnimalPerformancedata.csv")
data2$Sickness <- as.factor(data2$Sickness)
data2 = data2[-c(1),] #remove calf 1
str(data2)

# test normality 
shapiro.test(data2$ADG..kg.d.)
shapiro.test(data2$MR.Intake..pints..)
shapiro.test(data2$Starter.Intake..g.)

summary_ADG..kg.d. <-data2 %>%
  group_by(Sickness) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(ADG..kg.d.),
    SD = sd(ADG..kg.d.),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
test1 <- aov(ADG..kg.d. ~ Sickness, data = data2)
summary(test1)

A <- ggplot(data2, aes(x= Sickness, y = ADG..kg.d., fill =Sickness)) +
  geom_boxplot() + theme_bw() +
  #geom_errorbar(data=summary_ADG..kg.d., aes(x=Sickness, ymin=Mean - SD,
                                             #ymax=Mean + SD, y=NULL), color="black", width=0.2) +
  geom_jitter(color='gray', width = 0.25) +
  geom_point(data= summary_ADG..kg.d., aes(Sickness, y=Mean), color="black", shape= 17, size=3) +
  ylab("ADG (kg/d)") +xlab ("Sickness") +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text(size = 14, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 14), axis.text.y =
          element_text(color = "black", size = 14))

summary_MR.Intake..pints.. <-data2 %>%
  group_by(Sickness) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(MR.Intake..pints..),
    SD = sd(MR.Intake..pints..),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
test2<- wilcox.test(MR.Intake..pints..~ Sickness, data = data2)
test2

B <- ggplot(data2, aes(x= Sickness, y = MR.Intake..pints.., fill =Sickness)) +
  geom_boxplot() + theme_bw() +
  #geom_errorbar(data=summary_ADG..kg.d., aes(x=Sickness, ymin=Mean - SD,
  #ymax=Mean + SD, y=NULL), color="black", width=0.2) +
  geom_jitter(color='gray', width = 0.25) +
  geom_point(data= summary_MR.Intake..pints.., aes(Sickness, y=Mean), color="black", shape=17, size=3) +
  ylab("MR.Intake (Pints)") +xlab ("Sickness") +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text(size = 14, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 14), axis.text.y =
          element_text(color = "black", size = 14))

summary_Starter.Intake..g. <-data2 %>%
  group_by(Sickness) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(Starter.Intake..g.),
    SD = sd(Starter.Intake..g.),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE)
  )
test3<- wilcox.test(Starter.Intake..g.~ Sickness, data = data2)
test3

C <- ggplot(data2, aes(x= Sickness, y = Starter.Intake..g., fill =Sickness)) +
  geom_boxplot() + theme_bw() +
  #geom_errorbar(data=summary_ADG..kg.d., aes(x=Sickness, ymin=Mean - SD,
  #ymax=Mean + SD, y=NULL), color="black", width=0.2) +
  geom_jitter(color='gray', width = 0.25) +
  geom_point(data= summary_Starter.Intake..g., aes(Sickness, y=Mean), color="black", shape=17, size=3) +
  ylab("Started Intake (g)") +xlab ("Sickness") +
  theme(strip.text = element_text(size = 14, face = "bold")) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text(size = 14, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 14), axis.text.y =
          element_text(color = "black", size = 14))

ggarrange(A,B,C, labels = c("a", "b", "c"),
          ncol = 3, font.label = list(size = 20))
