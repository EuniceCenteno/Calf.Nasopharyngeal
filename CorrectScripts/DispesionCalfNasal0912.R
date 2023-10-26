#Betadisper Nasal Samples
library(vegan) 
library(ggplot2)
library(ggpubr)
library(dplyr) #tally function
library(data.table) #transpose function
library(car) ##for type III SS 
library(lme4)
library(afex)
library(emmeans)

setwd("~/Desktop/eunice/PhD/CalfNasal/Qiime/CorrectFiltered/") #sets new working directory for Windows systems (remember to replace … with your filepath)
metadata1 <-read.csv("CalfMetadata08.30.csv")
metadata1$calf <- as.factor(metadata1$calf)
metadata1$Sickness <- as.factor(metadata1$Sickness)
metadata1$day <- factor(metadata1$day)
metadata1$diagnostic <- factor(metadata1$diagnostic)
levels(metadata1$Sickness) <- list("Sick"="Got", "Healthy"="Never")
metadata1$Antibiotic <- factor(metadata1$Antibiotic)

##Subsetting healthy and sick animals
healthy1 <- subset(metadata1, Sickness=="Healthy")
str(healthy1)
healthy1$Sickness <- factor(healthy1$Sickness)
healthy1$diagnostic <- factor(healthy1$diagnostic)
healthy1$d <- factor(healthy1$d)
levels(healthy1$d) <- list("0"="0", "7"="7", "14"="14")
healthy1$calf <- factor(healthy1$calf)
rownames(healthy1) <- healthy1$ID2
#healthy$d <- as.numeric(as.character(healthy$d))
levels(healthy1$calf)
rownames(healthy1) <- healthy1$ID

#We need to import the distance matrices
DistWU = read.csv("WU.H.csv", header=TRUE)
row.names(DistWU) <- DistWU[,1]
DistWU = DistWU[,-c(1)]
distWU <- as.dist(DistWU)
#check homogeneity of variance
healthy1 <-healthy1[rownames(healthy1) %in% rownames(DistWU),]#47

time <- betadisper(distWU, type = c("centroid"), group = healthy1$d)
time
plot(time)
distances <- time[["distances"]]
distances <- as.data.frame(distances)
distances <- merge(distances, healthy1, by.x = 0, by.y=0) # we are combining the data and removing unneccesary data
distances <- distances[,-c(1, 5, 7:9)] # here I need to calculate by hand the mean of the groups
write.csv(distances, "distancesHWU.csv")

#now let's try a different approach which is calculated in correction factor
# you sum all the values for d0, d7 and d14 and divided by the total number of observations in each dat
metadata <- read.csv("betadisperCalfN.csv", header=TRUE)
x <- metadata$d0
x <- na.exclude(x) 
mean0 <- mean(x)
x <- sum(x)
nx <- 11 #total observations

y <- metadata$d7
y <- na.exclude(y) 
mean7 <- mean(y)
y <- sum(y)
ny <- 20

z <- metadata$d14
z <- na.exclude(z)
mean14 <- mean(z)
z <- sum(z)
nz <- 16
CF <- (x + y + z) / (nx + ny + nz)

#Now we get the SS of the model
# SScondition = the sum of N * (average distance of each day - the gran mean ) ^ 2
str(metadata)
metadata$average <- as.numeric(metadata$average) # the average is from each animal, not for day
GrandMean <- mean(metadata$average)
summary_table <-metadata %>%
  summarize(
    sqrtDev0 = (mean0 - GrandMean) ^ 2, ###SS of the conditions
    sqrtDev7 = (mean7 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev14 = (mean14 - GrandMean) ^ 2,##SS of the conditions
  ) 
#write.csv(summary_table,"summary_table.csv")
##We get the SS of the conditions (time points), you multiple by (n) number of subjects under each (ith) condition
metadata
SSconditions <- 20 * sum (summary_table$sqrtDev0 + summary_table$sqrtDev7  + 
                            summary_table$sqrtDev14)
SSconditions
SSconditions
##Now we get the SS of the within groups which is each individual variarion form the group mean
# Sum of the (individual value - mean of the group) ^ 2
summary_table2 <-metadata %>%
  summarize(
    sqrtDevW0 = (d0 - mean0) ^ 2, ###SS of the conditions
    sqrtDevW7 = (d7 - mean7) ^ 2, ###SS of the conditions
    sqrtDevW14 = (d14 - mean14) ^ 2, ###SS of the conditions
  )
##removing NAs
sqrtDevW0 <- as.data.frame(summary_table2$sqrtDevW0)
sqrtDevW0 <- as.data.frame(sqrtDevW0[-c(8,10,13,14,17:20),])
colnames(sqrtDevW0) <- c("sqrtDevW0")

sqrtDevW7 <- as.data.frame(summary_table2$sqrtDevW7)
sqrtDevW7 <- as.data.frame(sqrtDevW7[-c(12),])
colnames(sqrtDevW7) <- c("sqrtDevW7")

sqrtDevW14 <- as.data.frame(summary_table2$sqrtDevW14)
sqrtDevW14 <- as.data.frame(sqrtDevW14[-c(1,2,10,20),])
colnames(sqrtDevW14) <- c("sqrtDevW14")

SSwithin <- sum (sqrtDevW0$sqrtDevW0 + sqrtDevW7$sqrtDevW7  + 
                   sqrtDevW14$sqrtDevW14)

##Now we calculate the SS for each of the subjects
## Number of variables (K) * Sum of (individual average - grand mean) ^2
#K = number of variables: 5 time point
summary_table3 <-metadata %>%
  group_by(calf) %>%
  summarize(
    sqrtDevS = (average - GrandMean) ^ 2, ###SS of the conditions
  )


## Now the we have the square value for all the different time point of each animal, we need to sum, this will give us the SS within 
SSsubject <- 3 * sum (summary_table3$sqrtDevS)

##Now we have the SSwithin and SSmodel, we need to calcualte the SSerror
SSerror <- SSwithin - SSsubject

#Now we specify the degrees of freedom of the model and for the error term
#dfmodel = k - 1, k number of variables
#dferror = (k-1) * (n-1), n is the number of subject
dfmodel = 3 - 1 #we have 3 variable or points
dfSubject = 20 -1
dferror = dfSubject * dfmodel

##Now we process the calculate the MS of the model, MSerror and Fvalue
MSconditions = SSconditions / dfmodel ##SSconditions / (k-1)
MSerror = SSerror / dferror # SSerror / (n-1)(k-1)
Fvalue = MSconditions / MSerror
results <- data.frame(MSconditions,MSerror,Fvalue)

#Calculating pvalue
#we need the degrees of freedom and the F value
p.value <- pf(Fvalue, dfmodel, dferror, lower.tail = FALSE)

#creating the data frame
Response <- c("SSconditions", "SSerror")
SS <- c(SSconditions, SSerror)
DF <- c(dfmodel,dferror)
MS <- c(MSconditions, MSerror)
F.value <- c(Fvalue, "NA")
P <- c(p.value, "NA")

df <- data.frame(Response, SS, DF, MS, F.value, P)
print (df)




##now we'll run a mixed model to find a significance
#it works but data is not normal distributed
set_sum_contrasts() # important for afex
str(distances)
levels(distances$calf) #19
# we need to add the withing factor in the error term

m1 <- mixed(distances ~ d + (d||calf), data = distances,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(m1) ## Significant 
plot(m1$full_model)
qqnorm(residuals(m1$full_model))
qqline(residuals(m1$full_model))
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_1<- emmeans(m1, "d")
emm_1
update(pairs(emm_1), by = NULL, adjust = "holm")
with(distances, shapiro.test(distances))

#now using as factor time and performing permutest with strata function
time <- betadisper(distWU, type = c("centroid"), group = healthy1$d)
time

metadata <- read.csv("betadisperCalfN.csv", header=TRUE)

time <- betadisper(distWU, type = c("centroid"), group = healthy1$d)
time
aveTime <- print.default(format(tapply(time$distances, time$group, mean))) #getting the average distance from the centroids
mean0 <- 0.4131679
mean7 <- 0.4179798
mean14 <- 0.4353913
GrandMean <- (mean0 + mean7 + mean14) / 3 

##this is a test
metadata$average <- as.numeric(metadata$average)
GrandMean2 <- mean(metadata$average) #calculating the grand mean, it has to be numeric

x <- metadata$d0
x <- na.exclude(x) 
x <- mean(x)

y <- metadata$d7
y <- na.exclude(y) 
y <- mean(y)

z <- metadata$d14
z <- na.exclude(z) 
z <- mean(z)

X <- mean(x,y,z)
GrandMean2
#upload the formatted table
metadata <- read.csv("betadisperCalfN.csv", header=TRUE)
str(metadata)
metadata2 <- na.exclude(metadata) 

#Now we get the SS of the model
# SScondition = the sum of N * (average distance of each day - the gran mean ) ^ 2
str(metadata)
metadata$average <- as.numeric(metadata$average) # the average is from each animal, not for day
summary_table <-metadata %>%
  summarize(
    sqrtDev0 = (mean0 - GrandMean) ^ 2, ###SS of the conditions
    sqrtDev7 = (mean7 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev14 = (mean14 - GrandMean) ^ 2,##SS of the conditions
  ) 
#write.csv(summary_table,"summary_table.csv")
##We get the SS of the conditions (time points), you multiple by (n) number of subjects under each (ith) condition
metadata
s1 <- 12 * summary_table$sqrtDev0
s2 <- 19 * summary_table$sqrtDev7
s3 <- 16 * summary_table$sqrtDev14
SSconditions <- sum(s1,s2,s3)
SSconditions
##Now we get the SS of the within groups which is each individual variarion form the group mean
# Sum of the (individual value - mean of the group) ^ 2
summary_table2 <-metadata %>%
  summarize(
    sqrtDevW0 = (d0 - mean0) ^ 2, ###SS of the conditions
    sqrtDevW7 = (d7 - mean7) ^ 2, ###SS of the conditions
    sqrtDevW14 = (d14 - mean14) ^ 2, ###SS of the conditions
  )
##removing NAs
sqrtDevW0 <- as.data.frame(summary_table2$sqrtDevW0)
sqrtDevW0 <- as.data.frame(sqrtDevW0[-c(8,10,13,14,17:20),])
colnames(sqrtDevW0) <- c("sqrtDevW0")

sqrtDevW7 <- as.data.frame(summary_table2$sqrtDevW7)
sqrtDevW7 <- as.data.frame(sqrtDevW7[-c(12),])
colnames(sqrtDevW7) <- c("sqrtDevW7")

sqrtDevW14 <- as.data.frame(summary_table2$sqrtDevW14)
sqrtDevW14 <- as.data.frame(sqrtDevW14[-c(1,2,10,20),])
colnames(sqrtDevW14) <- c("sqrtDevW14")

SSwithin <- sum (sqrtDevW0$sqrtDevW0 + sqrtDevW7$sqrtDevW7  + 
                   sqrtDevW14$sqrtDevW14)

##Now we calculate the SS for each of the subjects
## Number of variables (K) * Sum of (individual average - grand mean) ^2
#K = number of variables: 5 time point
summary_table3 <-metadata %>%
  group_by(calf) %>%
  summarize(
    sqrtDevS = (average - GrandMean) ^ 2, ###SS of the conditions
  )


## Now the we have the square value for all the different time point of each animal, we need to sum, this will give us the SS within 
SSsubject <- 3 * sum (summary_table3$sqrtDevS)

##Now we have the SSwithin and SSmodel, we need to calcualte the SSerror
SSerror <- SSwithin - SSsubject

#Now we specify the degrees of freedom of the model and for the error term
#dfmodel = k - 1, k number of variables
#dferror = (k-1) * (n-1), n is the number of subject
dfmodel = 3 - 1 #we have 3 variable or points
dfSubject = 20 -1
dferror = dfSubject * dfmodel

##Now we process the calculate the MS of the model, MSerror and Fvalue
MSconditions = SSconditions / dfmodel ##SSconditions / (k-1)
MSerror = SSerror / dferror # SSerror / (n-1)(k-1)
Fvalue = MSconditions / MSerror
results <- data.frame(MSconditions,MSerror,Fvalue)

#Calculating pvalue
#we need the degrees of freedom and the F value
p.value <- pf(Fvalue, dfmodel, dferror, lower.tail = FALSE)

#creating the data frame
Response <- c("SSconditions", "SSerror")
SS <- c(SSconditions, SSerror)
DF <- c(dfmodel,dferror)
MS <- c(MSconditions, MSerror)
F.value <- c(Fvalue, "NA")
P <- c(p.value, "NA")

df <- data.frame(Response, SS, DF, MS, F.value, P)
print (df)

### ----------------------------------Performing a Tuker-Kramer pairwise comparisons
#we need the mean of each condition (timepoints)
#Sample size in each condition (N)
# Numerator and denominator degrees of freedom
# Q value (using the Standarized range Q table for the alpha used 0.05)
# Critical Value
# MSE of the model

healthy1 %>% tally()
healthy1 %>% count(day)
#0 12
#7 16
#14 19
#sample size is uneven 

#we already have the distance means from each group
mean0
mean7
mean14

#Let's specify the total number of combinations: k(k-1) / 2 
( 3 * (3-1)) / 2 #3 total combination between the different time points
#Make the table
d0d7= abs(mean0-mean7)
d0d14= abs(mean0-mean14)
d7d14= abs(mean7-mean14)

Comp <- c("d0-d7","d0-d14","d7-d14")
Abs <- c(d0d7, d0d14,d7d14)
data <- data.frame(Comp,Abs)
data
##Now we calculate the critical value 
#we need the degrees of freedom
#Numerator: (3 -1)= 2 
#denominator (k-1)(n-1) = (3-1)(20-1) = 38
Q = 2.863
# Q value: there is no 236 df in the table, so we'll use the closet one: 120, with an alpha of 0.05
healthy1 %>% count(day)
#0 12
#7 16
#14 19
#sample size is uneven 
Abs1 <- Q * sqrt((MSerror / 2) * ((1 / 12) + (1/ 16))) #0-7
Abs2 <- Q * sqrt((MSerror / 2) * ((1 / 12) + (1/ 19))) #0-14
Abs3 <- Q * sqrt((MSerror / 2) * ((1 / 16) + (1/ 19))) # they have the similar number of sample size 7-14

CD <- c(Abs1, Abs2,Abs3)
data <- data.frame(Comp,Abs, CD)
str(data)

data$significant <- ifelse(
  data$Abs > data$CD, 'Sig',
  ifelse(data$Abs < data$CD, 'NS',
         'Mid Change'))
data


### --------------------------------------------Bray curtis distances 
#We need to import the distance matrices
DistWU = read.csv("dist.brayH.csv", header=TRUE)
row.names(DistWU) <- DistWU[,1]
DistWU = DistWU[,-c(1)]
distWU <- as.dist(DistWU)
#check homogeneity of variance

timeD <- betadisper(distWU, type = c("centroid"), group = healthy1$calf)
timeD
distances <- timeD[["distances"]]
distances <- as.data.frame(distances)
distances <- merge(distances, healthy1, by.x = 0, by.y=0) # we are combining the data and removing unneccesary data
distances <- distances[,-c(1, 5, 7:9)]

average <- print.default(format(tapply(timeD$distances, timeD$group, mean))) #getting the average distance from the centroids
average <- as.data.frame(average)
average$calf <- row.names(average)
#write.csv(average,"average.csv")
distances <- merge(distances, average, by.x = "calf", by.y="calf")

#arraging the metadata
head(distances)
metadata <- distances[,c(1,2,6,5,3)]
#write.csv(metadata,"metadata2.csv") 

time <- betadisper(distWU, type = c("centroid"), group = healthy1$d)
time
aveTime <- print.default(format(tapply(time$distances, time$group, mean))) #getting the average distance from the centroids
mean0 <- 0.5964075
mean7 <- 0.6152064
mean14 <-0.6144381

#upload the formatted table
metadata <- read.csv("betadisperCalfWUH.csv", header=TRUE)
str(metadata)

#Now we get the SS of the model
# SScondition = the sum of N * (average distance of each day - the gran mean ) ^ 2
str(metadata)
metadata$average <- as.numeric(metadata$average)
GrandMean <- mean(metadata$average) #calculating the grand mean, it has to be numeric
summary_table <-metadata %>%
  summarize(
    sqrtDev0 = (mean0 - GrandMean) ^ 2, ###SS of the conditions
    sqrtDev7 = (mean7 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev14 = (mean14 - GrandMean) ^ 2,##SS of the conditions
  ) 
SSconditions <- 20 * sum (summary_table$sqrtDev0 + summary_table$sqrtDev7  + 
                            summary_table$sqrtDev14)
SSconditions
##Now we get the SS of the within groups which is each individual variarion form the group mean
# Sum of the (individual value - mean of the group) ^ 2
summary_table2 <-metadata %>%
  summarize(
    sqrtDevW0 = (d0 - mean0) ^ 2, ###SS of the conditions
    sqrtDevW7 = (d7 - mean7) ^ 2, ###SS of the conditions
    sqrtDevW14 = (d14 - mean14) ^ 2, ###SS of the conditions
  )
##removing NAs
sqrtDevW0 <- as.data.frame(summary_table2$sqrtDevW0)
sqrtDevW0 <- as.data.frame(sqrtDevW0[-c(8,10,13,14,17:20),])
colnames(sqrtDevW0) <- c("sqrtDevW0")

sqrtDevW7 <- as.data.frame(summary_table2$sqrtDevW7)
sqrtDevW7 <- as.data.frame(sqrtDevW7[-c(12),])
colnames(sqrtDevW7) <- c("sqrtDevW7")

sqrtDevW14 <- as.data.frame(summary_table2$sqrtDevW14)
sqrtDevW14 <- as.data.frame(sqrtDevW14[-c(1,2,10,20),])
colnames(sqrtDevW14) <- c("sqrtDevW14")

SSwithin <- sum (sqrtDevW0$sqrtDevW0 + sqrtDevW7$sqrtDevW7  + 
                   sqrtDevW14$sqrtDevW14)
##Now we calculate the SS for each of the subjects
## Number of variables (K) * Sum of (individual average - grand mean) ^2
#K = number of variables: 5 time point
summary_table3 <-metadata %>%
  group_by(calf) %>%
  summarize(
    sqrtDevS = (average - GrandMean) ^ 2, ###SS of the conditions
  )


## Now the we have the square value for all the different time point of each animal, we need to sum, this will give us the SS within 
SSsubject <- 3 * sum (summary_table3$sqrtDevS)

##Now we have the SSwithin and SSmodel, we need to calcualte the SSerror
SSerror <- SSwithin - SSsubject

#Now we specify the degrees of freedom of the model and for the error term
#dfmodel = k - 1, k number of variables
#dferror = (k-1) * (n-1), n is the number of subject
dfmodel = 3 - 1 #we have 3 variable or points
dfSubject = 20 -1
dferror = dfSubject * dfmodel

##Now we process the calculate the MS of the model, MSerror and Fvalue
MSconditions = SSconditions / dfmodel ##SSconditions / (k-1)
MSerror = SSerror / dferror # SSerror / (n-1)(k-1)
Fvalue = MSconditions / MSerror
resultsBC <- data.frame(MSconditions,MSerror,Fvalue)

#Calculating pvalue
#we need the degrees of freedom and the F value
p.value <- pf(Fvalue, dfmodel, dferror, lower.tail = FALSE)

#creating the data frame
Response <- c("SSconditions", "SSerror")
SS <- c(SSconditions, SSerror)
DF <- c(dfmodel,dferror)
MS <- c(MSconditions, MSerror)
F.value <- c(Fvalue, "NA")
P <- c(p.value, "NA")

df2 <- data.frame(Response, SS, DF, MS, F.value, P)
print (df2)

### ----------------------------------Performing a Tuker-Kramer pairwise comparisons
#we already have the distance means from each group ~bray curtis
mean0
mean7
mean14

#Let's specify the total number of combinations: k(k-1) / 2 
( 3 * (3-1)) / 2 #3 total combination between the different time points
#Make the table
d0d7= abs(mean0-mean7)
d0d14= abs(mean0-mean14)
d7d14= abs(mean7-mean14)

Comp <- c("d0-d7","d0-d14","d7-d14")
Abs <- c(d0d7, d0d14,d7d14)
data2 <- data.frame(Comp,Abs)
data2
##Now we calculate the critical value 
#we need the degrees of freedom
#Numerator: (3 -1)= 2 
#denominator (k-1)(n-1) = (3-1)(20-1) = 38
Q = 2.863
# Q value: there is no 236 df in the table, so we'll use the closet one: 120, with an alpha of 0.05
healthy1 %>% count(day)
#0 12
#7 16
#14 19
#sample size is uneven, we will use the same Critical values calculated before 
#Abs1, Abs2, Abs3 

CD <- c(Abs1, Abs2,Abs3)
data2 <- data.frame(Comp,Abs, CD)
str(data2)

data2$significance <- ifelse(
  data2$Abs > data2$CD, 'Sig',
  ifelse(data2$Abs < data2$CD, 'NS',
         'Mid Change'))
data2

###------------------------------------- Sick animals
ant <- read.table("CalfMetadataSickAnimals.txt", header=TRUE, sep="\t")
ant <- ant[,-c(11:14)]
#also using d as continuous factoe
rownames(ant) <- ant$ID2
ant$antibiotic <- as.factor(ant$antibiotic)
ant$sample <- ant$ID2
ant$day <- factor(ant$day, levels = c("s0","s1","s5", "s10"))
levels(ant$day ) <- list("0"="s0", "1"="s1", "5"="s5", "10"="s10")
X <- as.factor(ant$calf)


#We need to import the distance matrices
DistWU = read.csv("WU.S.csv", header=TRUE)
row.names(DistWU) <- DistWU[,1]
DistWU = DistWU[,-c(1)]
distWU <- as.dist(DistWU)
#check homogeneity of variance
ant <-ant[rownames(ant) %in% rownames(DistWU),]#32


timeD1 <- betadisper(distWU, type = c("centroid"), group = ant$calf)
timeD1
distances <- timeD1[["distances"]]
distances <- as.data.frame(distances)
distances <- merge(distances, ant, by.x = 0, by.y=0) # we are combining the data and removing unneccesary data
distances <- distances[,-c(1, 4,5, 7:9,11:13)]

average <- print.default(format(tapply(timeD1$distances, timeD1$group, mean))) #getting the average distance from the centroids
average <- as.data.frame(average)
average$calf <- row.names(average)
#write.csv(average,"average.csv")
distances <- merge(distances, average, by.x = "calf", by.y="calf")

#arraging the metadata
head(distances)
metadata2 <- distances[,c(1,2,5,3,4)]
head(metadata2)
#write.csv(metadata2,"metadata2S.csv") 

time <- betadisper(distWU, type = c("centroid"), group = ant$day)
time
aveTime <- print.default(format(tapply(time$distances, time$group, mean))) #getting the average distance from the centroids
mean0 <-0.4140507
mean1 <- 0.3178659
mean5 <- 0.3585800
mean10 <- 0.3871404

#upload the formatted table
metadata <- read.csv("betadisperCalfWUS.csv", header=TRUE)
str(metadata)

#Now we get the SS of the model
# SScondition = the sum of N * (average distance of each day - the gran mean ) ^ 2
str(metadata)
metadata$average <- as.numeric(metadata$average)
GrandMean <- mean(metadata$average) #calculating the grand mean, it has to be numeric
summary_table <-metadata %>%
  summarize(
    sqrtDev0 = (mean0 - GrandMean) ^ 2, ###SS of the conditions
    sqrtDev1 = (mean1 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev5 = (mean5 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev10 = (mean10 - GrandMean) ^ 2,##SS of the conditions
  ) 
#write.csv(summary_table,"summary_table.csv")
##We get the SS of the conditions (time points)
SSconditions <- 10 * sum (summary_table$sqrtDev0 + summary_table$sqrtDev1  + 
                            summary_table$sqrtDev5 + summary_table$sqrtDev10)
SSconditions
##Now we get the SS of the within groups which is each individual variarion form the group mean
# Sum of the (individual value - mean of the group) ^ 2
summary_table2 <-metadata %>%
  summarize(
    sqrtDevW0 = (s0 - mean0) ^ 2, ###SS of the conditions
    sqrtDevW1 = (s1 - mean1) ^ 2, ###SS of the conditions
    sqrtDevW5 = (s5 - mean5) ^ 2,
    sqrtDevW10 = (s10 - mean10) ^ 2,###SS of the conditions
  )
##removing NAs
sqrtDevW0 <- as.data.frame(summary_table2$sqrtDevW0)
sqrtDevW0 <- as.data.frame(sqrtDevW0[-c(3,7),])
colnames(sqrtDevW0) <- c("sqrtDevW0")

sqrtDevW1 <- as.data.frame(summary_table2$sqrtDevW1)
sqrtDevW1 <- as.data.frame(sqrtDevW1[-c(1,7),])
colnames(sqrtDevW1) <- c("sqrtDevW1")

sqrtDevW5 <- as.data.frame(summary_table2$sqrtDevW5)
sqrtDevW5 <- as.data.frame(sqrtDevW5[-c(1,2),])
colnames(sqrtDevW5) <- c("sqrtDevW5")

sqrtDevW10 <- as.data.frame(summary_table2$sqrtDevW10)
sqrtDevW10 <- as.data.frame(sqrtDevW10[-c(1),])
colnames(sqrtDevW10) <- c("sqrtDevW10")

SSwithin <- sum (sqrtDevW0$sqrtDevW0 + sqrtDevW1$sqrtDevW1  + sqrtDevW5$sqrtDevW5 +
                   sqrtDevW10$sqrtDevW10)

##Now we calculate the SS for each of the subjects
## Number of variables (K) * Sum of (individual average - grand mean) ^2
#K = number of variables: 5 time point
summary_table3 <-metadata %>%
  group_by(calf) %>%
  summarize(
    sqrtDevS = (average - GrandMean) ^ 2, ###SS of the conditions
  )


## Now the we have the square value for all the different time point of each animal, we need to sum, this will give us the SS within 
SSsubject <- 4 * sum (summary_table3$sqrtDevS)

##Now we have the SSwithin and SSmodel, we need to calcualte the SSerror
SSerror <- SSwithin - SSsubject

#Now we specify the degrees of freedom of the model and for the error term
#dfmodel = k - 1, k number of variables
#dferror = (k-1) * (n-1), n is the number of subject
dfmodel = 4 - 1 #we have 3 variable or points
dfSubject = 10 -1
dferror = dfSubject * dfmodel

##Now we process the calculate the MS of the model, MSerror and Fvalue
MSconditions = SSconditions / dfmodel ##SSconditions / (k-1)
MSerror = SSerror / dferror # SSerror / (n-1)(k-1)
Fvalue = MSconditions / MSerror
results <- data.frame(MSconditions,MSerror,Fvalue)

#Calculating pvalue
#we need the degrees of freedom and the F value
p.value <- pf(Fvalue, dfmodel, dferror, lower.tail = FALSE)

#creating the data frame
Response <- c("SSconditions", "SSerror")
SS <- c(SSconditions, SSerror)
DF <- c(dfmodel,dferror)
MS <- c(MSconditions, MSerror)
F.value <- c(Fvalue, "NA")
P <- c(p.value, "NA")

df3 <- data.frame(Response, SS, DF, MS, F.value, P)
print (df3)

###------------------ Performing a Tuker-Kramer pairwise comparisons
#we already have the distance means from each group ~bray curtis
mean0
mean1
mean5
mean10

#Let's specify the total number of combinations: k(k-1) / 2 
( 4 * (4-1)) / 2 #6 total combination between the different time points
#Make the table
d0d1= abs(mean0-mean1)
d0d5= abs(mean0-mean5)
d0d10= abs(mean0-mean10)
d1d5= abs(mean1-mean5)
d1d10= abs(mean1-mean10)
d5d10= abs(mean5-mean10)

Comp <- c("d0-d1","d0-d5","d0-d10","d1-d5","d1-d10","d5-d10")
Abs <- c(d0d1, d0d5,d0d10,d1d5, d1d10,d5d10)
data2 <- data.frame(Comp,Abs)
data2
##Now we calculate the critical value 
#we need the degrees of freedom
#Numerator: (4 -1)= 3
#denominator (k-1)(n-1) = (4-1)(10-1) = 27
Q = 3.506
# Q value: there is no 236 df in the table, so we'll use the closet one: 120, with an alpha of 0.05
ant %>% count(day)
#0 8
#1 8
#5 8
#10 9
#sample size is uneven, we will use the same Critical values calculated before 
#Abs1, Abs2, Abs3 
#sample size is uneven 
Abs1 <- Q * sqrt(MSerror / 8)  #0-1
Abs2 <- Q * sqrt(MSerror / 8) # they have the similar number of sample size 0-5
Abs3 <- Q * sqrt((MSerror / 2) * ((1 / 8) + (1/ 9))) #0-10
Abs4 <- Q * sqrt(MSerror / 8)  # they have the similar number of sample size 1-5
Abs5 <- Q * sqrt((MSerror / 2) * ((1 / 8) + (1/ 9)))  # they have the similar number of sample size 1-10
Abs6 <- Q * sqrt((MSerror / 2) * ((1 / 8) + (1/ 9))) # they have the similar number of sample size 5-10


CD <- c(Abs1, Abs2,Abs3,Abs4, Abs5,Abs6)
data2 <- data.frame(Comp,Abs, CD)
str(data2)

data2$significance <- ifelse(
  data2$Abs > data2$CD, 'Sig',
  ifelse(data2$Abs < data2$CD, 'NS',
         'Mid Change'))
data2

#--------We need to import the distance matrices
DistBC = read.csv("dist.brayS.csv", header=TRUE)
row.names(DistBC) <- DistBC[,1]
DistBC = DistBC[,-c(1)]
distBC <- as.dist(DistBC)
#check homogeneity of variance

timeD2 <- betadisper(distBC, type = c("centroid"), group = ant$calf)
timeD2
distances <- timeD2[["distances"]]
distances <- as.data.frame(distances)
distances <- merge(distances, ant, by.x = 0, by.y=0) # we are combining the data and removing unneccesary data
distances <- distances[,-c(1, 4,5, 7:12)]

average <- print.default(format(tapply(timeD2$distances, timeD2$group, mean))) #getting the average distance from the centroids
average <- as.data.frame(average)
average$calf <- row.names(average)
distances <- merge(distances, average, by.x = "calf", by.y="calf")

#arraging the metadata
head(distances)
metadata2 <- distances[,c(1,2,5,3,4)]
#write.csv(metadata2,"metadataS2.csv") 

time <- betadisper(distBC, type = c("centroid"), group = ant$day)
time
aveTime <- print.default(format(tapply(time$distances, time$group, mean))) #getting the average distance from the centroids
mean0 <-0.6101746
mean1 <- 0.4904151
mean5 <- 0.5904921
mean10 <- 0.6018749

#upload the formatted table
metadata <- read.csv("betadisperCalfBCS.csv", header=TRUE)
str(metadata)

#Now we get the SS of the model
# SScondition = the sum of N * (average distance of each day - the gran mean ) ^ 2
str(metadata)
metadata$average <- as.numeric(metadata$average)
GrandMean <- mean(metadata$average) #calculating the grand mean, it has to be numeric
summary_table <-metadata %>%
  summarize(
    sqrtDev0 = (mean0 - GrandMean) ^ 2, ###SS of the conditions
    sqrtDev1 = (mean1 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev5 = (mean5 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev10 = (mean10 - GrandMean) ^ 2,##SS of the conditions
  ) 
#write.csv(summary_table,"summary_table.csv")
##We get the SS of the conditions (time points)
SSconditions <- 10 * sum (summary_table$sqrtDev0 + summary_table$sqrtDev1  + 
                           summary_table$sqrtDev5 + summary_table$sqrtDev10)
SSconditions
##Now we get the SS of the within groups which is each individual variarion form the group mean
# Sum of the (individual value - mean of the group) ^ 2
summary_table2 <-metadata %>%
  summarize(
    sqrtDevW0 = (s0 - mean0) ^ 2, ###SS of the conditions
    sqrtDevW1 = (s1 - mean1) ^ 2, ###SS of the conditions
    sqrtDevW5 = (s5 - mean5) ^ 2, ###SS of the conditions
    sqrtDevW10 = (s10 - mean10) ^ 2, ###SS of the conditions
  )

##removing NAs
sqrtDevW0 <- as.data.frame(summary_table2$sqrtDevW0)
sqrtDevW0 <- as.data.frame(sqrtDevW0[-c(3,7),])
colnames(sqrtDevW0) <- c("sqrtDevW0")

sqrtDevW1 <- as.data.frame(summary_table2$sqrtDevW1)
sqrtDevW1 <- as.data.frame(sqrtDevW1[-c(1,7),])
colnames(sqrtDevW1) <- c("sqrtDevW1")

sqrtDevW5 <- as.data.frame(summary_table2$sqrtDevW5)
sqrtDevW5 <- as.data.frame(sqrtDevW5[-c(1,2),])
colnames(sqrtDevW5) <- c("sqrtDevW5")

sqrtDevW10 <- as.data.frame(summary_table2$sqrtDevW10)
sqrtDevW10 <- as.data.frame(sqrtDevW10[-c(1),])
colnames(sqrtDevW10) <- c("sqrtDevW10")

SSwithin <- sum (sqrtDevW0$sqrtDevW0 + sqrtDevW1$sqrtDevW1  + sqrtDevW5$sqrtDevW5 +
                   sqrtDevW10$sqrtDevW10)

##Now we calculate the SS for each of the subjects
## Number of variables (K) * Sum of (individual average - grand mean) ^2
#K = number of variables: 5 time point
summary_table3 <-metadata %>%
  group_by(calf) %>%
  summarize(
    sqrtDevS = (average - GrandMean) ^ 2, ###SS of the conditions
  )


## Now the we have the square value for all the different time point of each animal, we need to sum, this will give us the SS within 
SSsubject <- 4 * sum (summary_table3$sqrtDevS)

##Now we have the SSwithin and SSmodel, we need to calcualte the SSerror
SSerror <- SSwithin - SSsubject

#Now we specify the degrees of freedom of the model and for the error term
#dfmodel = k - 1, k number of variables
#dferror = (k-1) * (n-1), n is the number of subject
dfmodel = 4 - 1 #we have 3 variable or points
dfSubject = 10 -1
dferror = dfSubject * dfmodel

##Now we process the calculate the MS of the model, MSerror and Fvalue
MSconditions = SSconditions / dfmodel ##SSconditions / (k-1)
MSerror = SSerror / dferror # SSerror / (n-1)(k-1)
Fvalue = MSconditions / MSerror
results <- data.frame(MSconditions,MSerror,Fvalue)

#Calculating pvalue
#we need the degrees of freedom and the F value
p.value <- pf(Fvalue, dfmodel, dferror, lower.tail = FALSE)

#creating the data frame
Response <- c("SSconditions", "SSerror")
SS <- c(SSconditions, SSerror)
DF <- c(dfmodel,dferror)
MS <- c(MSconditions, MSerror)
F.value <- c(Fvalue, "NA")
P <- c(p.value, "NA")

df4 <- data.frame(Response, SS, DF, MS, F.value, P)
print (df4)

#PAIRWISE
#we already have the distance means from each group ~bray curtis
mean0
mean1
mean5
mean10

#Let's specify the total number of combinations: k(k-1) / 2 
( 4 * (4-1)) / 2 #6 total combination between the different time points
#Make the table
d0d1= abs(mean0-mean1)
d0d5= abs(mean0-mean5)
d0d10= abs(mean0-mean10)
d1d5= abs(mean1-mean5)
d1d10= abs(mean1-mean10)
d5d10= abs(mean5-mean10)

Comp <- c("d0-d1","d0-d5","d0-d10","d1-d5","d1-d10","d5-d10")
Abs <- c(d0d1, d0d5,d0d10,d1d5, d1d10,d5d10)
data2 <- data.frame(Comp,Abs)
data2
##Now we calculate the critical value 
#we need the degrees of freedom
#Numerator: (4 -1)= 3
#denominator (k-1)(n-1) = (4-1)(9-1) = 24
Q = 3.506
# Q value: there is no 236 df in the table, so we'll use the closet one: 120, with an alpha of 0.05
ant %>% count(day)
#0 8
#1 8
#5 8
#10 9
#sample size is uneven, we will use the same Critical values calculated before 
#Abs1, Abs2, Abs3 
#sample size is uneven 
Abs1 <- Q * sqrt(MSerror / 8) #0-1
Abs2 <- Q * sqrt(MSerror / 8) # they have the similar number of sample size 0-5
Abs3 <- Q * sqrt((MSerror / 2) * ((1 / 8) + (1/ 9))) #0-10
Abs4 <- Q * sqrt(MSerror / 8) # they have the similar number of sample size 1-5
Abs5 <- Q * sqrt((MSerror / 2) * ((1 / 8) + (1/ 9))) # they have the similar number of sample size 1-10
Abs6 <- Q * sqrt((MSerror / 2) * ((1 / 8) + (1/ 9))) # they have the similar number of sample size 5-10


CD <- c(Abs1, Abs2,Abs3,Abs4, Abs5,Abs6)
data2 <- data.frame(Comp,Abs, CD)
str(data2)

data2$significance <- ifelse(
  data2$Abs > data2$CD, 'Sig',
  ifelse(data2$Abs < data2$CD, 'NS',
         'Mid Change'))
data2
