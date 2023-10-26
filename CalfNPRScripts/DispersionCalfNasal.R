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
#write.csv(distances, "distancesHWU.csv")

#---Step 1: Within SS
#we need to calculate the group mean and then the grand mena
metadata <- read.csv("betadisperCalfN.csv", header=TRUE)
x <- metadata$d0
x <- na.exclude(x) 
mean0 <- mean(x)

y <- metadata$d7
y <- na.exclude(y) 
mean7 <- mean(y)

z <- metadata$d14
z <- na.exclude(z)
mean14 <- mean(z)
GrandMean <- (mean0 + mean7 + mean14) / 3

summary_table <-metadata %>%
  summarize(
    sqrtDevW0 = (d0 - mean0) ^ 2, ###SS of the conditions
    sqrtDevW7 = (d7 - mean7) ^ 2, ###SS of the conditions
    sqrtDevW14 = (d14 - mean14) ^ 2, ###SS of the conditions
  )
##removing NAs
sqrtDevW0 <- as.data.frame(summary_table$sqrtDevW0)
sqrtDevW0 <- as.data.frame(sqrtDevW0[-c(8,10,13,14,17:20),])
colnames(sqrtDevW0) <- c("sqrtDevW0")

sqrtDevW7 <- as.data.frame(summary_table$sqrtDevW7)
sqrtDevW7 <- as.data.frame(sqrtDevW7[-c(12),])
colnames(sqrtDevW7) <- c("sqrtDevW7")

sqrtDevW14 <- as.data.frame(summary_table$sqrtDevW14)
sqrtDevW14 <- as.data.frame(sqrtDevW14[-c(1,2,10,20),])
colnames(sqrtDevW14) <- c("sqrtDevW14")

SSwithin <- sum (sqrtDevW0$sqrtDevW0 + sqrtDevW7$sqrtDevW7  + 
                   sqrtDevW14$sqrtDevW14)

### Step 2: SS between (or conditions)
# the group sample size multiplied by the square different between the group's mean and grand mean
metadata
summary_table2 <-metadata %>%
  summarize(
    sqrtDev0 = (mean0 - GrandMean) ^ 2, ###SS of the conditions
    sqrtDev7 = (mean7 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev14 = (mean14 - GrandMean) ^ 2,##SS of the conditions
  ) 
#write.csv(summary_table,"summary_table.csv")
##We get the SS of the conditions (time points), you multiple by (n) number of subjects under each (ith) condition
metadata
s1 <- 12 * summary_table2$sqrtDev0
s2 <- 19 * summary_table2$sqrtDev7
s3 <- 16 * summary_table2$sqrtDev14
SSconditions <- sum(s1,s2,s3)
SSconditions

## Now, let's calculate the SStotal= SSbetween + SS within
SStotal <- SSwithin + SSconditions

## Step 3-- SSsubjects
#we need to calculate the mean of each subject across the 3 timepoints and then you subtract from the grand mean and you multiply by the total number of samples
summary_table3 <-metadata %>%
  group_by(calf, observations) %>%
  summarize(
    sqrtDevS = (average - GrandMean) ^ 2, ###SS of the conditions
  )
#now we multiply that number by the total number of observations per animal
SSsubject <- summary_table3 %>%
  group_by(calf) %>%
  summarize(
    sqrtDevS2 = sqrtDevS * observations, ###SS of the conditions, the animals have different number of observations
  ) 
SSsubject <- sum(SSsubject$sqrtDevS2)

##Step 4- SSErroc
SSerror <- SSwithin - SSsubject
#now the double check
SSwithin
total <- SSsubject + SSerror

###now we move into calculating the degrees of freedom
#Now we specify the degrees of freedom of the model and for the error term
#dfmodel = k - 1, k number of variables
#dferror = (k-1) * (n-1), n is the number of subject
dfbetween = 3 - 1 #we have 3 variable or points
dfwithin = 47 - 3 #(total number of data points -number of groups)
dfsubject = 20 - 1 #(total number of subjects -1)
dferror = dfwithin - dfsubject
dftotal = 47 - 1 #total number of data points - 1

##now we calculate the Mean Square (SS for each df)
MSconditions = SSconditions / dfbetween
MSwithin = SSwithin / dfwithin
MSsubject = SSsubject / dfsubject
MSerror = SSerror / dferror

### now we calculate the F-ratio
Fratio = MSconditions / MSerror
# we then use the F-distribution with the degrees of freedom of the between and error
# now we calculate the P.value
p.value <- pf(Fratio, dfbetween, dferror, lower.tail = FALSE)

#creating the data frame
Response <- c("SSconditions", "SSwithin", "SSsubject", "SSerror", "SStotal")
SS <- c(SSconditions, SSwithin, SSsubject, SSerror, SStotal)
DF <- c(dfbetween, dfwithin, dfsubject, dferror, "NA")
MS <- c(MSconditions, MSwithin, MSsubject, MSerror, "NA")
F.value <- c(Fratio, "NA","NA","NA","NA")
P <- c(p.value, "NA","NA","NA","NA")

df1 <- data.frame(Response, SS, DF, MS, F.value, P)
print (df1)

##----------Bray curtis
#We need to import the distance matrices
DistWU = read.csv("dist.brayH.csv", header=TRUE)
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
#write.csv(distances, "distancesHBC.csv")

#---Step 1: Within SS
#we need to calculate the group mean and then the grand mena
metadata <- read.csv("betadisperCalfBCH.csv", header=TRUE)
x <- metadata$d0
x <- na.exclude(x) 
mean0 <- mean(x)

y <- metadata$d7
y <- na.exclude(y) 
mean7 <- mean(y)

z <- metadata$d14
z <- na.exclude(z)
mean14 <- mean(z)
GrandMean <- (mean0 + mean7 + mean14) / 3

summary_table <-metadata %>%
  summarize(
    sqrtDevW0 = (d0 - mean0) ^ 2, ###SS of the conditions
    sqrtDevW7 = (d7 - mean7) ^ 2, ###SS of the conditions
    sqrtDevW14 = (d14 - mean14) ^ 2, ###SS of the conditions
  )
##removing NAs
sqrtDevW0 <- as.data.frame(summary_table$sqrtDevW0)
sqrtDevW0 <- as.data.frame(sqrtDevW0[-c(8,10,13,14,17:20),])
colnames(sqrtDevW0) <- c("sqrtDevW0")

sqrtDevW7 <- as.data.frame(summary_table$sqrtDevW7)
sqrtDevW7 <- as.data.frame(sqrtDevW7[-c(12),])
colnames(sqrtDevW7) <- c("sqrtDevW7")

sqrtDevW14 <- as.data.frame(summary_table$sqrtDevW14)
sqrtDevW14 <- as.data.frame(sqrtDevW14[-c(1,2,10,20),])
colnames(sqrtDevW14) <- c("sqrtDevW14")

SSwithin <- sum (sqrtDevW0$sqrtDevW0 + sqrtDevW7$sqrtDevW7  + 
                   sqrtDevW14$sqrtDevW14)

### Step 2: SS between (or conditions)
# the group sample size multiplied by the square different between the group's mean and grand mean
metadata
summary_table2 <-metadata %>%
  summarize(
    sqrtDev0 = (mean0 - GrandMean) ^ 2, ###SS of the conditions
    sqrtDev7 = (mean7 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev14 = (mean14 - GrandMean) ^ 2,##SS of the conditions
  ) 
#write.csv(summary_table,"summary_table.csv")
##We get the SS of the conditions (time points), you multiple by (n) number of subjects under each (ith) condition
metadata
s1 <- 12 * summary_table2$sqrtDev0
s2 <- 19 * summary_table2$sqrtDev7
s3 <- 16 * summary_table2$sqrtDev14
SSconditions <- sum(s1,s2,s3)
SSconditions

## Now, let's calculate the SStotal= SSbetween + SS within
SStotal <- SSwithin + SSconditions

## Step 3-- SSsubjects
#we need to calculate the mean of each subject across the 3 timepoints and then you subtract from the grand mean and you multiply by the total number of samples
summary_table3 <-metadata %>%
  group_by(calf, observations) %>%
  summarize(
    sqrtDevS = (average - GrandMean) ^ 2, ###SS of the conditions
  )
#now we multiply that number by the total number of observations per animal
SSsubject <- summary_table3 %>%
  group_by(calf) %>%
  summarize(
    sqrtDevS2 = sqrtDevS * observations, ###SS of the conditions, the animals have different number of observations
  ) 
SSsubject <- sum(SSsubject$sqrtDevS2)

##Step 4- SSErroc
SSerror <- SSwithin - SSsubject
#now the double check
SSwithin
total <- SSsubject + SSerror

###now we move into calculating the degrees of freedom
#Now we specify the degrees of freedom of the model and for the error term
#dfmodel = k - 1, k number of variables
#dferror = (k-1) * (n-1), n is the number of subject
dfbetween = 3 - 1 #we have 3 variable or points
dfwithin = 47 - 3 #(total number of data points -number of groups)
dfsubject = 20 - 1 #(total number of subjects -1)
dferror = dfwithin - dfsubject
dftotal = 47 - 1 #total number of data points - 1

##now we calculate the Mean Square (SS for each df)
MSconditions = SSconditions / dfbetween
MSwithin = SSwithin / dfwithin
MSsubject = SSsubject / dfsubject
MSerror = SSerror / dferror

### now we calculate the F-ratio
Fratio = MSconditions / MSerror
# we then use the F-distribution with the degrees of freedom of the between and error
# now we calculate the P.value
p.value <- pf(Fratio, dfbetween, dferror, lower.tail = FALSE)

#creating the data frame
Response <- c("SSconditions", "SSwithin", "SSsubject", "SSerror", "SStotal")
SS <- c(SSconditions, SSwithin, SSsubject, SSerror, SStotal)
DF <- c(dfbetween, dfwithin, dfsubject, dferror, "NA")
MS <- c(MSconditions, MSwithin, MSsubject, MSerror, "NA")
F.value <- c(Fratio, "NA","NA","NA","NA")
P <- c(p.value, "NA","NA","NA","NA")

df2 <- data.frame(Response, SS, DF, MS, F.value, P)
print (df2) #bray curtis

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


timeD1 <- betadisper(distWU, type = c("centroid"), group = ant$d)
timeD1
distances <- timeD1[["distances"]]
distances <- as.data.frame(distances)
distances <- merge(distances, ant, by.x = 0, by.y=0) # we are combining the data and removing unneccesary data
distances <- distances[,-c(1, 4,5, 7:9,11:13)]
#write.csv(distances, "distancesWUS.csv")

#upload the formatted table
metadata <- read.csv("betadisperCalfWUS.csv", header=TRUE)
str(metadata)
x <- metadata$s0
x <- na.exclude(x) 
mean0 <- mean(x)

y <- metadata$s1
y <- na.exclude(y) 
mean1 <- mean(y)

z <- metadata$s5
z <- na.exclude(z)
mean5 <- mean(z)

w <- metadata$s10
w <- na.exclude(w)
mean10 <- mean(w)

GrandMean <- (mean0 + mean1 + mean5 + mean10) / 4

summary_table <-metadata %>%
  summarize(
    sqrtDevW0 = (s0 - mean0) ^ 2, ###SS of the conditions
    sqrtDevW1 = (s1 - mean1) ^ 2, ###SS of the conditions
    sqrtDevW5 = (s5 - mean5) ^ 2,
    sqrtDevW10 = (s10 - mean10) ^ 2, ###SS of the conditions
  )
##removing NAs
sqrtDevW0 <- as.data.frame(summary_table$sqrtDevW0)
sqrtDevW0 <- as.data.frame(sqrtDevW0[-c(3,7),])
colnames(sqrtDevW0) <- c("sqrtDevW0")

sqrtDevW1 <- as.data.frame(summary_table$sqrtDevW1)
sqrtDevW1 <- as.data.frame(sqrtDevW1[-c(1,7),])
colnames(sqrtDevW1) <- c("sqrtDevW1")

sqrtDevW5 <- as.data.frame(summary_table$sqrtDevW5)
sqrtDevW5 <- as.data.frame(sqrtDevW5[-c(1,2),])
colnames(sqrtDevW5) <- c("sqrtDevW5")

sqrtDevW10 <- as.data.frame(summary_table$sqrtDevW10)
sqrtDevW10 <- as.data.frame(sqrtDevW10[-c(1),])
colnames(sqrtDevW10) <- c("sqrtDevW10")

SSwithin <- sum (sqrtDevW0$sqrtDevW0 + sqrtDevW1$sqrtDevW1  + sqrtDevW5$sqrtDevW5 +
                   sqrtDevW10$sqrtDevW10)

### Step 2: SS between (or conditions)
# the group sample size multiplied by the square different between the group's mean and grand mean
metadata
summary_table2 <-metadata %>%
  summarize(
    sqrtDev0 = (mean0 - GrandMean) ^ 2, ###SS of the conditions
    sqrtDev1 = (mean1 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev5 = (mean5 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev10 = (mean10 - GrandMean) ^ 2##SS of the conditions
  ) 
#write.csv(summary_table,"summary_table.csv")
##We get the SS of the conditions (time points), you multiple by (n) number of subjects under each (ith) condition
metadata
s1 <- 8 * summary_table2$sqrtDev0
s2 <- 8 * summary_table2$sqrtDev1
s3 <- 8 * summary_table2$sqrtDev5
s4 <- 9 * summary_table2$sqrtDev10
SSconditions <- sum(s1,s2,s3,s4)
SSconditions

## Now, let's calculate the SStotal= SSbetween + SS within
SStotal <- SSwithin + SSconditions

## Step 3-- SSsubjects
#we need to calculate the mean of each subject across the 3 timepoints and then you subtract from the grand mean and you multiply by the total number of samples
summary_table3 <-metadata %>%
  group_by(calf, observations) %>%
  summarize(
    sqrtDevS = (average - GrandMean) ^ 2, ###SS of the conditions
  )
#now we multiply that number by the total number of observations per animal
SSsubject <- summary_table3 %>%
  group_by(calf) %>%
  summarize(
    sqrtDevS2 = sqrtDevS * observations, ###SS of the conditions, the animals have different number of observations
  ) 
SSsubject <- sum(SSsubject$sqrtDevS2)

##Step 4- SSErroc
SSerror <- SSwithin - SSsubject
#now the double check
SSwithin
total <- SSsubject + SSerror

###now we move into calculating the degrees of freedom
#Now we specify the degrees of freedom of the model and for the error term
#dfmodel = k - 1, k number of variables
#dferror = (k-1) * (n-1), n is the number of subject
dfbetween = 4 - 1 #we have 3 variable or points
dfwithin = 33 - 4 #(total number of data points -number of groups)
dfsubject = 10 - 1 #(total number of subjects -1)
dferror = dfwithin - dfsubject
dftotal = 33 - 1 #total number of data points - 1

##now we calculate the Mean Square (SS for each df)
MSconditions = SSconditions / dfbetween
MSwithin = SSwithin / dfwithin
MSsubject = SSsubject / dfsubject
MSerror = SSerror / dferror

### now we calculate the F-ratio
Fratio = MSconditions / MSerror
# we then use the F-distribution with the degrees of freedom of the between and error
# now we calculate the P.value
p.value <- pf(Fratio, dfbetween, dferror, lower.tail = FALSE)

#creating the data frame
Response <- c("SSconditions", "SSwithin", "SSsubject", "SSerror", "SStotal")
SS <- c(SSconditions, SSwithin, SSsubject, SSerror, SStotal)
DF <- c(dfbetween, dfwithin, dfsubject, dferror, "NA")
MS <- c(MSconditions, MSwithin, MSsubject, MSerror, "NA")
F.value <- c(Fratio, "NA","NA","NA","NA")
P <- c(p.value, "NA","NA","NA","NA")

df3 <- data.frame(Response, SS, DF, MS, F.value, P)
print (df3) #bray curtis

#We need to import the distance matrices
DistWU = read.csv("dist.brayS.csv", header=TRUE)
row.names(DistWU) <- DistWU[,1]
DistWU = DistWU[,-c(1)]
distWU <- as.dist(DistWU)
#check homogeneity of variance
ant <-ant[rownames(ant) %in% rownames(DistWU),]#32

timeD1 <- betadisper(distWU, type = c("centroid"), group = ant$d)
timeD1
distances <- timeD1[["distances"]]
distances <- as.data.frame(distances)
distances <- merge(distances, ant, by.x = 0, by.y=0) # we are combining the data and removing unneccesary data
distances <- distances[,-c(1, 4,5, 7:9,11:13)]
#write.csv(distances, "distancesBCS.csv")

#upload the formatted table
metadata <- read.csv("betadisperCalfBCS.csv", header=TRUE)
str(metadata)
x <- metadata$s0
x <- na.exclude(x) 
mean0 <- mean(x)

y <- metadata$s1
y <- na.exclude(y) 
mean1 <- mean(y)

z <- metadata$s5
z <- na.exclude(z)
mean5 <- mean(z)

w <- metadata$s10
w <- na.exclude(w)
mean10 <- mean(w)

GrandMean <- (mean0 + mean1 + mean5 + mean10) / 4

summary_table <-metadata %>%
  summarize(
    sqrtDevW0 = (s0 - mean0) ^ 2, ###SS of the conditions
    sqrtDevW1 = (s1 - mean1) ^ 2, ###SS of the conditions
    sqrtDevW5 = (s5 - mean5) ^ 2,
    sqrtDevW10 = (s10 - mean10) ^ 2, ###SS of the conditions
  )
##removing NAs
sqrtDevW0 <- as.data.frame(summary_table$sqrtDevW0)
sqrtDevW0 <- as.data.frame(sqrtDevW0[-c(3,7),])
colnames(sqrtDevW0) <- c("sqrtDevW0")

sqrtDevW1 <- as.data.frame(summary_table$sqrtDevW1)
sqrtDevW1 <- as.data.frame(sqrtDevW1[-c(1,7),])
colnames(sqrtDevW1) <- c("sqrtDevW1")

sqrtDevW5 <- as.data.frame(summary_table$sqrtDevW5)
sqrtDevW5 <- as.data.frame(sqrtDevW5[-c(1,2),])
colnames(sqrtDevW5) <- c("sqrtDevW5")

sqrtDevW10 <- as.data.frame(summary_table$sqrtDevW10)
sqrtDevW10 <- as.data.frame(sqrtDevW10[-c(1),])
colnames(sqrtDevW10) <- c("sqrtDevW10")

SSwithin <- sum (sqrtDevW0$sqrtDevW0 + sqrtDevW1$sqrtDevW1  + sqrtDevW5$sqrtDevW5 +
                   sqrtDevW10$sqrtDevW10)

### Step 2: SS between (or conditions)
# the group sample size multiplied by the square different between the group's mean and grand mean
metadata
summary_table2 <-metadata %>%
  summarize(
    sqrtDev0 = (mean0 - GrandMean) ^ 2, ###SS of the conditions
    sqrtDev1 = (mean1 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev5 = (mean5 - GrandMean) ^ 2,##SS of the conditions
    sqrtDev10 = (mean10 - GrandMean) ^ 2##SS of the conditions
  ) 
#write.csv(summary_table,"summary_table.csv")
##We get the SS of the conditions (time points), you multiple by (n) number of subjects under each (ith) condition
metadata
s1 <- 8 * summary_table2$sqrtDev0
s2 <- 8 * summary_table2$sqrtDev1
s3 <- 8 * summary_table2$sqrtDev5
s4 <- 9 * summary_table2$sqrtDev10
SSconditions <- sum(s1,s2,s3,s4)
SSconditions

## Now, let's calculate the SStotal= SSbetween + SS within
SStotal <- SSwithin + SSconditions

## Step 3-- SSsubjects
#we need to calculate the mean of each subject across the 3 timepoints and then you subtract from the grand mean and you multiply by the total number of samples
summary_table3 <-metadata %>%
  group_by(calf, observations) %>%
  summarize(
    sqrtDevS = (average - GrandMean) ^ 2, ###SS of the conditions
  )
#now we multiply that number by the total number of observations per animal
SSsubject <- summary_table3 %>%
  group_by(calf) %>%
  summarize(
    sqrtDevS2 = sqrtDevS * observations, ###SS of the conditions, the animals have different number of observations
  ) 
SSsubject <- sum(SSsubject$sqrtDevS2)

##Step 4- SSErroc
SSerror <- SSwithin - SSsubject
#now the double check
SSwithin
total <- SSsubject + SSerror

###now we move into calculating the degrees of freedom
#Now we specify the degrees of freedom of the model and for the error term
#dfmodel = k - 1, k number of variables
#dferror = (k-1) * (n-1), n is the number of subject
dfbetween = 4 - 1 #we have 3 variable or points
dfwithin = 33 - 4 #(total number of data points -number of groups)
dfsubject = 10 - 1 #(total number of subjects -1)
dferror = dfwithin - dfsubject
dftotal = 33 - 1 #total number of data points - 1

##now we calculate the Mean Square (SS for each df)
MSconditions = SSconditions / dfbetween
MSwithin = SSwithin / dfwithin
MSsubject = SSsubject / dfsubject
MSerror = SSerror / dferror

### now we calculate the F-ratio
Fratio = MSconditions / MSerror
# we then use the F-distribution with the degrees of freedom of the between and error
# now we calculate the P.value
p.value <- pf(Fratio, dfbetween, dferror, lower.tail = FALSE)

#creating the data frame
Response <- c("SSconditions", "SSwithin", "SSsubject", "SSerror", "SStotal")
SS <- c(SSconditions, SSwithin, SSsubject, SSerror, SStotal)
DF <- c(dfbetween, dfwithin, dfsubject, dferror, "NA")
MS <- c(MSconditions, MSwithin, MSsubject, MSerror, "NA")
F.value <- c(Fratio, "NA","NA","NA","NA")
P <- c(p.value, "NA","NA","NA","NA")

df4 <- data.frame(Response, SS, DF, MS, F.value, P)
print (df4) #bray curtis

##respons
df1
df2
df3
df4
