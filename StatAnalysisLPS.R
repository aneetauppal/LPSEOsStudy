#Aneeta Uppal
#May 2018
#Code for statistical Analysis of RT-qPCR data

library(multcomp)

#input files contained fold changes from both donors. Each input file was seperated by
#inflammatory marker that was measured. Data frame was 10 columns x 10 rows. Each column was the treatment application
#rows 1:5 were measurments from donor 2 according to each application. Rows 6:10 were measurements from donor 1 per application

COX2_raw <- read.csv("/Users/aneetauppal/Graduate_PhD/DissertationResearch/Inflammation_Study/COX2/COX2FoldChangeD1D2Raw.csv")
CXCL10_raw <- read.csv("/Users/aneetauppal/Graduate_PhD/DissertationResearch/Inflammation_Study/CXCL10/CXCL10FoldChangeD2D1Raw.csv")
IL1B_raw<- read.csv("/Users/aneetauppal/Graduate_PhD/DissertationResearch/Inflammation_Study/IL1B/IL1BFoldChangeD2D1Raw.csv")
IL8_raw<-read.csv("/Users/aneetauppal/Graduate_PhD/DissertationResearch/Inflammation_Study/IL8/IL8FoldChangeD2D1Raw.csv")

#log transform data
COX2<-log10(COX2_raw)
CXCL10<-log10(CXCL10_raw)
IL1B <- log10(IL1B_raw)
IL8<-log10(IL8_raw)


#check for normality for donor 2
for (i in seq(1, 10, 1)){
  print(shapiro.test(COX2[1:5,i]))
  print(shapiro.test(CXCL10[1:5,i]))
  print(shapiro.test(IL1B[1:5,i]))
  print(shapiro.test(IL8[1:5,i]))
}

#check for normality for donor 1
for (i in seq(1, 10, 1)){
  print(shapiro.test(COX2[6:10,i]))
  print(shapiro.test(CXCL10[6:10,i]))
  print(shapiro.test(IL1B[6:10,i]))
  print(shapiro.test(IL8[6:10,i]))
}


#check for outliers using studentized residuals, remove any +3/-3
#change name labels according to what application is being looked at against the LPS controls
#change values according the application column that is being looked at
names = c(rep("LPS",5), rep("FCOsys",5))
value = c(CXCL10[1:5,2],CXCL10[1:5,3])
value = c(COX2[1:5,2],COX2[1:5,3])
value = c(IL8[1:5,2],IL8[1:5,3])
value = c(IL1B[1:5,2],IL1B[1:5,3])

#put data into linear model before look at rstudent values
datalm = data.frame(names,value)
summary(lm(value ~ names, data = datalm))
rstudent(lm(value ~ names, data = datalm))
datalm

#####
#Dunnetts test for COX2  D1
treatment = c(rep("LPS",5), rep("COS",5), rep("PS",5), rep("CS",5), rep("BS",5), rep("COT", 5), rep("PT",5), rep("CT",5), rep("BT", 5))
value = c(data[6:10,2],data[6:10,3], data[6:10,4], data[6:10,5], data[6:10,6], data[6:10,7], data[6:10,8], data[6:10,9], data[6:10,10])

datalm = data.frame(treatment,value)
levels(datalm$treatment)

datalm$treatment <- ordered(datalm$treatment,
                            levels = c( "LPS", "CS", "PS", "COS", "BS",  "CT",  "PT", "COT", "BT"))

datasumm <- aov(value ~ treatment, data = datalm)
summary(datasumm)
dunnett <- (glht(datasumm,  linfct=mcp(treatment="Dunnett")))
COXD1Results <- summary(dunnett)

models = lm(value ~ treatment, data = datalm)
plotted = resid(models)

plot(density(resid(models)))
qqnorm(resid(models))
qqline(resid(models))

#####
#Dunnetts test for COX2 D2
treatment = c(rep("LPS",5), rep("COS",5), rep("PS",5), rep("CS",5), rep("BS",5), rep("COT", 5), rep("PT",5), rep("CT",5), rep("BT", 5))
value = c(data[1:5,2],data[1:5,3], data[1:5,4], data[1:5,5], data[1:5,6], data[1:5,7], data[1:5,8], data[1:5,9], data[1:5,10])

datalm = data.frame(treatment,value)
levels(datalm$treatment)

datalm$treatment <- ordered(datalm$treatment,
                            levels = c( "LPS", "CS", "PS", "COS", "BS",  "CT",  "PT", "COT", "BT"))


datasumm <- aov(value ~ treatment, data = datalm)
summary(datasumm)
dunnett <- (glht(datasumm,  linfct=mcp(treatment="Dunnett")))
CoxD2Results<-summary(dunnett)

#####
#Dunnetts test for IL1B D2
treatment = c(rep("LPS",5), rep("COS",5), rep("PS",5), rep("CS",5), rep("BS",5), rep("COT", 5), rep("PT",5), rep("CT",5), rep("BT", 5))
value = c(IL1B[1:5,2],IL1B[1:5,3], IL1B[1:5,4],IL1B[1:5,5], IL1B[1:5,6], IL1B[1:5,7], IL1B[1:5,8], IL1B[1:5,9], IL1B[1:5,10])

datalm = data.frame(treatment,value)
levels(datalm$treatment)

datalm$treatment <- ordered(datalm$treatment,
                            levels = c( "LPS", "CS", "PS", "COS", "BS",  "CT",  "PT", "COT", "BT"))


datasumm <- aov(value ~ treatment, data = datalm)
summary(datasumm)
dunnett <- (glht(datasumm,  linfct=mcp(treatment="Dunnett")))
IL1BD2Results<-summary(dunnett)

#####
#Dunnetts test for IL1B D1
treatment = c(rep("LPS",5), rep("COS",5), rep("PS",5), rep("CS",5), rep("BS",5), rep("COT", 5), rep("PT",5), rep("CT",5), rep("BT", 5))
value = c(IL1B[6:10,2],IL1B[6:10,3], IL1B[6:10,4],IL1B[6:10,5], IL1B[6:10,6], IL1B[6:10,7], IL1B[6:10,8], IL1B[6:10,9], IL1B[6:10,10])

datalm = data.frame(treatment,value)
levels(datalm$treatment)

datalm$treatment <- ordered(datalm$treatment,
                            levels = c( "LPS", "CS", "PS", "COS", "BS",  "CT",  "PT", "COT", "BT"))


datasumm <- aov(value ~ treatment, data = datalm)
summary(datasumm)
dunnett <- (glht(datasumm,  linfct=mcp(treatment="Dunnett")))
IL1BD1Results<-summary(dunnett)

####
#Dunnetts test for CXCL10 D2
treatment = c(rep("LPS",5), rep("COS",5), rep("PS",5), rep("CS",5), rep("BS",5), rep("COT", 5), rep("PT",5), rep("CT",5), rep("BT", 5))
value = c(CXCL10[1:5,2],CXCL10[1:5,3], CXCL10[1:5,4],CXCL10[1:5,5],CXCL10[1:5,6], CXCL10[1:5,7], CXCL10[1:5,8], CXCL10[1:5,9], CXCL10[1:5,10])

datalm = data.frame(treatment,value)
levels(datalm$treatment)

datalm$treatment <- ordered(datalm$treatment,
                            levels = c( "LPS", "CS", "PS", "COS", "BS",  "CT",  "PT", "COT", "BT"))


datasumm <- aov(value ~ treatment, data = datalm)
summary(datasumm)
dunnett <- (glht(datasumm,  linfct=mcp(treatment="Dunnett")))
CXCL10D2Results<-summary(dunnett)

####
#Dunnetts test for CXCL10 D1
treatment = c(rep("LPS",5), rep("COS",5), rep("PS",5), rep("CS",5), rep("BS",5), rep("COT", 5), rep("PT",5), rep("CT",5), rep("BT", 5))
value = c(CXCL10[6:10,2],CXCL10[6:10,3], CXCL10[6:10,4],CXCL10[6:10,5],CXCL10[6:10,6], CXCL10[6:10,7], CXCL10[6:10,8], CXCL10[6:10,9], CXCL10[6:10,10])

datalm = data.frame(treatment,value)
levels(datalm$treatment)

datalm$treatment <- ordered(datalm$treatment,
                            levels = c( "LPS", "CS", "PS", "COS", "BS",  "CT",  "PT", "COT", "BT"))


datasumm <- aov(value ~ treatment, data = datalm)
summary(datasumm)
dunnett <- (glht(datasumm,  linfct=mcp(treatment="Dunnett")))
CXCL10D1Results <- summary(dunnett)

####
#Dunnetts test for IL8 D1
treatment = c(rep("LPS",5), rep("COS",5), rep("PS",5), rep("CS",5), rep("BS",5), rep("COT", 5), rep("PT",5), rep("CT",5), rep("BT", 5))
value = c(IL8[6:10,2],IL8[6:10,3], IL8[6:10,4],IL8[6:10,5],IL8[6:10,6], IL8[6:10,7], IL8[6:10,8], IL8[6:10,9], IL8[6:10,10])

datalm = data.frame(treatment,value)
levels(datalm$treatment)

datalm$treatment <- ordered(datalm$treatment,
                            levels = c( "LPS", "CS", "PS", "COS", "BS",  "CT",  "PT", "COT", "BT"))


datasumm <- aov(value ~ treatment, data = datalm)
summary(datasumm)
dunnett <- (glht(datasumm,  linfct=mcp(treatment="Dunnett")))
IL8D1Results<-summary(dunnett)

####
#Dunnetts test for IL8 D2
treatment = c(rep("LPS",5), rep("COS",5), rep("PS",5), rep("CS",5), rep("BS",5), rep("COT", 5), rep("PT",5), rep("CT",5), rep("BT", 5))
value = c(IL8[1:5,2],IL8[1:5,3], IL8[1:5,4],IL8[1:5,5],IL8[1:5,6], IL8[1:5,7], IL8[1:5,8], IL8[1:5,9], IL8[1:5,10])

datalm = data.frame(treatment,value)
levels(datalm$treatment)

datalm$treatment <- ordered(datalm$treatment,
                            levels = c( "LPS", "CS", "PS", "COS", "BS",  "CT",  "PT", "COT", "BT"))


datasumm <- aov(value ~ treatment, data = datalm)
summary(datasumm)
dunnett <- (glht(datasumm,  linfct=mcp(treatment="Dunnett")))
IL8D2Results<-summary(dunnett)
