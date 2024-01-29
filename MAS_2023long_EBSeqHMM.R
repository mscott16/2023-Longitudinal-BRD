#load libraries
library(ggplot2)
library(ggvis)
library(dplyr)
library(stringr)
library(tibble)
library(locfit)
library(scales)
library(edgeR)
library(plyr)
library(viridis)
library(sva)
library(glmmSeq)
library(openxlsx)
library(PCAtools)
library(Biobase)
library(EBSeq)  
library(EBSeqHMM)

#working directory
setwd("C:/Users/mscot/Dropbox/Research/2022 Longitudinal")

#load in raw count data
raw_data = read.csv("long22_notag_no14_v2.csv",header=TRUE, check.names = TRUE)
data = ddply(raw_data,"gene_id", numcolwise(sum))
data.w.rownames <- data.matrix(data.frame(data[,-1], row.names=data[,1]))

#load in metadata
traitdata = read.csv("Metadata_no14.csv")
traitdata$Platform
traitdata$GR

#just to make life easier later on...
day<-factor(rep(c(traitdata$Day)))
severity<-factor(rep(c(traitdata$Severity)))
platform<-factor(rep(c(traitdata$Platform)))
year<-factor(rep(c(traitdata$Year)))
animalID=traitdata$A_ID
Timing<-factor(rep(c(traitdata$BRD_Risk)))

###use of ComBat-Seq to account for sequencing batch effect
##ComBat-seq takes untransformed, raw count matrix as input (RTM :-) ); will pre-process afterwards
adjusted<-ComBat_seq(data.w.rownames, batch = traitdata$Platform, group = traitdata$Severity)
adjusted
#pre-processing low read counts
keepexper<- filterByExpr(adjusted, min.count = 0, min.total.count = 100, group = traitdata$Severity)
x<-adjusted[keepexper,]
dim(x)

sum(data.w.rownames)
sum(adjusted)
sum(x)

keep<-rowSums(cpm(x)>0.5) >= 12
y<-x[keep,]
dim(y)

levels(day)


##break down data set by severity; will need to pull quantile values

normed<-QuantileNorm(y,.75)
#Upper-Quantile Normalization in Bullard et al. (2010))
normed
write.csv(normed, "normed.csv")
write.csv(y, "y.csv")

#From here, the normalized count table is split by severity group (H, T1, T2)
#Not elegant, but it works...
H_days<-factor(rep(c("0", "28", "28", "28", "63", "63", "63", "28", "28", "0", 
                     "28", "28", "28", "28", "28", "28", "28", "28", "28", "63", 
                     "63", "63", "63", "63", "63", "63", "63", "63", "63", "63", 
                     "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0")))
T1_days<-factor(rep(c("63", "28", "28", "0", "28", "28", "28", "28", "28", 
                      "28", "28", "28", "28", "63", "63", "63", "63", "63", 
                      "63", "63", "63", "63", "63", "0", "0", "0", "0", "0",
                      "0", "0", "0", "0", "0")))
T2_days<-factor(rep(c("28", "28", "0", "63", "63", "63", "63", "0", "0", "0", 
                      "0", "28", "28", "28", "28", "28", "0", "28", "63", "63",
                      "0", "0", "0", "0", "0")))

H_counts<-read.csv("H_adjust.csv",header=TRUE, check.names = TRUE)
T1_counts<-read.csv("T1_adjust.csv",header=TRUE, check.names = TRUE)
T2_counts<-read.csv("T2_adjust.csv",header=TRUE, check.names = TRUE)

H_rownames <- data.matrix(data.frame(H_counts[,-1], row.names=H_counts[,1]))
T1_rownames <- data.matrix(data.frame(T1_counts[,-1], row.names=T1_counts[,1]))
T2_rownames <- data.matrix(data.frame(T2_counts[,-1], row.names=T2_counts[,1]))


H_normed <-read.csv("H_normed.csv",header=TRUE, check.names = TRUE)
T1_normed <-read.csv("T1_normed.csv",header=TRUE, check.names = TRUE)
T2_normed <-read.csv("T2_normed.csv",header=TRUE, check.names = TRUE)


HH<-as.numeric(H_normed$x)
T1T1<-as.numeric(T1_normed$x)
T2T2<-as.numeric(T2_normed$x)

GeneNormData <- GetNormalizedMat(T2_rownames, T2T2)
PlotExp(GeneNormData, T2_days, Name="ALOX15")

H_med<-MedianNorm(H_rownames, alternative = FALSE)
T1_med<-MedianNorm(T1_rownames, alternative = FALSE)
T2_med <- MedianNorm(T2_rownames, alternative = FALSE)


#T2 group

T2s_out<-EBSeqHMMTest(Data = T2_rownames, sizeFactors = T2_med, Conditions = T2_days, UpdateRd = 100,
                      UpdatePI = FALSE, FCV = 2.4)

GeneDECalls <- GetDECalls(T2s_out, FDR=.01)
head(GeneDECalls)
GeneDECalls
GeneDECalls["TARP",]

T2s_out$LLSum

#https://groups.google.com/g/ebseq-users/c/2mgW4zTuLEY/m/f65XW1mFBgAJ
#"You may compare the log likelihood between different runs with
#different numbers of iterations. e.g., look at EBSeqHMMGeneOut$LLSum
#If the LL doesn't improve much by increasing the number of iterations,
#then that means the model fitting based on the smaller number of
#iterations is 'good enough'. On the other side, if it improves a lot,
#you may consider further increasing the number of iterations."

QQP(T2s_out,GeneLevel=TRUE)
DenNHist(T2s_out, GeneLevel=TRUE)

PlotExp(GeneNormData, T2_days, Name="TARP")


GeneConfCalls <- GetConfidentCalls(T2s_out, FDR=.01,cutoff=.5, OnlyDynamic=TRUE)
head(GeneConfCalls$Overall)
str(GeneConfCalls$EachPathNames)

GeneConfCalls$EachPathNames$`Up-Down`
write.csv(GeneConfCalls$EachPathNames$`Up-Down`, "T2_Up-Down.csv")

GeneConfCalls$EachPathNames$`Up-Up`
write.csv(GeneConfCalls$EachPathNames$`Up-Up`, "T2_Up-Up.csv")

GeneConfCalls$EachPathNames$`Down-Down`
write.csv(GeneConfCalls$EachPathNames$`Down-Down`, "T2_Down-Down.csv")

GeneConfCalls$EachPathNames$
write.csv(GeneConfCalls$EachPathNames$`Down-Up`, "T2_Down-Up.csv")



#T1 group

T1s_out<-EBSeqHMMTest(Data = T1_rownames, sizeFactors = T1_med, Conditions = T1_days, UpdateRd = 100,
                      UpdatePI = FALSE, FCV = 2.4)

QQP(T1s_out,GeneLevel=TRUE)
DenNHist(T1s_out, GeneLevel=TRUE)

GeneDECalls_T1 <- GetDECalls(T1s_out, FDR=.01)
head(GeneDECalls_T1)
GeneDECalls_T1
GeneDECalls_T1["CFP",]

T1s_out$LLSum

GeneConfCalls_T1 <- GetConfidentCalls(T1s_out, FDR=.01,cutoff=.5, OnlyDynamic=TRUE)
head(GeneConfCalls_T1$Overall)
str(GeneConfCalls_T1$EachPathNames)

GeneConfCalls_T1$EachPathNames$`Up-Down`
write.csv(GeneConfCalls_T1$EachPathNames$`Up-Down`, "T1_Up-Down.csv")

GeneConfCalls_T1$EachPathNames$`Up-Up`
write.csv(GeneConfCalls_T1$EachPathNames$`Up-Up`, "T1_Up-Up.csv")

GeneConfCalls_T1$EachPathNames$`Down-Down`
write.csv(GeneConfCalls_T1$EachPathNames$`Down-Down`, "T1_Down-Down.csv")

GeneConfCalls_T1$EachPathNames$`Down-Up`
write.csv(GeneConfCalls_T1$EachPathNames$`Down-Up`, "T1_Down-Up.csv")


#Healthy group

H_out<-EBSeqHMMTest(Data = H_rownames, sizeFactors = H_med, Conditions = H_days, UpdateRd = 100,
                    UpdatePI = FALSE, FCV = 2.4)

GeneDECalls_H <- GetDECalls(H_out, FDR=.01)
head(GeneDECalls_H)
GeneDECalls_H
GeneDECalls_H["HPGD",]
H_out$LLSum

QQP(H_out,GeneLevel=TRUE)
DenNHist(H_out, GeneLevel=TRUE)

GeneConfCalls_H <- GetConfidentCalls(H_out, FDR=.01,cutoff=.5, OnlyDynamic=TRUE)
head(GeneConfCalls_H$NumEach)
str(GeneConfCalls_H$EachPathNames)

GeneConfCalls_H$EachPathNames$`Up-Down`
write.csv(GeneConfCalls_H$EachPathNames$`Up-Down`, "H_Up-Down.csv")

GeneConfCalls_H$EachPathNames$`Up-Up`
write.csv(GeneConfCalls_H$EachPathNames$`Up-Up`, "H_Up-Up.csv")

GeneConfCalls_H$EachPathNames$`Down-Down`
write.csv(GeneConfCalls_H$EachPathNames$`Down-Down`, "H_Down-Down.csv")

GeneConfCalls_H$EachPathNames$`Down-Up`
write.csv(GeneConfCalls_H$EachPathNames$`Down-Up`, "H_Down-Up.csv")

write.csv(GeneConfCalls_H$Overall, "H_all.csv")
write.csv(GeneConfCalls_T1$Overall, "T1_all.csv")
write.csv(GeneConfCalls$Overall, "T2_all.csv")

