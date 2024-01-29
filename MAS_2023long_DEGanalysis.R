#load libraries
library(ggplot2)
library(ggvis)
library(plyr)
library(dplyr)
library(stringr)
library(tibble)
library(locfit)
library(scales)
library(edgeR)
library(DESeq2)
library(viridis)
library(sva)
library(glmmSeq)
library(openxlsx)
library(PCAtools)
library(Biobase)
library(ggpubr)

#working directory
setwd("C:/Users/mscot/Dropbox/Research/2022 Longitudinal")

#load in raw count data
##had to add data.matrix with data.w.rownames; kept running into data matrix error...
raw_data = read.csv("2023_long_BRD.csv",header=TRUE, check.names = TRUE)
data.w.rownames <- data.matrix(data.frame(raw_data[,-1], row.names=raw_data[,1]))

#load in metadata
traitdata = read.csv("metadata_2023long.csv")

#just to make life easier later on...
day<-factor(rep(c(traitdata$Day)))
severity<-factor(rep(c(traitdata$Severity)))
platform_all<-factor(rep(c(traitdata$Platform)))
year<-factor(rep(c(traitdata$Year)))
animalID=traitdata$A_ID
Timing<-factor(rep(c(traitdata$BRD_Risk)))

###use of ComBat-Seq to account for sequencing batch effect
##ComBat-seq takes untransformed, raw count matrix as input (RTM :-) ); will pre-process afterwards
adjusted<-ComBat_seq(data.w.rownames, batch = traitdata$Platform, group = traitdata$Severity)
head(data.w.rownames)
head(adjusted)

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

#create read matrix
dgedata <- DGEList(counts=y, group=traitdata$Severity)
dgedata$samples

#evaluating library sizes by sample_ID
barplot(dgedata$samples$lib.size, names=colnames(dgedata),las=2, cex.names=.55, col=plasma(3)[severity])
title("Total Library Count per Sample")
abline(h=median(dgedata$samples$lib.size),col="black")

###evaluate for DEGs via glmmSeq###
dgedata_glmm <- DGEList(counts=y)

# glmmSeq -----------------------------------------------------------------
#calculate normalization factors via edgeR TMM method
glmm_normal <- calcNormFactors(y, method = c("TMM"))
glmm_normal

#estimate tagwise dispersion (edgeR method)
edisp2<-setNames(edgeR::estimateDisp(dgedata_glmm)$tagwise.dispersion, rownames(dgedata_glmm))
edisp2

#running glmm model - adjust cores as necessary
##primary interest is in Day:Severity
results <- glmmSeq(~ Day * Severity + (1 | actual_ID) + (1 | year) + (1 | Vax),
                   id = "actual_ID",
                   countdata = dgedata$counts,
                   metadata = traitdata,
                   dispersion = edisp2,
                   sizeFactors = glmm_normal,
                   removeDuplicatedMeasures = FALSE,
                   removeSingles=FALSE,
                   progress=TRUE,
                   cores = 8)

###runtime was ~15 minutes
names(attributes(results))
results@stats

#calculate Qvalues; pi0 adjusted to 1 to enact BH procedure
qvals<-glmmQvals(results, pi0 = 1)
qvals@stats

predict = data.frame(results@predict)
stats = data.frame(results@stats)
glmdegs = data.frame(qvals@stats)

#write analysis table to csv file
write.csv(stats, file = "glm_stats.csv")
write.csv(glmdegs, file = "glm_degs.csv")

#starting to plot genes of interest
pairedPlot(glmmResult=qvals,
           geneName = "CATHL2",
           x1Label = "Day",
           x2Label="Severity",
           xTitle="Day",
           yTitle = "log10 Gene Expression",
           IDColumn = "actual_ID",
           graphics = "ggplot",
           colours = magma(3, begin = 0, end = 0.75),
           modelColour = mako(9, direction = -1),
           modelLineColour = mako(9, begin = 0, end = 0.75, direction = -1),
           fontSize=10,
           x2Offset = 8,
           logTransform=TRUE,
           addViolin = TRUE,
           pairedOnly = FALSE)


modelPlot(glmmResult=qvals,
          "CATHL1",
          x1Label="Day",
          x2Label="Severity",
          xTitle="Day",
          yTitle = "log10 Gene Expression",
          fontSize=8,
          x2Offset=1,
          addErrorbars = FALSE,
          overlap=TRUE,
          graphics="ggplot",
          logTransform = TRUE,
          colours = magma(3, begin = 0, end = 0.75))

modelPlot(glmmResult=qvals,
          "CATHL1",
          x1Label="Day",
          x2Label="Severity",
          xTitle="Day",
          yTitle = "log10 Gene Expression",
          fontSize=8,
          x2Offset=1,
          addErrorbars = TRUE,
          overlap=TRUE,
          graphics="ggplot",
          logTransform = TRUE,
          colours = magma(3, begin = 0, end = 0.75))



p1 = pairedPlot(glmmResult=qvals,
                geneName = "CATHL1",
                x1Label = "Day",
                x2Label="Severity",
                xTitle="Day",
                yTitle = "log10 Gene Expression",
                IDColumn = "actual_ID",
                graphics = "ggplot",
                colours = magma(3, begin = 0, end = 0.75),
                modelColour = mako(9, direction = -1),
                modelLineColour = mako(9, begin = 0, end = 0.75, direction = -1),
                fontSize=10,
                x2Offset = 8,
                logTransform=TRUE,
                addViolin = TRUE,
                pairedOnly = FALSE)

p2 = modelPlot(glmmResult=qvals,
               "CATHL2",
               x1Label="Day",
               x2Label="Severity",
               xTitle="Day",
               yTitle = "log10 Gene Expression",
               fontSize=8,
               x2Offset=1,
               addErrorbars = FALSE,
               overlap=TRUE,
               graphics="ggplot",
               logTransform = TRUE,
               colours = magma(3, begin = 0, end = 0.75))

ggarrange(p1, p2, ncol=2, common.legend = T, legend="bottom")

labels = c('MS4A2', 'ALOX15', 'IL16', 'CFB', 'GZMB',
           'CATHL2', 'CATHL3', 'CXCL10')
maPlots <- maPlot(results,
                  x1Label="Day",
                  x2Label="Severity",
                  colours=c('grey', 'midnightblue',
                            'mediumseagreen', 'goldenrod'),
                  labels = labels,
                  graphics="ggplot")

maPlots$combined



# edgeR -------------------------------------------------------------------
###utilize quasi-likelihood F-tests as a means of post-hoc comparative analysis
##defining experimental factors within reduced model
#Setting up the matrix to block for paired samples
dgedatanormal <- calcNormFactors(dgedata, method = c("TMM"))
platform_all<-traitdata$Platform
design<-model.matrix(~0+animalID+year)

#not elegant, but effective; setting up for analysis between cohorts
healthy_d0 = severity=='Healthy' & day=='0'
healthy_d28 <- severity=='Healthy' & day=='28'
healthy_d63 <- severity=='Healthy' & day=='63'
T1_d0 <- severity=='T1' & day=='0'
T1_d28 <- severity=='T1' & day=='28'
T1_d63 <- severity=='T1' & day=='63'
T2_d0 <- severity=='T2' & day=='0'
T2_d28 <- severity=='T2' & day=='28'
T2_d63 <- severity=='T2' & day=='63'

design<-cbind(design, healthy_d0, healthy_d28, healthy_d63, T1_d0, T1_d28,
              T1_d63, T2_d0, T2_d28, T2_d63)
design

rownames(design) <- colnames(dgedatanormal)
design
colnames(design)

#estimate genewise dispersion for abundance trends; over 6 replicates per cohort 
#so should not need robust function
#"Since the NB dispersion under the QL framework reflects the overall biological
#variability, it does not make sense to use the tagwise dispersions." (2.9.4)
total <- estimateDisp(dgedatanormal, design)
plotBCV(total)
total$common.dispersion

#QL dispersions are estimated and visualized
fit<-glmQLFit(total, design)
plotQLDisp(fit)
fit$offset
fit$fitted.values

#makeContrasts function for pairwise evaluation
BRD_contrasts <- makeContrasts(
  hd0vT1d0 = healthy_d0-T1_d0,
  hd28vT1d28 = healthy_d28-T1_d28,
  hd63vT1d63 = healthy_d63-T1_d63,
  hd0vT2d0 = healthy_d0-T2_d0,
  hd28vT2d28 = healthy_d28-T2_d28,
  hd63vT2d63 = healthy_d63-T2_d63,
  T1d0vT2d0 = T1_d0-T2_d0,
  T1d28vT2d28 = T1_d28-T2_d28,
  T1d63vT2d63 = T1_d63-T2_d63,
  levels = design)

###QLF differential expression analysis; FDR a priori at 0.10

#Healthy vs T1
qlfhd0vT1d0<-glmQLFTest(fit, contrast = BRD_contrasts[,"hd0vT1d0"])
topTags(qlfhd0vT1d0)
test_hd0vT1d0<-topTags(qlfhd0vT1d0, n=55000)
summary(decideTests.DGELRT(qlfhd0vT1d0,p.value = 0.1))
test_hd0vT1d0

qlfhd28vT1d28<-glmQLFTest(fit, contrast = BRD_contrasts[,"hd28vT1d28"])
topTags(qlfhd28vT1d28)
test_hd28vT1d28<-topTags(qlfhd28vT1d28, n=55000)
summary(decideTests.DGELRT(qlfhd28vT1d28,p.value = 0.1))
test_hd28vT1d28

qlfhd63vT1d63<-glmQLFTest(fit, contrast = BRD_contrasts[,"hd63vT1d63"])
topTags(qlfhd63vT1d63)
test_hd63vT1d63<-topTags(qlfhd63vT1d63, n=55000)
summary(decideTests.DGELRT(qlfhd63vT1d63,p.value = 0.1))
test_hd63vT1d63

wb_HvT1<- createWorkbook("HvT1")
addWorksheet(wb_HvT1, "hd0vT1d0")
addWorksheet(wb_HvT1, "hd28vT1d28")
addWorksheet(wb_HvT1, "hd63vT1d63")

writeData(wb_HvT1, sheet = "hd0vT1d0", test_hd0vT1d0$table, rowNames = TRUE)
writeData(wb_HvT1, sheet = "hd28vT1d28", test_hd28vT1d28$table, rowNames = TRUE)
writeData(wb_HvT1, sheet = "hd63vT1d63", test_hd63vT1d63$table, rowNames = TRUE)
saveWorkbook(wb_HvT1, "HvT1_QLF.xlsx")

#Healthy vs T2
qlfhd0vT2d0<-glmQLFTest(fit, contrast = BRD_contrasts[,"hd0vT2d0"])
topTags(qlfhd0vT2d0)
test_hd0vT2d0<-topTags(qlfhd0vT2d0, n=55000)
summary(decideTests.DGELRT(qlfhd0vT2d0,p.value = 0.1))
test_hd0vT2d0

qlfhd28vT2d28<-glmQLFTest(fit, contrast = BRD_contrasts[,"hd28vT2d28"])
topTags(qlfhd28vT2d28)
test_hd28vT2d28<-topTags(qlfhd28vT2d28, n=55000)
summary(decideTests.DGELRT(qlfhd28vT2d28,p.value = 0.1))
test_hd28vT2d28

qlfhd63vT2d63<-glmQLFTest(fit, contrast = BRD_contrasts[,"hd63vT2d63"])
topTags(qlfhd63vT2d63)
test_hd63vT2d63<-topTags(qlfhd63vT2d63, n=55000)
summary(decideTests.DGELRT(qlfhd63vT2d63,p.value = 0.1))
test_hd63vT2d63

wb_HvT2<- createWorkbook("HvT2")
addWorksheet(wb_HvT2, "hd0vT2d0")
addWorksheet(wb_HvT2, "hd28vT2d28")
addWorksheet(wb_HvT2, "hd63vT2d63")

writeData(wb_HvT2, sheet = "hd0vT2d0", test_hd0vT2d0$table, rowNames = TRUE)
writeData(wb_HvT2, sheet = "hd28vT2d28", test_hd28vT2d28$table, rowNames = TRUE)
writeData(wb_HvT2, sheet = "hd63vT2d63", test_hd63vT2d63$table, rowNames = TRUE)
saveWorkbook(wb_HvT2, "HvT2_QLF.xlsx")

#T1 vs T2
qlfT1d0vT2d0<-glmQLFTest(fit, contrast = BRD_contrasts[,"T1d0vT2d0"])
topTags(qlfT1d0vT2d0)
test_T1d0vT2d0<-topTags(qlfT1d0vT2d0, n=55000)
summary(decideTests.DGELRT(qlfT1d0vT2d0,p.value = 0.1))
test_T1d0vT2d0

qlfT1d28vT2d28<-glmQLFTest(fit, contrast = BRD_contrasts[,"T1d28vT2d28"])
topTags(qlfT1d28vT2d28)
test_T1d28vT2d28<-topTags(qlfT1d28vT2d28, n=55000)
summary(decideTests.DGELRT(qlfT1d28vT2d28,p.value = 0.1))
test_T1d28vT2d28

qlfT1d63vT2d63<-glmQLFTest(fit, contrast = BRD_contrasts[,"T1d63vT2d63"])
topTags(qlfT1d63vT2d63)
test_T1d63vT2d63<-topTags(qlfT1d63vT2d63, n=55000)
summary(decideTests.DGELRT(qlfT1d63vT2d63,p.value = 0.1))
test_T1d63vT2d63

wb_T1VT2<- createWorkbook("T1vT2")
addWorksheet(wb_T1VT2, "T1d0vT2d0")
addWorksheet(wb_T1VT2, "T1d28vT2d28")
addWorksheet(wb_T1VT2, "T1d63vT2d63")

writeData(wb_T1VT2, sheet = "T1d0vT2d0", test_T1d0vT2d0$table, rowNames = TRUE)
writeData(wb_T1VT2, sheet = "T1d28vT2d28", test_T1d28vT2d28$table, rowNames = TRUE)
writeData(wb_T1VT2, sheet = "T1d63vT2d63", test_T1d63vT2d63$table, rowNames = TRUE)
saveWorkbook(wb_T1VT2, "T1vT2_QLF.xlsx")

