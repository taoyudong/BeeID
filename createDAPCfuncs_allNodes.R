
setwd("C:/Users/donth/OneDrive/Kiran_projects/bee_SNP_analysis/genotyping/fluidigm/UIUC/All_plates/KnowBase_Agglom_clus")
getwd()


library(adegenet)
options(stringsAsFactors = F, scipen = 999)

trainPerc <- 0.9 # percentage of individuals for training per group
metaCol <- 8 # Number of columns with meta data

######## Root as base
root_data <- read.csv("Plates1to9Test_impRound_NoFrance1panama_w_metadata_root.csv")
root_gtypes <- root_data[, -c(1:metaCol)]
rownames(root_gtypes) <- root_data$SampleID
root_Xdapc<-xvalDapc(root_gtypes, root_data$AggloClusLab, center=TRUE, training.set = trainPerc, 
                     result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                    scale=FALSE)


 # To get contingency table
table(root_Xdapc$DAPC$assign,root_Xdapc$DAPC$grp)

# Identify samples that got misclassified by xvalDAPC and exclude them from running DAPC
# To get misclassified samples from xvalDapc
missCl <- root_Xdapc$DAPC$assign != root_Xdapc$DAPC$grp 
correctCl <- root_Xdapc$DAPC$assign == root_Xdapc$DAPC$grp
root_gtypesFil <- root_gtypes[correctCl,] # this now retains only samples that got classified as expected
root_dataFil <- root_data[correctCl,]
root_gtypesFil <- na.omit(root_gtypesFil)
root_mcSam <- root_data$SampleID[missCl]

optPcaNum <- as.numeric(root_Xdapc$`Number of PCs Achieving Lowest MSE`)
root_dapc <- dapc(root_gtypesFil, root_dataFil$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum,
                  scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                  pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)
table(root_dapc$assign,root_dapc$grp)
# Get misclassified samples from the DAPC run
missCl2 <- root_dapc$assign != root_dapc$grp 
correctCl2 <- root_dapc$assign == root_dapc$grp
root_gtypesFil2 <- root_gtypesFil[correctCl2,]
root_gtypesFil2 <- na.omit(root_gtypesFil2)
root_dataFil2 <- root_dataFil[correctCl2,]
root_mcSam <- c(root_mcSam, root_dataFil$SampleID[missCl2])


# re run xvaldapc and dapc after excluding misclassified samples
root_Xdapc2<-xvalDapc(root_gtypesFil2, root_dataFil2$AggloClusLab, center=TRUE, training.set = trainPerc, 
                     result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                     scale=FALSE)
table(root_Xdapc2$DAPC$assign,root_Xdapc2$DAPC$grp)
root_dapc2 <- dapc(root_gtypesFil2, root_dataFil2$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum,
                  scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                  pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)
table(root_dapc2$assign,root_dapc2$grp)



######## A_1 as base
library(dplyr)
A_1_data <- read.csv("Plates1to9Test_impRound_NoFrance1panama_w_metadata_A_1.csv")

# Exclude root misclassified samples from the downstream steps
A_1_data %>% filter(!SampleID %in% root_mcSam) -> A_1_dataRem

A_1_gtypes <- A_1_dataRem[, -c(1:metaCol)]

A_1_Xdapc<-xvalDapc(A_1_gtypes, A_1_dataRem$AggloClusLab, center=TRUE, training.set = trainPerc, 
                     result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                     scale=FALSE)
# To get contingency table
table(A_1_Xdapc$DAPC$assign,A_1_Xdapc$DAPC$grp)

## No samples got misclassified for A_1 xvaldapc

optPcaNum <- as.numeric(A_1_Xdapc$`Number of PCs Achieving Lowest MSE`)
A_1_dapc <- dapc(A_1_gtypes, A_1_dataRem$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum,
                  scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                  pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)
table(A_1_dapc$assign,A_1_dapc$grp)

## No samples got misclassified for A_1 dapc


######## A2 as base
A_2_data <- read.csv("Plates1to9Test_impRound_NoFrance1panama_w_metadata_A_2.csv")

# Exclude root misclassified samples from the downstream steps
A_2_data %>% filter(!SampleID %in% root_mcSam) -> A_2_dataRem

A_2_gtypes <- A_2_dataRem[, -c(1:metaCol)]

A_2_Xdapc<-xvalDapc(A_2_gtypes, A_2_dataRem$AggloClusLab, center=TRUE, training.set = trainPerc, 
                     result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                     scale=FALSE)


# To get contingency table
table(A_2_Xdapc$DAPC$assign,A_2_Xdapc$DAPC$grp)

# Identify samples that got misclassified by xvalDAPC and exclude them from running DAPC
# To get misclassified samples from xvalDapc
missCl <- A_2_Xdapc$DAPC$assign != A_2_Xdapc$DAPC$grp 
correctCl <- A_2_Xdapc$DAPC$assign == A_2_Xdapc$DAPC$grp
A_2_gtypesFil <- A_2_gtypes[correctCl,] # this now retains only samples that got classified as expected
A_2_dataFil <- A_2_dataRem[correctCl,]
A_2_mcSam <- A_2_dataRem$SampleID[missCl]

optPcaNum <- as.numeric(A_2_Xdapc$`Number of PCs Achieving Lowest MSE`)
A_2_dapc <- dapc(A_2_gtypesFil, A_2_dataFil$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum,
                  scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                  pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)
table(A_2_dapc$assign,A_2_dapc$grp)
# No misclassified samples from the DAPC run

# re run xvaldapc after excluding misclassified samples
A_2_Xdapc2<-xvalDapc(A_2_gtypesFil, A_2_dataFil$AggloClusLab, center=TRUE, training.set = trainPerc, 
                      result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                      scale=FALSE)
table(A_2_Xdapc2$DAPC$assign,A_2_Xdapc2$DAPC$grp)

# add up all misclassified samples so far
total_mcSam <- c(root_mcSam,A_2_mcSam)


####### B_1_2 as base
B_1_2_data <- read.csv("Plates1to9Test_impRound_NoFrance1panama_w_metadata_B_1_2.csv")

# Exclude root misclassified samples from the B_1_2_data
B_1_2_data %>% filter(!SampleID %in% root_mcSam) -> B_1_2_dataRem

B_1_2_gtypes <- B_1_2_dataRem[, -c(1:metaCol)]

B_1_2_Xdapc<-xvalDapc(B_1_2_gtypes, B_1_2_dataRem$AggloClusLab, center=TRUE, training.set = trainPerc, 
                     result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                     scale=FALSE)


# To get contingency table
table(B_1_2_Xdapc$DAPC$assign,B_1_2_Xdapc$DAPC$grp)

# No samples got misclassified by xvalDAPC

optPcaNum <- as.numeric(B_1_2_Xdapc$`Number of PCs Achieving Lowest MSE`)
B_1_2_dapc <- dapc(B_1_2_gtypes, B_1_2_dataRem$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum,
                  scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                  pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)
table(B_1_2_dapc$assign,B_1_2_dapc$grp)
# No misclassified samples from the DAPC run


######## B_2_1 as base
B_2_1_data <- read.csv("Plates1to9Test_impRound_NoFrance1panama_w_metadata_B_2_1.csv")

# Exclude misclassified samples from the B_2_1_data
B_2_1_data %>% filter(!SampleID %in% total_mcSam) -> B_2_1_dataRem

B_2_1_gtypes <- B_2_1_dataRem[, -c(1:metaCol)]

B_2_1_Xdapc<-xvalDapc(B_2_1_gtypes, B_2_1_dataRem$AggloClusLab, center=TRUE, training.set = trainPerc, 
                     result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                     scale=FALSE)

# To get contingency table
table(B_2_1_Xdapc$DAPC$assign,B_2_1_Xdapc$DAPC$grp)

# Identify samples that got misclassified by xvalDAPC and exclude them from running DAPC
# To get misclassified samples from xvalDapc
missCl <- B_2_1_Xdapc$DAPC$assign != B_2_1_Xdapc$DAPC$grp 
correctCl <- B_2_1_Xdapc$DAPC$assign == B_2_1_Xdapc$DAPC$grp
B_2_1_gtypesFil <- B_2_1_gtypes[correctCl,] # this now retains only samples that got classified as expected
B_2_1_dataFil <- B_2_1_dataRem[correctCl,]
B_2_1_mcSam <- B_2_1_dataRem$SampleID[missCl]

optPcaNum <- as.numeric(B_2_1_Xdapc$`Number of PCs Achieving Lowest MSE`)
B_2_1_dapc <- dapc(B_2_1_gtypesFil, B_2_1_dataFil$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum,
                  scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                  pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)
table(B_2_1_dapc$assign,B_2_1_dapc$grp)
# No misclassified samples from the DAPC run


# re run xvaldapc after excluding misclassified samples
B_2_1_Xdapc2<-xvalDapc(B_2_1_gtypesFil, B_2_1_dataFil$AggloClusLab, center=TRUE, training.set = trainPerc, 
                      result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                      scale=FALSE)
table(B_2_1_Xdapc2$DAPC$assign,B_2_1_Xdapc2$DAPC$grp)

# add up all misclassified samples so far
total_mcSam <- c(total_mcSam,B_2_1_mcSam)


######## B_2_2 as base
B_2_2_data <- read.csv("Plates1to9Test_impRound_NoFrance1panama_w_metadata_B_2_2.csv")
# Exclude root misclassified samples from the downstream steps
B_2_2_data %>% filter(!SampleID %in% total_mcSam) -> B_2_2_dataRem

B_2_2_gtypes <- B_2_2_dataRem[, -c(1:metaCol)]

B_2_2_Xdapc<-xvalDapc(B_2_2_gtypes, B_2_2_dataRem$AggloClusLab, center=TRUE, training.set = trainPerc, 
                     result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                     scale=FALSE)


# To get contingency table
table(B_2_2_Xdapc$DAPC$assign,B_2_2_Xdapc$DAPC$grp)


# Identify samples that got misclassified by xvalDAPC and exclude them from running DAPC
# To get misclassified samples from xvalDapc
missCl <- B_2_2_Xdapc$DAPC$assign != B_2_2_Xdapc$DAPC$grp 
correctCl <- B_2_2_Xdapc$DAPC$assign == B_2_2_Xdapc$DAPC$grp
B_2_2_gtypesFil <- B_2_2_gtypes[correctCl,] # this now retains only samples that got classified as expected
B_2_2_dataFil <- B_2_2_dataRem[correctCl,]
B_2_2_mcSam <- B_2_2_dataRem$SampleID[missCl]


optPcaNum <- as.numeric(B_2_2_Xdapc$`Number of PCs Achieving Lowest MSE`)
B_2_2_dapc <- dapc(B_2_2_gtypesFil, B_2_2_dataFil$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum,
                  scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                  pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)
table(B_2_2_dapc$assign,B_2_2_dapc$grp)
# No misclassified samples from the DAPC run

# re run xvaldapc after excluding misclassified samples
B_2_2_Xdapc2<-xvalDapc(B_2_2_gtypesFil, B_2_2_dataFil$AggloClusLab, center=TRUE, training.set = trainPerc, 
                       result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                       scale=FALSE)
table(B_2_2_Xdapc2$DAPC$assign,B_2_2_Xdapc2$DAPC$grp)

# add up all misclassified samples so far
total_mcSam <- c(total_mcSam,B_2_2_mcSam)




####### C_1_2_2 as base
C_1_2_2_data <- read.csv("Plates1to9Test_impRound_NoFrance1panama_w_metadata_C_1_2_2.csv")
# Exclude total misclassified samples from the downstream steps
C_1_2_2_data %>% filter(!SampleID %in% total_mcSam) -> C_1_2_2_dataRem
C_1_2_2_gtypes <- C_1_2_2_dataRem[, -c(1:metaCol)]


C_1_2_2_Xdapc<-xvalDapc(C_1_2_2_gtypes, C_1_2_2_dataRem$AggloClusLab, center=TRUE, training.set = trainPerc, 
                     result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                     scale=FALSE)


# To get contingency table
table(C_1_2_2_Xdapc$DAPC$assign,C_1_2_2_Xdapc$DAPC$grp)

# No misclassified samples from xvalDapc

optPcaNum <- as.numeric(C_1_2_2_Xdapc$`Number of PCs Achieving Lowest MSE`)
C_1_2_2_dapc <- dapc(C_1_2_2_gtypes, C_1_2_2_dataRem$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum,
                  scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                  pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)
table(C_1_2_2_dapc$assign,C_1_2_2_dapc$grp)
# No misclassified samples from the DAPC run



######## C_2_1_2 as base
C_2_1_2_data <- read.csv("Plates1to9Test_impRound_NoFrance1panama_w_metadata_C_2_1_2.csv")
# Exclude total misclassified samples from the downstream steps
C_2_1_2_data %>% filter(!SampleID %in% total_mcSam) -> C_2_1_2_dataRem
C_2_1_2_gtypes <- C_2_1_2_dataRem[, -c(1:metaCol)]

C_2_1_2_Xdapc<-xvalDapc(C_2_1_2_gtypes, C_2_1_2_dataRem$AggloClusLab, center=TRUE, training.set = trainPerc, 
                     result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                     scale=FALSE)

# To get contingency table
table(C_2_1_2_Xdapc$DAPC$assign,C_2_1_2_Xdapc$DAPC$grp)

# No misclassified samples from xvalDapc

optPcaNum <- as.numeric(C_2_1_2_Xdapc$`Number of PCs Achieving Lowest MSE`)
C_2_1_2_dapc <- dapc(C_2_1_2_gtypes, C_2_1_2_dataRem$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum,
                  scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                  pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)
table(C_2_1_2_dapc$assign,C_2_1_2_dapc$grp)
# No misclassified samples from the DAPC run



######## D_1_2_2_1 as base
D_1_2_2_1_data <- read.csv("Plates1to9Test_impRound_NoFrance1panama_w_metadata_D_1_2_2_1.csv")
# Exclude total misclassified samples from the downstream steps
D_1_2_2_1_data %>% filter(!SampleID %in% total_mcSam) -> D_1_2_2_1_dataRem
D_1_2_2_1_gtypes <- D_1_2_2_1_dataRem[, -c(1:metaCol)]
D_1_2_2_1_Xdapc<-xvalDapc(D_1_2_2_1_gtypes, D_1_2_2_1_dataRem$AggloClusLab, center=TRUE, training.set = trainPerc, 
                     result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                     scale=FALSE)


# To get contingency table
table(D_1_2_2_1_Xdapc$DAPC$assign,D_1_2_2_1_Xdapc$DAPC$grp)


# No misclassified samples from xvalDapc

optPcaNum <- as.numeric(D_1_2_2_1_Xdapc$`Number of PCs Achieving Lowest MSE`)
D_1_2_2_1_dapc <- dapc(D_1_2_2_1_gtypes, D_1_2_2_1_data$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum,
                  scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                  pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)
table(D_1_2_2_1_dapc$assign,D_1_2_2_1_dapc$grp)

# No misclassified samples from the DAPC run

############################################################
######## TIER II ###########################################
############################################################

# for scutellata and capensis


####### create DAPC function for B_1_1 as base
setwd("C:/Users/donth/OneDrive/Kiran_projects/bee_SNP_analysis/genotyping/fluidigm/UIUC/All_plates/KnowBase_Agglom_clus")
B_1_1_data <- read.csv("Plates1to9Test_impRound_NoFrance1panama_w_metadata_B_1_1.csv")

# Exclude root misclassified samples from the B_1_2_data
B_1_1_data %>% filter(!SampleID %in% root_mcSam) -> B_1_1_dataRem

B_1_1_gtypes <- B_1_1_dataRem[, -c(1:metaCol)]
rownames(B_1_1_gtypes) <- B_1_1_dataRem$SampleID
B_1_1_Xdapc<-xvalDapc(B_1_1_gtypes, B_1_1_dataRem$AggloClusLab, center=TRUE, training.set = trainPerc, 
                      result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                      scale=FALSE)


# To get contingency table
table(B_1_1_Xdapc$DAPC$assign,B_1_1_Xdapc$DAPC$grp)

# Identify samples that got misclassified by xvalDAPC and exclude them from running DAPC
# To get misclassified samples from xvalDapc
missCl <- B_1_1_Xdapc$DAPC$assign != B_1_1_Xdapc$DAPC$grp 
correctCl <- B_1_1_Xdapc$DAPC$assign == B_1_1_Xdapc$DAPC$grp
B_1_1_gtypesFil <- B_1_1_gtypes[correctCl,] # this now retains only samples that got classified as expected
B_1_1_dataFil <- B_1_1_dataRem[correctCl,]
B_1_1_mcSam <- B_1_1_dataRem$SampleID[missCl]


optPcaNum <- as.numeric(B_1_1_Xdapc$`Number of PCs Achieving Lowest MSE`)
B_1_1_dapc <- dapc(B_1_1_gtypesFil, B_1_1_dataFil$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum, scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                   pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)

table(B_1_1_dapc$assign,B_1_1_dapc$grp)

# Get misclassified samples from the DAPC run
missCl2 <- B_1_1_dapc$assign != B_1_1_dapc$grp 
correctCl2 <- B_1_1_dapc$assign == B_1_1_dapc$grp
B_1_1_gtypesFil2 <- B_1_1_gtypesFil[correctCl2,]
B_1_1_gtypesFil2 <- na.omit(B_1_1_gtypesFil2)
B_1_1_dataFil2 <- B_1_1_dataFil[correctCl2,]
B_1_1_mcSam <- c(B_1_1_mcSam, B_1_1_dataFil$SampleID[missCl2])

# re run xvaldapc and dapc after excluding misclassified samples
B_1_1_Xdapc2<-xvalDapc(B_1_1_gtypesFil2, B_1_1_dataFil2$AggloClusLab, center=TRUE, training.set = trainPerc, 
                      result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                      scale=FALSE)
table(B_1_1_Xdapc2$DAPC$assign,B_1_1_Xdapc2$DAPC$grp)
# Get misclassified samples from the above XDAPC run
missCl <- B_1_1_Xdapc2$DAPC$assign != B_1_1_Xdapc2$DAPC$grp 
correctCl <- B_1_1_Xdapc2$DAPC$assign == B_1_1_Xdapc2$DAPC$grp
B_1_1_gtypesFil3 <- B_1_1_gtypesFil2[correctCl,] # this now retains only samples that got classified as expected
B_1_1_dataFil3 <- B_1_1_dataFil2[correctCl,]
B_1_1_mcSam <- c(B_1_1_mcSam, B_1_1_dataFil3$SampleID[missCl])


B_1_1_dapc2 <- dapc(B_1_1_gtypesFil3, B_1_1_dataFil3$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum,
                   scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                   pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)
table(B_1_1_dapc2$assign,B_1_1_dapc2$grp)




####### create DAPC function for C_1_2_1 as base
setwd("C:/Users/donth/OneDrive/Kiran_projects/bee_SNP_analysis/genotyping/fluidigm/UIUC/All_plates/KnowBase_Agglom_clus")
C_1_2_1_data <- read.csv("Plates1to9Test_impRound_NoFrance1panama_w_metadata_C_1_2_1.csv")

# Exclude root misclassified samples from the C_1_2_1_data
C_1_2_1_data %>% filter(!SampleID %in% root_mcSam) -> C_1_2_1_dataRem

C_1_2_1_gtypes <- C_1_2_1_dataRem[, -c(1:metaCol)]
rownames(C_1_2_1_gtypes) <- C_1_2_1_dataRem$SampleID
C_1_2_1_Xdapc<-xvalDapc(C_1_2_1_gtypes, C_1_2_1_dataRem$AggloClusLab, center=TRUE, training.set = trainPerc, 
                      result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                      scale=FALSE)


# To get contingency table
table(C_1_2_1_Xdapc$DAPC$assign,C_1_2_1_Xdapc$DAPC$grp)

# No misclassified samples from xvalDapc

optPcaNum <- as.numeric(C_1_2_1_Xdapc$`Number of PCs Achieving Lowest MSE`)
C_1_2_1_dapc <- dapc(C_1_2_1_gtypes, C_1_2_1_dataRem$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum, scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                   pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)

table(C_1_2_1_dapc$assign,C_1_2_1_dapc$grp)

# No misclassified samples from the DAPC run


####### create DAPC function for D_1_2_2_2 as base
D_1_2_2_2_data <- read.csv("Plates1to9Test_impRound_NoFrance1panama_w_metadata_D_1_2_2_2.csv")

# Exclude root misclassified samples from the C_1_2_1_data
D_1_2_2_2_data %>% filter(!SampleID %in% total_mcSam) -> D_1_2_2_2_dataRem

D_1_2_2_2_gtypes <- D_1_2_2_2_dataRem[, -c(1:metaCol)]
rownames(D_1_2_2_2_gtypes) <- D_1_2_2_2_dataRem$SampleID
D_1_2_2_2_Xdapc<-xvalDapc(D_1_2_2_2_gtypes, D_1_2_2_2_dataRem$AggloClusLab, center=TRUE, training.set = trainPerc, 
                        result = c("groupMean", "overall"), xval.plot = TRUE, n.pca = NULL,
                        scale=FALSE)


# To get contingency table
table(D_1_2_2_2_Xdapc$DAPC$assign,D_1_2_2_2_Xdapc$DAPC$grp)

# No misclassified samples from xvalDapc

optPcaNum <- as.numeric(D_1_2_2_2_Xdapc$`Number of PCs Achieving Lowest MSE`)
D_1_2_2_2_dapc <- dapc(D_1_2_2_2_gtypes, D_1_2_2_2_dataRem$AggloClusLab, center=TRUE, n.pca = optPcaNum, n.da = optPcaNum, scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                     pca.select=c("nbEig","percVar"), perc.pca=NULL, dudi=NULL)

table(D_1_2_2_2_dapc$assign,D_1_2_2_2_dapc$grp)

# No misclassified samples from the DAPC run




#save.image("dapc_funcs_AllNodes.RData")
load("dapc_funcs_AllNodes.RData")

