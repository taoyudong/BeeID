## This script plots genotype dosage for 
## selected terminal nodes in the AHB SNP Panel

#----------------------------------------------------------------------------
# Load analysis options and libraries
#----------------------------------------------------------------------------

# libraries
pkgs <- c("tidyverse", "magrittr", "RColorBrewer", "data.table", "MASS", "car", "BiocParallel")
invisible(lapply(pkgs, library, character.only = T))
rm(pkgs); gc()

# options
options(stringsAsFactors = F, scipen = 9999)
set.seed(12345)

#----------------------------------------------------------------------------
# Data Structuring
#----------------------------------------------------------------------------

# read in the files
gt <- fread("./data_sup_tab4.csv", sep = ",", data.table = F)

# create an apply loop that consolidates proportion of genotypes
nt.prop <- tapply(X     = 1:nrow(gt),
                  INDEX = as.factor(gt$KBN_code),
                  FUN   = function(i, data){
                  tmp   = data[i, ]
                  prop  = data.frame(REF = colSums(tmp == 0) / nrow(tmp),
                                     HET = colSums(tmp == 1) / nrow(tmp),
                                     ALT = colSums(tmp == 2) / nrow(tmp))
                  },
                  data  = as.matrix(gt[, -c(1:10)])
                 )

# plot each node
pdf(file = "./terminal_nodes.pdf", width = 6.5, height = 2)
lapply(1:length(nt.prop), function(i, data){
  print(i)
  barplot(t(data[[i]]), col = c("green4", "grey80", "magenta3"),
          border = NA, axes = F, xaxt = "n", ylab = names(data)[i])
}, data = nt.prop)
dev.off()

# repeat this process for the head nodes
# create the index
h.idx <- setNames(
    nm     = unique(gt$KBN_code),
    object = c("A_2", "A_1", "A_2", "A_2", "A_1", "A_1", "A_1", "A_1", "A_2", "A_2")
)

#apply the same summary to the head nodes
h.prop <- tapply(X     = 1:nrow(gt),
                 INDEX = as.factor(as.character(h.idx[gt$KBN_code])),
                 FUN   = function(i, data){
                 tmp   = data[i, ]
                 prop  = data.frame(REF = colSums(tmp == 0) / nrow(tmp),
                                    HET = colSums(tmp == 1) / nrow(tmp),
                                    ALT = colSums(tmp == 2) / nrow(tmp))
                 },
                 data  = as.matrix(gt[, -c(1:10)])
                )

# plot the head nodes
pdf(file = "./head_nodes.pdf", width = 6.5, height = 2)
lapply(1:length(h.prop), function(i, data){
  print(i)
  barplot(t(data[[i]]), col = c("green4", "grey80", "magenta3"),
          border = NA, axes = F, xaxt = "n", ylab = names(data)[i])
}, data = h.prop)
dev.off()

