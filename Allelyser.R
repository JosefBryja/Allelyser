#####################################################
##        ___    ____     __                       ##
##       /   |  / / /__  / /_  __________  _____   ##
##      / /| | / / / _ \/ / / / / ___/ _ \/ ___/   ##
##     / ___ |/ / /  __/ / /_/ (__  )  __/ /       ##
##    /_/  |_/_/_/\___/_/\__, /____/\___/_/        ##
##                      /____/                     ##
#####################################################

######################
###### Packages ######
######################

library(tidyverse)
library(readxl)
library(here)


##########################
###### Data loading ######
##########################

ChooseFile <- function(path){
  # Reading the excel file
  data <- read_xlsx(path)
  # Selecting a SNP
  a <- menu(colnames(data[,-1]), title="Please, select the SNP of your interest:") + 1
  SNP <- cbind(data[, 1], data[, a])
  print(paste(colnames(SNP)[2], "selected"))
  
  # Cleaning the SNP table
  SNP <- SNP[SNP[,2] %in% c("AA", "GG", "CC","TT", "AT", "AC", "AG", "CT", "CG", "GT", NA), ]
  
  # Factorizing the data
  SNP[, 1] <- as.factor(SNP[, 1])
  SNP[, 2] <- as.factor(SNP[, 2])
  
  # Return the SNP table
  return(SNP)
}


####################################################
###### Testing for Hardy-Weinberg equilibrium ######
####################################################

HWE <- function(SNPtable, alpha = 0.05){
  # Observed genotype frequencies in the dataset
  HWTable <- table(SNPtable[, 2])

    # Observed allele frequencies in the dataset
  HWsum <- 2*(HWTable[[1]] + HWTable[[2]] + HWTable[[3]])
  p <- (2*HWTable[[1]] + HWTable[[2]])/HWsum
  q <- (2*HWTable[[3]] + HWTable[[2]])/HWsum
  
  # Theoretical probabilities of genotype frequencies
  HWexp <- c(p^2, 2*p*q, q^2)
  
  # Adding the expected genotype frequencies to the table
  HWTable <- rbind(HWTable, HWexp*sum(HWTable))
  rownames(HWTable) <- c("Observed", "Expected")
  
  # Printing the table with expected genotype counts
  print(HWTable)
  
  # Finally, we can perform the chi squared test
  # We use just observed values for the testing. To determine the theoretical allele frequencies,
  # we use the 'p' argument of the chisq.test() function.
  pvalue <- chisq.test(HWTable[1,], p = HWexp)$p.value
  
  # Output of the function
  if(pvalue < alpha){
    print(paste0("The population is not at Hardy-Weinberg equilibrium (p-value = ", round(pvalue, 4), ")"))
  }else{
    print(paste0("The population is at Hardy-Weinberg equilibrium (p-value = ", round(pvalue, 4), ")"))
  }
  return(pvalue)
}


##############################################################################
###### Testing if there is association among the phenotype and genotype ######
##############################################################################

ChiSq <- function(data, significant_only = FALSE, alpha = 0.05){
  
  # Creating and printing contingency table
  conTable <- table(data)
  
  # For each diagnosis, we want to test if there is difference in the genotype frequencies
  # Empty data frame
  pData <- data.frame(matrix(nrow = 0, ncol = 2))

  # Chi-square for each diagnosis
  for(x in 1:length(rownames(conTable))){
    
    # Single diagnosis
    t <- conTable[x, ]
    
    # Chi-square test
    if(sum(t) > 0){
      pData <- rbind(pData, c(rownames(conTable)[x], round(chisq.test(t)$p.value, 6)))
    }else{
      pData <- rbind(pData, c(rownames(conTable)[x], NA))
    }
  }
  
  # Renaming colums
  colnames(pData) <- c("Diagnosis", "p-value")
  
  # Changing p-values to numeric
  pData[,2] <- as.numeric(pData[,2])
    
  # Return the output table
  if(significant_only == TRUE){
    return(na.omit(pData[pData[,2] < alpha,]))
  }else{
    return(pData)
  }
}


###########################
###### Visualization ######
###########################

SNPHeatmap <- function(data, scaled = TRUE){
  
  # Create data for visualization
  name <- colnames(data)[2]
  t <- table(data)
  t <- t[rowSums(t) > 0,]
  
  # Plot the data
  if(scaled == TRUE){
    heatmap(t, scale = "row", main = name)
  }else if(scaled == FALSE){
    heatmap(t, scale = "none")
  }
}

######################################
###### Analysis of multiple SNP ######
######################################

# Load the data
SNPdata <- read_xlsx(here("APCgenotypesAnonym.xlsx"))
data_cl <- as.data.frame(lapply(SNPdata[,-1], function(x) ifelse(!x %in% c("AA", "GG", "CC","TT", "AT", "AC", "AG", "CT", "CG", "GT", NA), NA, x)))
data_cl <- cbind(SNPdata[, 1], data_cl)

# Extract rsIDs
SNP_ID <- str_extract(colnames(SNPdata)[-1], "(?<=_|\\b)rs\\d+")

# Chisq test for all the SNPs
ChiSqAll <- function(data){
  chi_sq <- data.frame(Diagnosis = levels(as.factor(data[,1])))
  
  # Single SNP
  for(i in 2:ncol(data)){
    ct <- table(data[,1], data[,i])
    part_chi <- data.frame(matrix(nrow = 0, ncol = 2))
    
    for(x in 1:nrow(ct)){
      # Single diagnosis
      t <- ct[x, ]
      
      # Chi-square test
      if(sum(t) > 0){
        part_chi <- rbind(part_chi, c(rownames(ct)[x], round(chisq.test(t)$p.value, 6)))
      }else{
        part_chi <- rbind(part_chi, c(rownames(ct)[x], NA))
      }
    }
    colnames(part_chi) <- c("Diagnosis", colnames(data_cl)[i])
    chi_sq <- merge(chi_sq, part_chi, by = "Diagnosis")
  }
  return(chi_sq)
}


#####################
###### LD plot ######
#####################

calculate_ld <- function(snp1, snp2) {
  snp1_factor <- factor(snp1)
  snp2_factor <- factor(snp2)
  contingency_table <- table(snp1_factor, snp2_factor)
  chi_square <- chisq.test(contingency_table)$statistic
  sqrt(chi_square / length(snp1))
}

ld_matrix <- outer(data_cl[,-1], data_cl[,-1], Vectorize(calculate_ld))
ld_matrix[lower.tri(ld_matrix)] <- NA
rownames(ld_matrix) <- SNP_ID
colnames(ld_matrix) <- SNP_ID

ggplot(data = as.data.frame(as.table(ld_matrix)), aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "cornflowerblue", high = "darkred", na.value = "transparent", guide = "legend") +  
  theme_minimal() +
  labs(x = "", y = "", title = "LD Plot") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  coord_fixed() +
  theme(legend.position = "right")


##########################
###### Presentation ######
##########################

# Single polymorphisms
SNP <- ChooseFile(here("APCgenotypesAnonym.xlsx"))
HWE(SNP)
table(SNP)
ChiSq(SNP)
SNPHeatmap(SNP, scaled = TRUE)

Chi_sq_table <- ChiSqAll(data_cl)


