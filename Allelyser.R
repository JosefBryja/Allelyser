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
  SNP <<- cbind(data[, 1], data[, a])
  print(paste(colnames(SNP)[2], "selected"))
  
  # Cleaning the SNP table
  SNP <<- SNP[SNP[,2] %in% c("AA", "GG", "AG","TT", "CC", "CT", NA), ]
  
  # Factorizing the data
  SNP[, 1] <<- as.factor(SNP[, 1])
  SNP[, 2] <<- as.factor(SNP[, 2])
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
    
  # Print out the output table
  if(significant_only == TRUE){
    print(na.omit(pData[pData[,2] < alpha,]))
  }else{
    print(pData)
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


##########################
###### Presentation ######
##########################

ChooseFile(here("APCgenotypesAnonym.xlsx"))
HWE(SNP)
ChiSq(SNP, significant_only = TRUE)
SNPHeatmap(SNP, scaled = TRUE)





