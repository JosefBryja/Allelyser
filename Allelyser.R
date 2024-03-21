#####################################################
##        ___    ____     __                       ##
##       /   |  / / /__  / /_  __________  _____   ##
##      / /| | / / / _ \/ / / / / ___/ _ \/ ___/   ##
##     / ___ |/ / /  __/ / /_/ (__  )  __/ /       ##
##    /_/  |_/_/_/\___/_/\__, /____/\___/_/        ##
##                      /____/                     ##
#####################################################

# Packages
library(tidyverse)
library(readxl)
library(here)

##########################
###### Data loading ######
##########################

# Reading the excel file
data <- read_xlsx(here("APCgenotypesAnonym.xlsx"))

# Selecting a SNP
a <- menu(colnames(data[,-1]), title="Please, select the SNP of your interest:") + 1
SNP <- cbind(data[, 1], data[, a])
print(paste(colnames(SNP)[2], "selected"))

# Cleaning of the SNP table
SNP <- SNP[SNP[,2] %in% c("AA", "GG", "AG","TT", "CC", "CT", NA), ]

# Factorizing the data
SNP[, 1] <- as.factor(SNP[, 1])
SNP[, 2] <- as.factor(SNP[, 2])


####################################################
###### Testing for Hardy-Weinberg equilibrium ######
####################################################

# Observed genotype frequencies in the dataset
HWTable <- table(SNP[, 2])

# Observed allele frequencies in the dataset
p <- (2*HWTable[[1]] + HWTable[[2]])/(2*(HWTable[[1]] + HWTable[[2]] + HWTable[[3]]))
q <- (2*HWTable[[3]] + HWTable[[2]])/(2*(HWTable[[1]] + HWTable[[2]] + HWTable[[3]]))

# Theoretical probabilities of genotype frequencies
HW_exp <- c(p^2, 2*p*q, q^2)

# Adding the expected genotype frequencies to the table
HWTable <- rbind(HWTable, HW_exp*sum(HWTable))
rownames(HWTable) <- c("Observed", "Expected")

# Printing the table with expected genotype counts
print(HWTable)

# Finally, we can perform the chi squared test
# We use just observed values for the testing. To determine the theoretical allele frequencies,
# we use the 'p' argument of the chisq.test() function.
chisq.test(HWTable[1,], p = HW_exp)


##############################################################################
###### Testing if there is association among the phenotype and genotype ######
##############################################################################

# Creating and printing contingency table
conTable <- table(SNP)






