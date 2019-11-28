# 1. packages via checkpoint snapshot
# the following packages will be loaded by checkpoint (all library instances)
library(pec) 
library(glmnet) 
library(prioritylasso) 
library(CoxBoost) 
library(ipflasso)
library(randomForestSRC) 
library(ranger)
library(dplyr)
library(data.table) 
library(survival)
library(GRridge) 
library(batchtools)
library(devtools)
library(SGL)
library(mboost)
library(blockForest)
library(tuneRanger)
library(survAUC)

# 2. packages not covered by checkpoint
# some of the packages could not be snapshoted via checkpoint, so we cloned them

# !!!
# mlr branch with advanced survial methods not included in the CRAN version
# !!!
devtools::install_github("HerrMo/mlr@survival_probabilities_ibrier",
                         dependencies = FALSE) 
devtools::install_github("HerrMo/GRridge") # is not coverd by checkpoint since bioconductor package
install.packages(ipflasso) # version with two-step ipf-lasso was not available on CRAN at the time of experiment

# additional clones
# devtools::install_github("HerrMo/blockForest") 
# devtools::install_github("HerrMo/tuneRanger")
# devtools::install_github("HerrMo/mboost") 

# If problems with devtools occure (under windows)
# library(devtools)             
#
# assignInNamespace("version_info",                                              
#                   c(devtools:::version_info, 
#                     list("3.5" = list(version_min = "3.3.0", 
#                                       version_max = "99.99.99",
#                                       path = "bin"))), 
#                   "devtools")
#                                                                                
# install_github("mlr-org/mlr@survival_probabilities_ibrier",                    
#                dependencies = FALSE)