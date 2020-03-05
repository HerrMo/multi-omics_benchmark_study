### merge the single benchmarking results

# 1. Loading ancillary code ---------------------------------------------------
library(data.table)
library(mlr)
library(dplyr)

source("R/ancillary_code_merge.R")

# 2. merging bench results to df ----------------------------------------------

count_bmr <- 234
bmr_res <- vector(mode = "list", length = count_bmr)

for (bmr in 1:count_bmr) {
  bmr_res[[bmr]] <- readRDS(paste0("data/bmr-results/", bmr, ".rds"))
}

erg2 <- bmr_res

# Attention: 
# bmr_res needs little adjustment for mergeBenchmarkResults to work
#  1. task description of grridge/BRCA flawed, using other BRCA desc to adjust
#     erg2[[52]]$results$BRCA$grridge$task.desc <- erg2[[43]]$results[[1]][[1]]$task.desc
#  2. resample description of grridge/UCEC flawed, using other UCEC desc to adjust
#     erg2[[234]]$results$UCEC$grridge$pred$instance$desc <- erg2[[52]]$results$BRCA$grridge$pred$instance$desc
# nicenames has to be adjusted and plots and tables accordingly

erg2[[52]]$results$BRCA$grridge$task.desc <- erg2[[43]]$results[[1]][[1]]$task.desc
erg2[[234]]$results$UCEC$grridge$pred$instance$desc <- erg2[[52]]$results$BRCA$grridge$pred$instance$desc

# check for fully observed data sets
res_names <- unlist(lapply(erg2, function(x) names(x[[1]])))
res <- as.vector(table(res_names))
names(res) <- sort(unique(res_names))
datsets <- names(res[res == 13]) # fully observed data sets

# make list of merged benchmark results for every data set
bmr_lists <- vector("list", length(datsets))
names(bmr_lists) <- datsets
for (task_id in datsets) {
  lis <- list()
  for (i in 1:234) {
    if (names(erg2[[i]][[1]][1]) == task_id) {
      lis = c(lis, list(erg2[[i]]))
    }
  }
  bmr_lists[[task_id]] <- mergeBenchmarkResults(lis)
}

# impute missings values
bmr_df_lists <- lapply(bmr_lists, as.data.frame)
bmr_df_lists <- lapply(bmr_df_lists, impute_missings2, na_rat = 0.2)

df_res <- do.call(rbind.data.frame, bmr_df_lists)
df_res <- as.data.table(df_res)
df_res[learner.id == "Kaplanmeier" & cindex.uno == 0, "cindex.uno"] <- 0.5 # set cindex to 0.5 for kaplanmeier, since mlr doesn't do this
df_res[learner.id == "CoxBoost" & cindex.uno == 0, "cindex.uno"] <- 0.5 # for mstop = 0, set cindex to 0.5

nicenam <- c("Kaplan-Meier", "Lasso", "glmboost", "CoxBoost", "Clinical only",
             "prioritylasso", "prioritylasso favoring", "CoxBoost favoring", 
             "grridge", "blockForest", "rfsrc", "ranger", "ipflasso")

# short_nicenam <- c("km", "lasso", "glmb", "coxb", "clin", "priorl", "priorl_f", 
#                    "coxb_f", "grrdige", "blockf", "rfsrc", "ranger", "ipfl")

levels(df_res$learner.id) <- nicenam
head(df_res)

# save(df_res, file = "data/merged-results.RData")
