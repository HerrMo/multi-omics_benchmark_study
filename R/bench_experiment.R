### 
### Benchmarking experiment ###
###

# Pipeline: OpenML for datasets, mlr for benchmark, batchtools for parallelisat.

# TODO
# !!! Take care of packages version, especially mlr (s. packages)!!!
# !!! Take care of makeClusterFunction !!!

# Setup ---------------------------------------------------

library(checkpoint) # get package snapshots used packages
checkpoint("2019-04-17", project = getwd())

library(batchtools)

packs <- c("mlr", "pec", "glmnet", "mboost", "prioritylasso", "CoxBoost", 
           "ipflasso", "randomForestSRC", "ranger", "dplyr", "data.table", "survival",
           "GRridge", "tuneRanger", "SGL", "blockForest", "survAUC", "riskRegression")

# unlink("bench_exp", recursive = TRUE)
regis <- makeExperimentRegistry("bench_exp", 
                                packages = packs,
                                  source = "ancillary_code_bench.R")


# 0. Partial reproduction -----------------------------------------------------

# If only parts of the benchmark experiment are to be run, e.g. for reproducing 
# single results due to computation times, change settings as desired 
# (here full bench experiment). 
# Important: Change seed to randomly choose the datasets/learners for which
# the results should be reproduced.

# n_lrns <- 14  
# n_datsets <- 18
# seed <- 124
# 
# partition <- sample_reproduction(n_lrns = n_lrns, 
#                                  n_datsets = n_datsets, 
#                                  seed = seed)

# 1. Tasks and problems--------------------------------------------------------

# OpenML dataset ids to querry

load("data/datset_ids.RData")

nams <- c("LAML", "BLCA", "LGG",  "BRCA", "COAD", "ESCA", 
          "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", 
          "OV", "PAAD", "SARC", "SKCM", "STAD", "UCEC")

# For partial reproduction:
# nams <- nams[partition$datsets] 

for (nam in nams) {

  # download dataset
  dat_part1 <- getOMLDataSet(datset_ids[[nam]][[1]])
  dat_part2 <- getOMLDataSet(datset_ids[[nam]][[2]])

  dat <- cbind.data.frame(dat_part1, dat_part2)
  
  if (nam = "BRCA") {
    dat_part3 <- getOMLDataSet(datset_ids[[nam]][[3]])
    dat <- cbind.data.frame(dat, dat_part3)
  }
  
  task <- convertOMLDataSetToMlr(obj = dat,
                                 task.type = "Survival analysis",
                                 target = c("time", "status"),
                                 fix.colnames = FALSE)
  # convert to mlr task
  task <- makeSurvTask(id = nam, 
                       data = task[, -1], # delet patient code
                       target = c("time", "status"))

  # adding task as batch problem
  addProblem(name = getTaskId(task), data = task)
}


# ------------------
# get data from disk

# for (nam in nams) {
#   load("data/", nam, ".RData"))
#   task <- get(nam)[, -1]
#   task <- makeSurvTask(id = nam, data = task, target = c("time", "status"))
#   addProblem(name = getTaskId(task), data = task)
# }
# ------------------

rm(list = c(nams)) # the loaded datasets occupy a lot memory, so delete


# 2. Algorithms and learners --------------------------------------------------

# For each learner an algorithm is defined to make paralleliz. on learner-level
# possible.
# Generic learners do not depend on dataset specific parameters (group structure)
# Task specific learners depend on dataset specific parameters (group structure)
#
# In general, lrnrs use default settings. Find other lrnr-configurations
# in the following list. 
# Task specific arguments cannot be set here and are added in make_spec_lrns(). 
# They are listed here for convenience overview and set to NULL.

# defining lrn params in advance and using for loop to DRY
l_lrn_args = list("lrn_km" = list(cl = "surv.kaplanmeier",
                                  id = "Kaplanmeier",
                                  predict.type = "prob"),
                  "lrn_lasso" = list(cl = "surv.cv.glmnet2",
                                     id = "Lasso",
                                     s = "lambda.min",
                                     predict.type = "prob"),
                  "lrn_glmboost" = list(cl = "surv.cv.glmboost",
                                        id = "glmBoost",
                                        use.formula = FALSE,
                                        mstop = 1L,
                                        predict.type = "prob"),
                  "lrn_CoxBoost" = list(cl = "surv.cv.CoxBoost2",
                                        id = "CoxBoost",
                                        predict.type = "prob"),
                  "lrn_rfsrc" = list(cl = "surv.randomForestSRC",
                                     id = "rfsrc",
                                     predict.type = "prob"),
                  "lrn_ranger" = list(cl = "surv.tuneMtryFast2",
                                      id = "ranger",
                                      predict.type = "prob",
                                      write.forest = TRUE),
                  "lrn_clin_ref" = list(cl = "surv.clinic_reference",
                                        id = "Clinical only",
                                        predict.type = "prob",
                                        clinicals = NULL,
                                        nfolds = 10), # only used if p_clin > n_train (which is not the case)
                  "lrn_ts_prior" = list(cl = "surv.ts.priorlasso",
                                        id = "Prioritylasso",
                                        blocks = NULL,
                                        predict.type = "prob",
                                        favoring = FALSE),
                  "lrn_ts_prior_fav" = list(cl = "surv.ts.priorlasso",
                                            id = "Prioritylasso favoring",
                                            blocks = NULL,
                                            predict.type = "prob",
                                            favoring = TRUE),
                  "lrn_tsipf" = list(cl = "surv.ts.ipflasso",
                                     id = "IPF-Lasso",
                                     blocks = NULL,
                                     nfolds = 10,
                                     ncv = 1,
                                     predict.type = "prob"),
                  "lrn_cv_coxboost_unpen" = list(cl = "surv.cv.CoxBoost2",
                                                 id = "Coxboost favoring",
                                                 predict.type = "prob",
                                                 unpen.index = NULL),
                  "lrn_blockForest" = list(cl = "surv.blockForest",
                                           id = "blockForest",
                                           predict.type = "prob",
                                           blocks = NULL),
                  "lrn_SGL" = list(cl = "surv.cvSGL",
                                   id = "SGL",
                                   predict.type = "prob",
                                   index = NULL,
                                   nfold = 10),
                  "lrn_grridge" = list(cl = "surv.grridge",
                                       id = "grridge",
                                       partitions = NULL,
                                       predict.type = "prob",
                                       innfold = 10,
                                       standardizeX = TRUE,
                                       selectionEN = TRUE,
                                       maxsel = c(1000))
)

# You can also choose specific lrns (random selection will only take place if
# no manual selection is done)

# all_lrns <- names(l_lrn_args) # default is not to select manually (i.e., all lrns)
# l_lrn_args <- l_lrn_args[all_lrns]

# if (length(l_lrn_args) == 14) {
#   l_lrn_args <- l_lrn_args[partition$datsets] # check if $datsets correct argument
# }


addAlgorithm("learner", fun = function(job, data, instance, lrns, ...) {
  par.vals = list(...)
  
  lrnr = lrns
  
  # 1. Task -----------------------------------------------------------------   all equal
  task = data
  task_id = getTaskId(task)
  
  # 2. Learner --------------------------------------------------------------   main diffs here
  if (lrnr %in% c("lrn_clin_ref", "lrn_ts_prior", "lrn_ts_prior_fav", 
                  "lrn_tsipf","lrn_cv_coxboost_unpen", "lrn_blockForest",
                  "lrn_SGL", "lrn_grridge")) {
    lrn = make_spec_lrns(task, args = l_lrn_args, lrnr = lrnr)
  } else {
    lrn = list(assign(lrnr, do.call(makeLearner, l_lrn_args[[lrnr]])))
  }
  
  
  # 3. Measures -------------------------------------------------------------   all equal
  mrs = list(timetrain, cindex.uno, ibrier, featselc_default, featselc_clin, 
             featselc_cnv, featselc_mirna, featselc_mutation, featselc_rna) 
  
  # 4. Resampling -----------------------------------------------------------   all equal
  if (task_id %in% c("BRCA", "LUAD", "LUSC", "HNSC", "LGG", "UCEC", "BLCA")) {
    rdesc = makeResampleDesc("RepCV", reps = 5, folds = 5, stratify = TRUE)
  } else {
    rdesc = makeResampleDesc("RepCV", reps = 10, folds = 5, stratify = TRUE)    
  }
  
  set.seed(124) 
  rin = makeResampleInstance(rdesc, task = task)
  if (lrnr == "lrn_cv_coxboost_unpen") {
    save(rin, file = paste0(task_id, "_rin.RData"))
  }
  
  # 5. Benchmarking ---------------------------------------------------------   all equal 
  configureMlr(on.learner.error = "warn", show.learner.output = FALSE)
  set.seed(421)
  
  # keep.pred must be TRUE to be able to merge the results later
  bmr = benchmark(lrn, task, rin, mrs, 
                  keep.pred = TRUE, models = FALSE, show.info = TRUE)
  bmr
})


algo.designs <- list(learner = data.table(lrns = c(names(l_lrn_args))))

addExperiments(algo.designs = algo.designs)
summarizeExperiments()

# kernels to use (depends on system, here linux)
regis$cluster.functions <- makeClusterFunctionsMulticore(6)                    # change according to system (?makeClusterFunctionsMulticore)


# group jobs by learners with fast, normal, and slow run time and submit sequentially
# fast <- findExperiments(
#   algo.name = "learner",
#   algo.pars = (
#     lrns == "lrn_km" |
#       lrns == "lrn_ranger" |
#       lrns == "lrn_clin_ref"
#   )
# )

# normal <-
#   findExperiments(
#     algo.name = "learner",
#     algo.pars = (
#       lrns != "lrn_blockForest" &
#         lrns != "lrn_grridge" & 
#         lrns != "lrn_SGL" &
#         lrns != "lrn_km" &
#         lrns != "lrn_ranger" &
#         lrns != "lrn_clin_ref"
#     )
#   )

# slow <-
#   findExperiments(
#     algo.name = "learner",
#     algo.pars = (
#       lrns == "lrn_blockForest" |
#       lrns == "lrn_grridge"
#     )
# )

submitJobs()
getStatus()

done = findDone()

erg <- reduceResultsList(findExperiments())
errors = getErrorMessages()

save(erg, errors, file = "ergebnis.RData")
