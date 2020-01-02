### producing tables and figures 
### ancillary code

# 1. impute missing values -------------------------------------------------------

impute_missings2 <- function(df_bmr, na_rat) {
  # max_errs <- apply(df_bmr[,-c(1:3, 5)], 2, max, na.rm = TRUE)
  # min_cindex <- min(df_bmr[df_bmr$learner.id != "Kaplanmeier", "cindex.uno"], na.rm = TRUE) 
  
  cindex <- 0.5
  ibrier <- 0.25
  
  errs <- c(cindex, ibrier)
  names(errs) <- c("cindex.uno", "ibrier")
  
  meas <- names(df_bmr[,-c(1:3)])
  lrns <- unique(df_bmr$learner.id)
  
  df_out <- df_bmr
  
  for (lrn in lrns) {
    df_tmp <- df_bmr[df_bmr$learner.id == lrn, ]
    if (any(is.na(df_tmp))) {
      n_fails <- sum(apply(df_tmp, 1, FUN = function(x) any(is.na(x))))
      r_fails <- n_fails/nrow(df_tmp)
      if (r_fails > na_rat) {
        for (mea in meas) {
          if (mea %in% c("cindex.uno", "ibrier")) {
            df_tmp[is.na(df_tmp[, mea]), mea] <- errs[names(errs) == mea]
          } else {
            mean_err <- mean(df_tmp[, mea], na.rm = TRUE)
            if (mea != "timetrain") {
              df_tmp[is.na(df_tmp[, mea]), mea] <- round(mean_err, 0)
            } else {
              df_tmp[is.na(df_tmp[, mea]), mea] <- mean_err
            }
          }
        }
      } else {
        mean_errs <- apply(df_tmp[,-c(1:3)], 2, mean, na.rm = TRUE)
        names(mean_errs) <- names(df_tmp[,-c(1:3)])
        for (mea2 in names(mean_errs)) {
          if (!mea2 %in% c("cindex.uno", "ibrier", "timetrain")) {
            df_tmp[is.na(df_tmp[,mea2]), mea2] <- round(mean_errs[names(mean_errs) == mea2], 2)
          } else {
            df_tmp[is.na(df_tmp[,mea2]), mea2] <- mean_errs[names(mean_errs) == mea2]
          }
        }
      }
      df_out[df_bmr$learner.id == lrn, ] <- df_tmp
    }
  }
  df_out
}

# 2. plot functions --------------------------------------------------------------

my.ranks <- function (df, bmr, measure = NULL, ties.method = "average", aggregation = "default") {
  checkmate::assertClass(bmr, "BenchmarkResult")
  measure = mlr:::checkBMRMeasure(measure, bmr)
  checkmate::assertChoice(aggregation, c("mean", "default"))
  if (aggregation == "mean") {
    dt = setDT(df)
    dt = dt[, list(x = mean(get(measure$id))), by = c("task.id", 
                                                      "learner.id")]
  }
  else if (aggregation == "default") {
    aggr.meas = measureAggrName(measure)
    df = setDT(getBMRAggrPerformances(bmr, as.df = TRUE))
    df = df[, c("task.id", "learner.id", aggr.meas), with = FALSE]
    setnames(df, aggr.meas, "x")
  }
  if (!measure$minimize) 
    dt$x = -dt$x
  dt[, `:=`("alg.rank", rank(.SD$x, ties.method = ties.method)), 
     by = "task.id"]
  dt = melt(setDF(dt), c("task.id", "learner.id"), "alg.rank")
  dt = dcast(dt, learner.id ~ task.id)
  task.id.names = setdiff(colnames(dt), "learner.id")
  mat = as.matrix(dt[, task.id.names])
  rownames(mat) = dt$learner.id
  colnames(mat) = task.id.names
  return(mat)
}

my.plot.ranks <- function (df, bmr, measure = NULL, ties.method = "average", aggregation = "default", 
                           pos = "stack", order.lrns = NULL, order.tsks = NULL, pretty.names = TRUE) {
  checkmate::assertClass(bmr, "BenchmarkResult")
  measure = mlr:::checkBMRMeasure(measure, bmr)
  checkmate::assertChoice(pos, c("tile", "stack", "dodge"))
  df = as.data.frame(my.ranks(df, bmr, measure, ties.method = ties.method, 
                              aggregation = aggregation))
  df$learner.id = factor(rownames(df))
  setDT(df)
  df = melt(df, id.vars = "learner.id")
  setnames(df, c("variable", "value"), c("task.id", "rank"))
  df = mlr:::orderBMRLrns(bmr, df, order.lrns)
  df = mlr:::orderBMRTasks(bmr, df, order.tsks)
  if (pretty.names) {
    learner.ids = getBMRLearnerIds(bmr)
    learner.short.names = getBMRLearnerShortNames(bmr)
    checkDuplicatedLearnerNames(learner.short.names)
    names(learner.short.names) = learner.ids
    if (is.null(order.lrns)) {
      learner.short.names = learner.short.names[sort(learner.ids)]
    }
    else {
      learner.short.names = learner.short.names[order.lrns]
    }
    levels(df$learner.id) = learner.short.names
  }
  df$rank = as.factor(df$rank)
  setDF(df)
  if (pos == "tile") {
    p = ggplot(df, aes_string("rank", "task.id", fill = "learner.id"))
    p = p + geom_tile()
    p = p + ylab(NULL)
  }
  else {
    p = ggplot(df, aes_string("rank", fill = "learner.id"))
    p = p + geom_bar(position = pos)
    p = p + ylab(NULL)
  }
  return(p)
}

my.boxplots <- function (df, bmr, measure = NULL, style = "box", order.lrns = NULL, 
                         order.tsks = NULL, pretty.names = TRUE, facet.wrap.nrow = NULL, 
                         facet.wrap.ncol = NULL) {
  checkmate::assertClass(bmr, "BenchmarkResult")
  measure = mlr:::checkBMRMeasure(measure, bmr)
  checkmate::assertChoice(style, c("box", "violin"))
  df = mlr:::orderBMRLrns(bmr, df, order.lrns)
  df = mlr:::orderBMRTasks(bmr, df, order.tsks)
  if (pretty.names) {
    learner.short.names = getBMRLearnerShortNames(bmr)
    checkDuplicatedLearnerNames(learner.short.names)
    if (!is.null(order.lrns)) {
      learner.ids = getBMRLearnerIds(bmr)
      names(learner.short.names) = learner.ids
      learner.short.names = learner.short.names[order.lrns]
    }
    levels(df$learner.id) = learner.short.names
  }
  p = ggplot(df, aes_string("learner.id", measure$id)) + theme_bw()
  p = p + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = -45, hjust = 0), legend.position =  "none", plot.margin = margin(0.5, 1.5, 0.5, 0.5, "cm"))
  p = p + facet_wrap(~task.id, nrow = facet.wrap.nrow, ncol = facet.wrap.ncol)
  if (measure$name == "Uno's Concordance index") 
    p = p + ylab("cindex")
  if (style == "box") 
    p = p + geom_boxplot(aes(fill = learner.id)) 
    

  else p = p + geom_violin() + stat_summary(fun.ymin = median, 
                                            fun.ymax = median, fun.y = median, geom = "crossbar")
  return(p)
}

# 3. aggregation functions -------------------------------------------------------

# 3.1 aggregation for section using multi omics -----------

# aggregation over naive learners and sturcture learners
get_median_cindex <- function(cancer) {
  bmr_df_lists[[cancer]] %>%
    group_by(learner.id) %>%
    summarise(median(cindex.uno))
}

# cindex
get_groupmeans_cindex <- function(cancer, as_df = TRUE) {
  cans <- cancer
  res_means <- vector(mode = "list", length = length(cans))
  names(res_means) <- cans
  for (canc in cans) {
    # browser()
    means <- 
      bmr_df_lists[[canc]] %>%
      group_by(learner.id) %>%
      summarise(mean(cindex.uno))
    
    naiv <- mean(means[[2]][2:6])
    struc <- mean(means[[2]][8:11])
    res_means[[canc]] <- c("naiv" = naiv, "structure" = struc)
  }
  if (as_df) {
    as.data.frame(res_means)
  } else {
    res_means
  } 
}

# ibrier
get_groupmeans_ibrier <- function(cancer, as_df = TRUE) {
  cans <- cancer
  res_means <- vector(mode = "list", length = length(cans))
  names(res_means) <- cans
  for (canc in cans) {
    # browser()
    means <- 
      bmr_df_lists[[canc]] %>%
      group_by(learner.id) %>%
      summarise(mean(ibrier))
    
    naiv <- mean(means[[2]][2:6])
    struc <- mean(means[[2]][8:11])
    res_means[[canc]] <- c("naiv" = naiv, "structure" = struc)
  }
  if (as_df) {
    as.data.frame(res_means)
  } else {
    res_means
  } 
}


# 3.2 aggregation for added value -------------------------

learner_means <- function(df, cancer, meas = "cindex.uno", arg = FALSE) {
  meas <- sym(meas)
  df_temp <- 
    df %>%
    filter(task.id == cancer) %>%
    group_by(learner.id, add = TRUE) %>%
    summarize_(mean = expr(mean(!!meas)))
  if (arg) df_temp <- df_temp %>% arrange(desc(mean))
  df_temp
}


learner_medians <- function(df = df_add_val, cancer, meas = "cindex") {
  if (meas == "cindex") {
    df %>%
      filter(task.id == cancer) %>%
      group_by(learner.id, add = TRUE) %>%
      summarize(median = median(cindex.uno)) %>%
      arrange(desc(median))
  } else {
    df %>%
      filter(task.id == cancer) %>%
      group_by(learner.id, add = TRUE) %>%
      summarize(median = median(ibrier)) %>%
      arrange(median)
  }
}

# 3.3 aggregation over data sets --------------------------
# (mean of cv means)

get_mean_df <- function(dsets, lrns, base_df, meas = "cindex.uno") {
  l_means <- lapply(dsets, learner_means, df = base_df, meas = meas)
  names(l_means) <- dsets
  l_means <- lapply(l_means, function(x) x[[2]])
  df_means <- as.data.frame(l_means)
  row.names(df_means) <- lrns
  df_means
}
