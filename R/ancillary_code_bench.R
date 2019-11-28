#### 
#### sources & code  ####
#### 

# 1. custom learners ----------------------------------------------------------

#  1.1 surv.clinicalreference -------------------------------------------------
makeRLearner.surv.clinic_reference = function() {
  makeRLearnerSurv(
    cl = "surv.clinic_reference",
    package = "survival",
    par.set = makeParamSet(
      makeDiscreteLearnerParam(
        id = "model",
        default = "coxph",
        values = c("coxph", "lasso", "ridge")
      ), 
      makeUntypedLearnerParam(id = "clinicals", default = NULL), 
      makeIntegerLearnerParam(
        id = "maxit",
        default = 100000,
        lower = 10^4,
        upper = 10^6
      ), 
      makeIntegerLearnerParam(id = "nfolds", default = 10L, lower = 3L)
    ),
    properties = c("numerics", "factors", "weights", "prob"),
    name = "Reference model for survival prediction",
    short.name = "surv.clinic_ref",
    note = "Computes a reference models for survival analysis. This is most 
    usefull in (multi-)omics contexts. The reference is a model with 
    clinical-variables only. 
    The clinical-only model can be fit as a standard coxph model 
    (survival::coxph) or with regularisation via lasso or ridge 
    regression (glmnet::cv.glmnet).
    It is required to provide a vector of indices in the 'clinicals' 
    argument, indicating the position of clinical variables in the 
    feature matrix.",
    callees = c("coxph", "surv")
  )
}

trainLearner.surv.clinic_reference = function(.learner, .task, .subset, .weights = NULL, ...) {
  if (is.null(getLearnerParVals(.learner)$model)) {
    mod <- "coxph"
  } else {
    mod <- getLearnerParVals(.learner)$model
  }
  
  clinicals <- getLearnerParVals(.learner)$clinicals + 2
  data <- getTaskData(.task, subset = .subset)[, c(1, 2, clinicals)]
  
  if (is.null(clinicals)) {
    stop("No index vector for clinical variables supplied") 
  } 
  
  if ((getTaskSize(.task) <= length(clinicals)) & mod == "coxph") {
    warning("There are less observations than clinical variables. Fitting coxph 
            model is not possible. Lasso model fitted instead." )
    mod = "lasso"
  }
  
  if (mod != "coxph") {
    alpha = ifelse(mod == "lasso", 1, 0)
    args = list(x = as.matrix(data[, -c(1:2)]), 
                y = Surv(data$time, data$status), 
                family = "cox",
                type.measure = "deviance", 
                alpha = alpha, 
                maxit = getLearnerParVals(.learner)$maxit,
                nfolds = getLearnerParVals(.learner)$nfolds)
    model <- do.call(glmnet::cv.glmnet, args)
  } else {
    model <- coxph(Surv(time, status) ~ ., data = data, x = TRUE)
  }
  if (.learner$predict.type == "response") {
    model
  } else {
    if (mod != "coxph") {
      # check for NA coeffs and set to zero 
      if (any(is.na(model$coefficients))) {
        model$coefficients[is.na(model$coefficients)] = 0
      }
      # get index of nonzero coeffs
      coeffs = as.vector(coef(model, s = model$lambda.min))
      index_nz = get_nz(coeffs)
      
      cph = compute_coxph_mod(coeffs, index_nz, data = data) 
      model = mlr:::attachTrainingTime(model, .task, data)
      model$cph = cph
      model$index_nz = index_nz
      model
    } else {
      model = mlr:::attachTrainingTime(model, .task, data)
      model
    }
  }
}

predictLearner.surv.clinic_reference = function(.learner, .model, .newdata, ...) {
  coxreg <- inherits(.model$learner.model, "coxph")
  newdata = .newdata[, getLearnerParVals(.learner)$clinicals]
  
  if (coxreg) {
    preds = predict(.model$learner.model,
                    newdata = newdata,
                    type = "lp",
                    ...)
  } else {
    preds = as.numeric(predict(
      .model$learner.model,
      newx = as.matrix(newdata),
      type = "link",
      ...
    ))
  }
  
  if (.learner$predict.type == "response") {
    preds 
  } else {
    train.times = sort(unique(c(0, .model$learner.model$times)))
    if (coxreg) {
      if (any(is.na(.model$learner.model$coefficients))) {
        .model$learner.model$coefficients[is.na(.model$learner.model$coefficients)] = 0
      }
      probs = pec::predictSurvProb(.model$learner.model, newdata = newdata, times = train.times)
      list(preds = preds,
           probs = probs,
           train.times = train.times)
    } else {
      probs = pec::predictSurvProb(.model$learner.model$cph,
                                   newdata = newdata[, .model$learner.model$index_nz, drop = FALSE],
                                   times = train.times)
      list(preds = preds,
           probs = probs,
           train.times = train.times)
    }
  }  
}

#  1.2 surv.cv.glmnet (alias for lasso) -------------------------------------------
makeRLearner.surv.cv.glmnet2 = function() {
  makeRLearnerSurv(
    cl = "surv.cv.glmnet2",
    package = "glmnet",
    par.set = makeParamSet(
      makeNumericLearnerParam(id = "alpha", default = 1, lower = 0, upper = 1),
      makeIntegerLearnerParam(id = "nfolds", default = 10L, lower = 3L),
      makeDiscreteLearnerParam(id = "type.measure", values = "deviance", default = "deviance"),
      makeLogicalLearnerParam(id = "grouped", default = TRUE),
      makeDiscreteLearnerParam(id = "s", values = c("lambda.1se", "lambda.min"), default = "lambda.1se", when = "predict"),
      makeIntegerLearnerParam(id = "nlambda", default = 100L, lower = 1L),
      makeNumericLearnerParam(id = "lambda.min.ratio", lower = 0, upper = 1),
      makeLogicalLearnerParam(id = "standardize", default = TRUE),
      makeLogicalLearnerParam(id = "intercept", default = TRUE),
      makeNumericLearnerParam(id = "thresh", default = 1e-07, lower = 0),
      makeIntegerLearnerParam(id = "dfmax", lower = 0L),
      makeIntegerLearnerParam(id = "pmax", lower = 0L),
      makeIntegerVectorLearnerParam(id = "exclude", lower = 1L), # this is like a subset?
      makeNumericVectorLearnerParam(id = "penalty.factor", lower = 0, upper = 1),
      makeNumericVectorLearnerParam(id = "lower.limits", upper = 0),
      makeNumericVectorLearnerParam(id = "upper.limits", lower = 0),
      makeIntegerLearnerParam(id = "maxit", default = 100000L, lower = 1L),
      makeNumericLearnerParam(id = "fdev", default = 1.0e-5, lower = 0, upper = 1),
      makeNumericLearnerParam(id = "devmax", default = 0.999, lower = 0, upper = 1),
      makeNumericLearnerParam(id = "eps", default = 1.0e-6, lower = 0, upper = 1),
      makeNumericLearnerParam(id = "big", default = 9.9e35),
      makeIntegerLearnerParam(id = "mnlam", default = 5, lower = 1),
      makeNumericLearnerParam(id = "pmin", default = 1.0e-9, lower = 0, upper = 1),
      makeNumericLearnerParam(id = "exmx", default = 250.0),
      makeNumericLearnerParam(id = "prec", default = 1e-10),
      makeIntegerLearnerParam(id = "mxit", default = 100, lower = 1),
      makeLogicalLearnerParam(id = "factory", default = FALSE)
    ),
    properties = c("numerics", "factors", "ordered", "weights", "prob"),
    name = "GLM with Regularization (Cross Validated Lambda)",
    short.name = "cvglmnet",
    note = "Factors automatically get converted to dummy columns, ordered factors to integer.",
    callees = c("cv.glmnet", "glmnet", "glmnet.control", "predict.glmnet")
  )
}

trainLearner.surv.cv.glmnet2 = function(.learner, .task, .subset, .weights = NULL,  ...) {
  
  d = getTaskData(.task, 
                  subset = .subset, 
                  target.extra = TRUE, 
                  recode.target = "surv")
  
  info = mlr:::getFixDataInfo(d$data, 
                              factors.to.dummies = TRUE, 
                              ordered.to.int = TRUE)
  
  args = c(list(x = as.matrix(mlr:::fixDataForLearner(d$data, info)), 
                y = d$target, 
                family = "cox", 
                parallel = FALSE), 
           list(...))
  
  
  if (!is.null(.weights))
    args$weights = .weights
  
  saved.ctrl = glmnet::glmnet.control()
  is.ctrl.arg = names(args) %in% names(saved.ctrl)
  
  if (any(is.ctrl.arg)) {
    on.exit(do.call(glmnet::glmnet.control, saved.ctrl))
    do.call(glmnet::glmnet.control, args[is.ctrl.arg])
    args = args[!is.ctrl.arg]
  }
  
  model = mlr:::attachTrainingInfo(do.call(glmnet::cv.glmnet, args), info)  
  
  if (.learner$predict.type == "response") {  
    model
  } else {
    dat = getTaskData(.task, subset = .subset)
    
    # computing baseline survival
    coeffs = as.vector(coef(model, s = model$lambda.min))
    index_nz = get_nz(coeffs)
    
    cph = compute_coxph_mod(coeffs, index_nz, data = dat) 
    model = mlr:::attachTrainingTime(model, .task, dat)
    model$cph = cph
    model$index_nz <- index_nz
    model
  }
}

predictLearner.surv.cv.glmnet2 = function(.learner, .model, .newdata, ...) {
  info = mlr:::getTrainingInfo(.model)
  .newdata = mlr:::fixDataForLearner(.newdata, info)
  
  preds = as.numeric(predict(.model$learner.model, newx = as.matrix(.newdata), type = "link", ...))
  
  if (.learner$predict.type == "response") {
    preds
  } else {
    train.times = sort(unique(c(0, .model$learner.model$times)))
    probs = pec::predictSurvProb(.model$learner.model$cph,
                                 newdata = .newdata[, .model$learner.model$index_nz, drop = FALSE],
                                 times = train.times)
    list(preds = preds, probs = probs, train.times = train.times)
  }
}

#  1.3 surv.ts.ipflasso -------------------------------------------------------
makeRLearner.surv.ts.ipflasso = function() {
  makeRLearnerSurv(
    cl = "surv.ts.ipflasso",
    package = "ipflasso",
    par.set = makeParamSet(
      makeLogicalLearnerParam(id = "standardize", default = TRUE),
      makeNumericLearnerParam(id = "alpha", default = 1, lower = 0, upper = 1),
      makeUntypedLearnerParam(id = "blocks"),
      # makeNumericVectorLearnerParam(id = "pf", tunable = FALSE),
      makeIntegerLearnerParam(id = "nfolds", default = 5L, lower = 3L),
      makeIntegerLearnerParam(id = "ncv", default = 10L)
    ),
    properties = c("numerics", "factors", "prob"),
    name = "Two Step IPF-Lasso with prespecified penalty factors",
    short.name = "tsipflasso",
    note = "The penalty factors get prespecified via ridge regression on the training data. The block 
    structure has to be supplied."
  )
}

trainLearner.surv.ts.ipflasso = function(.learner, .task, .subset, ...) {
  
  # data, info and arguments
  d = getTaskData(
    .task,
    subset = .subset,
    target.extra = TRUE,
    recode.target = "surv"
  )
  
  info = mlr:::getFixDataInfo(d$data,
                              factors.to.dummies = TRUE,
                              ordered.to.int = TRUE)
  
  # blocks <- getLearnerParVals(.learner)$blocks
  
  # new args for cvr.adaptive.lasso
  args <- c(
    list(
      X = as.matrix(d$data),
      Y = d$target,
      family = "cox", # target variable datatype, i.e. survival
      type.measure = "deviance", # CV-criterion to optimize in cvr.ipflasso
      standardize = TRUE, # 
      alpha = 0, # Ridge regression in the first step
      type.step1 = "sep" # separate models in the first step
    ),
    list(...)  # blocks, ncv and nfolds are provided within the bench experiment routine (ncv = 1, nfolds = 10)
  )
  rm(d)
  
  # model fitting
  mod <- do.call(ipflasso::cvr.adaptive.ipflasso, args) # 1. and 2. step
  
  
  if (.learner$predict.type == "response") {
    if (sum(mod$coeff[ , mod$ind.bestlambda] != 0) == 0) { 
      class(mod) <- "FailureModel"
    } else {
      class(mod) <- "cvr.ipflasso" 
    }
    mod 
  } else {
    coeffs <- mod$coeff[-1, mod$ind.bestlambda] # -1 to remove intercept (which is NA for Cox!)
    index_nz = get_nz(coeffs) # if all coefs are zero, failure model results
    
    dat <- getTaskData(.task, subset = .subset)
    cph <- compute_coxph_mod(coeffs, index_nz, dat)
    mod$cph <- cph
    mod$index_nz <- index_nz
    mod = mlr:::attachTrainingTime(mod, .task, dat) 
    class(mod) <- "cvr.ipflasso"
    mod 
  }
}

predictLearner.surv.ts.ipflasso = function(.learner, .model, .newdata, ...) {
  
  mod = .model$learner.model
  preds = as.matrix(.newdata) %*% mod$coeff[-1, mod$ind.bestlambda]
  
  if (.learner$predict.type == "response") {
    preds
  } else {
    train.times = sort(unique(c(0, .model$learner.model$times)))
    probs = pec::predictSurvProb(.model$learner.model$cph,
                                 newdata = .newdata[, .model$learner.model$index_nz, drop = FALSE],
                                 times = train.times)
    list(preds = preds, probs = probs, train.times = train.times)
  }
}

#  1.4 surv.ts.priorlasso -----------------------------------------------------
makeRLearner.surv.ts.priorlasso = function() {
  makeRLearnerSurv(
    cl = "surv.ts.priorlasso",
    package = "prioritylasso",
    par.set = makeParamSet(
      makeUntypedLearnerParam(id = "blocks"),
      makeLogicalLearnerParam(id = "block1.penalization", default = TRUE),
      makeDiscreteLearnerParam(
        id = "lambda.type",
        default = "lambda.min",
        values = c("lambda.min", "lambda.1se"),
        tunable = FALSE
      ),
      makeLogicalLearnerParam(id = "standardize", default = TRUE),
      makeIntegerLearnerParam(
        id = "nfolds",
        default = 10L,
        lower = 2L
      ),
      makeLogicalLearnerParam(id = "cvoffset", default = FALSE),
      makeIntegerLearnerParam(
        id = "cvoffsetnfolds",
        default = 10L,
        lower = 2L,
        requires = quote(cvoffset == TRUE)
      ),
      makeLogicalLearnerParam(id = "favoring", default = FALSE)
    ),
    properties = c("numerics", "factors", "prob"),
    name = "Two step Priority Lasso",
    short.name = "ts.priorlasso",
    note = "Prioritylasso variante determining the priority datadriven in an separte step. Similiar to the 1. step
    in tsipf-lasso in the 1. step penalty factors for each modality are fit. According to these penalty factors the
    priortiy order is determined. The smaller (near zero) the penalty factor is, the greater is the importance of the
    associated modality. Ordering the penalty factors implies the priority order of the modalities."
  )
}

trainLearner.surv.ts.priorlasso = function(.learner, .task, .subset, ...) {
  
  d = getTaskData(
    .task,
    subset = .subset,
    target.extra = TRUE,
    recode.target = "surv"
  )
  
  info = mlr:::getFixDataInfo(d$data,
                              factors.to.dummies = TRUE,
                              ordered.to.int = TRUE)
  
  blocks <- getLearnerParVals(.learner)$blocks
  favoring <- getLearnerParVals(.learner)$favoring
  
  if (is.null(names(blocks))) {
    warning("The list of bocks is not named. Naming the blocks by 'block_i' 
            with i in {1, ..., number of blocks}")
    names(blocks) <- paste("block", seq(1:length(blocks)), sep = "_")
  }
  
  l_pf <- prespecify_pf(
    X = as.matrix(d$data),
    Y = d$target,
    blocks = blocks,
    separate = TRUE,
    alpha = 0
  )
  
  if(length(blocks) != length(l_pf$pf)) {
    stop("Length of penalty factors unequal to number of blocks.")
  }
  print("Penalty factors successfully computed. Continuing with model fitting")
  
  args <- list(
    X = as.matrix(d$data),
    Y = d$target,
    family = "cox",
    type.measure = "deviance"
  )
  
  # should clinical variables be favored
  if (favoring) {
    priority_order <- blocks[c(names(l_pf$pf[1]), names(sort(l_pf$pf[-1])))]
    args$blocks <- priority_order
    args$block1.penalization <- FALSE
  } else {
    priority_order <- blocks[names(sort(l_pf$pf))]
    args$blocks <- priority_order
  }
  
  print(c("Priortiy order: ", names(sort(l_pf$pf))), quote = FALSE)
  
  model <- do.call(prioritylasso::prioritylasso, args)
  model$priority_order <-  names(priority_order)
  
  if (.learner$predict.type == "response") {
    model
  } else {
    coeffs <- model$coefficients
    index_nz = get_nz(coeffs)
    
    dat <- getTaskData(.task, subset = .subset)
    cph <- compute_coxph_mod(coeffs, index_nz, dat)
    model$cph = cph
    model$index_nz <- index_nz
    model = mlr:::attachTrainingTime(model, .task, dat) 
    model
  }
}

predictLearner.surv.ts.priorlasso = function(.learner, .model, .newdata, ...) {
  preds <- predict(.model$learner.model, newdata = .newdata, type = "link", ...)
  
  if (.learner$predict.type == "response") {
    preds
  } else {
    train.times = sort(unique(c(0, .model$learner.model$times)))
    probs = pec::predictSurvProb(.model$learner.model$cph,
                                 newdata = .newdata[, .model$learner.model$index_nz, drop = FALSE],
                                 times = train.times)
    list(preds = preds, probs = probs, train.times = train.times)
  }
}

#  1.5 surv.cv.glmboost -------------------------------------------------------
makeRLearner.surv.cv.glmboost = function() {
  makeRLearnerSurv(
    cl = "surv.cv.glmboost",
    package = c("!survival", "mboost"),
    par.set = makeParamSet(
      makeDiscreteLearnerParam(
        id = "family",
        default = "CoxPH",
        values = c(
          "CoxPH",
          "Weibull",
          "Loglog",
          "Lognormal",
          "Gehan",
          "custom.family"
        )
      ),
      makeNumericVectorLearnerParam(
        id = "nuirange",
        default = c(0, 100),
        requires = quote(family %in% c("Weibull", "Loglog", "Lognormal"))
      ),
      makeUntypedLearnerParam(
        id = "custom.family.definition",
        requires = quote(family == "custom.family")),
      makeIntegerLearnerParam(
        id = "mstop",
        default = 100L,
        lower = 1L
      ),
      makeNumericLearnerParam(
        id = "nu",
        default = 0.1,
        lower = 0,
        upper = 1
      ),
      makeDiscreteLearnerParam(
        id = "risk",
        values = c("inbag", "oobag", "none")
      ),
      makeLogicalLearnerParam(
        id = "stopintern",
        default = FALSE
      ),
      # 'risk' and 'stopintern' will be kept for completeness sake
      makeLogicalLearnerParam(
        id = "center",
        default = FALSE
      ),
      makeLogicalLearnerParam(
        id = "use.formula",
        default = TRUE,
        when = "both"
      ),
      makeLogicalLearnerParam(
        id = "trace",
        default = FALSE,
        tunable = FALSE
      )
    ), 
    par.vals = list(
      family = "CoxPH",
      use.formula = TRUE
    ),
    properties = c("numerics", "factors", "ordered", "weights", "prob"),
    name = "Gradient Boosting with Componentwise Linear Models and automatically 
    tuned mstop.",
    short.name = "cv.glmboost",
    note = "`family` has been set to `CoxPH()` by default.",
    callees = c("glmboost", "mboost_fit", "boost_control", "CoxPH", 
                "Weibull", "Loglog", "Lognormal", "Gehan", "cvrisk")
  )
}

trainLearner.surv.cv.glmboost = function(.learner,
                                         .task,
                                         .subset,
                                         .weights = NULL,
                                         nuirange = c(0, 100),
                                         family,
                                         custom.family.definition,
                                         mstop,
                                         nu,
                                         risk,
                                         stopintern,
                                         trace,
                                         use.formula,
                                         ...) {
  ctrl = learnerArgsToControl(
    mboost::boost_control, 
    mstop, 
    nu, 
    risk, 
    trace, 
    stopintern
  )
  family = switch(
    family,
    CoxPH = mboost::CoxPH(),
    Weibull = mboost::Weibull(nuirange = nuirange),
    Loglog = mboost::Loglog(nuirange = nuirange),
    Lognormal = mboost::Lognormal(nuirange = nuirange),
    Gehan = mboost::Gehan(),
    custom.family = custom.family.definition
  )
  
  if (use.formula) {
    f = getTaskFormula(.task)
    model = if (is.null(.weights)) {
      mboost::glmboost(
        f,
        data = getTaskData(.task, subset = .subset, recode.target = "surv"),
        control = ctrl,
        family = family,
        ...
      )
    } else  {
      mboost::glmboost(
        f,
        data = getTaskData(.task, subset = .subset, recode.target = "surv"), 
        control = ctrl,
        weights = .weights,
        family = family,
        ...
      )
    }
  } else {
    data = getTaskData(
      .task,
      subset = .subset,
      target.extra = TRUE,
      recode.target = "surv"
    )
    info = mlr:::getFixDataInfo(data$data,
                                factors.to.dummies = TRUE,
                                ordered.to.int = TRUE)
    data$data = as.matrix(mlr:::fixDataForLearner(data$data, info))
    data$data = cbind(intercept = 1, data$data)
    
    model = if (is.null(.weights)) {
      mboost::glmboost(
        x = data$data,
        y = data$target,
        control = ctrl,
        family = family,
        ...
      )
    } else {
      mboost::glmboost(
        x = data$data,
        y = data$target,
        control = ctrl,
        weights = .weights,
        family = family,
        ...
      )
    }
  }
  
  # compute optimal number of boosting iterations and adjust model accordingly
  cv10f <- cv(model.weights(model), type = "kfold")
  cvm = cvrisk(model, folds = cv10f, grid = 0L:1000L, papply = lapply) 
  mstop(model) <- mstop(cvm)
  model = mlr:::attachTrainingInfo(model, info)
  
  if (.learner$predict.type == "response") {
    model
  } else {
    coeffs = coef(model)
    features = names(coeffs)
    dat = cbind(getTaskTargets(.task), data$data)
    sub_data = dat[, c("time", "status", features)]
    cph = coxph(
      Surv(time, status) ~ .,
      data = sub_data,
      init = coeffs,
      iter.max = 0,
      x = TRUE
    )
    model$cph = cph
    model$features = features
    model = mlr:::attachTrainingTime(model, .task, dat) 
    model
  }
}

predictLearner.surv.cv.glmboost = function(.learner,
                                           .model,
                                           .newdata,
                                           use.formula,
                                           ...) {
  if (!use.formula) {
    info = mlr:::getTrainingInfo(.model)
    newdata = mlr:::fixDataForLearner(.newdata, info)
    newdata = cbind(intercept = 1, newdata)
  }
  preds = predict(
    .model$learner.model,
    newdata = as.matrix(newdata),
    type = "link"
  )
  
  if (.learner$predict.type == "response") {
    preds
  } else {
    train.times = sort(unique(c(0, .model$learner.model$times)))
    probs = pec::predictSurvProb(
      .model$learner.model$cph,
      newdata = newdata[, .model$learner.model$features, drop = FALSE],
      times = train.times
    )
    list(preds = preds,
         probs = probs,
         train.times = train.times)
  }
}
#  1.6 surv.cv.coxboost -------------------------------------------------------
makeRLearner.surv.cv.CoxBoost2 = function() {
  makeRLearnerSurv(
    cl = "surv.cv.CoxBoost2",
    package = "!CoxBoost",
    par.set = makeParamSet(
      makeIntegerLearnerParam(id = "maxstepno", 
                              default = 100L, 
                              lower = 0L),
      makeIntegerLearnerParam(id = "K",
                              default = 10L,
                              lower = 1L), 
      makeDiscreteLearnerParam(id = "type", 
                               default = "verweij", 
                               values = c("verweij", "naive")),
      makeLogicalLearnerParam(id = "parallel",
                              default = FALSE,
                              tunable = FALSE), 
      makeLogicalLearnerParam(id = "upload.x",
                              default = FALSE,
                              tunable = FALSE), 
      makeLogicalLearnerParam(id = "multicore",
                              default = FALSE,
                              tunable = FALSE), 
      makeIntegerVectorLearnerParam(id = "unpen.index"),
      makeLogicalLearnerParam(id = "standardize", default = TRUE), 
      makeNumericLearnerParam(id = "penalty", lower = 0),
      makeDiscreteLearnerParam(id = "criterion", 
                               default = "pscore", 
                               values = c("pscore", "score", "hpscore", "hscore")),
      makeNumericLearnerParam(id = "stepsize.factor",
                              default = 1,
                              lower = 0), 
      makeLogicalLearnerParam(id = "trace",
                              default = FALSE,
                              tunable = FALSE)
    ),
    properties = c("numerics", "factors", "weights", "prob"),
    name = "Cox Proportional Hazards Model with Componentwise Likelihood based 
    Boosting, tuned for the optimal number of boosting steps",
    short.name = "cv.CoxBoost",
    note = "Factors automatically get converted to dummy columns, 
    ordered factors to integer.",
    callees = c("cv.CoxBoost", "CoxBoost")
  )
}

trainLearner.surv.cv.CoxBoost2 = function(.learner,
                                          .task,
                                          .subset,
                                          .weights = NULL,
                                          penalty = NULL,
                                          unpen.index = NULL,
                                          ...) {
  # preparing data and params -------------------------------------------------
  data = getTaskData(.task,
                     subset = .subset,
                     target.extra = TRUE,
                     recode.target = "surv")
  
  info = mlr:::getFixDataInfo(data$data,
                              factors.to.dummies = TRUE,
                              ordered.to.int = TRUE)
  
  unpen.index = getLearnerParVals(.learner)$`unpen.index`
  if (is.null(penalty))
    penalty = 9 * sum(data$target[, 2L])
  
  pars = c(list(time = data$target[, 1L],
                status = data$target[, 2L],
                x = as.matrix(mlr:::fixDataForLearner(data$data, info)),
                penalty = penalty,
                weights = .weights,
                unpen.index = unpen.index), 
           list(...))
  rm(data)
  
  # tuning --------------------------------------------------------------------
  # determining the optimal number of boosting iterations
  res = do.call(CoxBoost::cv.CoxBoost, pars)
  res$optimal.step
  if (res$optimal.step == 0L)
    warning("Could not determine the optimal step number in cv.CoxBoost")
  
  pars = BBmisc:::insert(pars, list(stepno = res$optimal.step))
  pars$maxstepno = NULL
  
  # model fitting -------------------------------------------------------------
  # with optimal number of boosting iterations
  model = mlr:::attachTrainingInfo(do.call(CoxBoost::CoxBoost, pars), info)
  
  if (.learner$predict.type == "response") {
    model
  } else {
    dat = getTaskData(.task, subset = .subset)
    model = mlr:::attachTrainingTime(model, .task, dat)
    model
  }
}

predictLearner.surv.cv.CoxBoost2 = function(.learner, .model, .newdata, ...) {
  info = mlr:::getTrainingInfo(.model)
  .newdata = as.matrix(mlr:::fixDataForLearner(.newdata, info))
  preds = as.numeric(predict(.model$learner.model, 
                             newdata = .newdata, 
                             type = "lp"))
  
  if (.learner$predict.type == "response") {
    preds
  } else {
    train.times = sort(unique(c(0, .model$learner.model$times)))
    probs = predict(.model$learner.model, 
                    .newdata, type = "risk", 
                    times = train.times)
    
    list(preds = preds,
         probs = probs,
         train.times = train.times)
  }
}

#  1.7 surv.blockForest -------------------------------------------------------
makeRLearner.surv.blockForest = function() {
  makeRLearnerSurv(
    cl = "surv.blockForest",
    package = "blockForest",
    par.set = makeParamSet(
      makeUntypedLearnerParam(id = "blocks"),
      makeIntegerLearnerParam(id = "num.trees", default = 2000L, lower = 1L),
      makeIntegerLearnerParam(id = "mtry", default = NULL, special.vals = list(NULL)),
      makeIntegerLearnerParam(id = "nsets", default = 300, lower = 1),
      makeIntegerLearnerParam(id = "num.trees.pre", default = 1500, lower = 1),
      makeIntegerLearnerParam(id = "num.threads", default = NULL, special.vals = list(NULL)),
      makeDiscreteLearnerParam(id = "splitrule", default = "extratrees", values = c("extratrees", "logrank", "C", "maxstat"))
    ),
    properties = c("missings", "numerics", "factors", "weights", "prob", "rcens"),
    name = "Random Forest variants for block-structured covariate data",
    short.name = "blockForest",
    note = ""
  )
}

trainLearner.surv.blockForest = function (.learner, .task, .subset, .weights = NULL, ...) {
  d = getTaskData(.task, subset = .subset, target.extra = TRUE)
  info = mlr:::getFixDataInfo(d$data,
                              factors.to.dummies = TRUE,
                              ordered.to.int = TRUE)
  
  args = c(list(X = as.matrix(d$data), 
                y = as.matrix(d$target),
                block.method = "BlockForest"), 
           list(...))
  
  model = do.call(blockForest::blockfor, args)
  dat = getTaskData(.task, subset = .subset)
  model = mlr:::attachTrainingTime(model, .task, dat)
  model
}

predictLearner.surv.blockForest = function (.learner, .model, .newdata, ...) {
  p = predict(object = .model$learner.model$forest, data = .newdata)
  preds = rowMeans(p$chf)
  
  if (.learner$predict.type == "response") {
    preds
  } else {
    train.times = sort(unique(c(0, .model$learner.model$times)))
    ptemp = p$survival
    pos = prodlim::sindex(jump.times = p$unique.death.times, eval.times = train.times)
    probs = cbind(1, ptemp)[, pos + 1, drop = FALSE]
    list(preds = preds, probs = probs, train.times = train.times)
  }
}

#  1.8 surv.cvSGL ---------------------------------------------------------------
makeRLearner.surv.cvSGL = function() {
  makeRLearnerSurv(
    cl = "surv.cvSGL",
    package = "SGL",
    par.set = makeParamSet(
      makeUntypedLearnerParam(id = "index", default = NULL),
      makeIntegerLearnerParam(id = "maxit", default = 1000),
      makeNumericLearnerParam(id = "thresh", default = 0.001),
      makeNumericLearnerParam(id = "min.frac", default = 0.1),
      makeIntegerLearnerParam(id = "nlam", default = 20),
      makeNumericLearnerParam(id = "gamma", default = 0.8, lower = 0, upper = 1),
      makeIntegerLearnerParam(id = "nfold", default = 10L, lower = 3),
      makeLogicalLearnerParam(id = "standardize", default = TRUE),
      makeLogicalLearnerParam(id = "verbose", default = FALSE),
      makeNumericLearnerParam(id = "step", default = 1, lower = 0, upper = 1),
      makeIntegerLearnerParam(id = "reset", default = 10),
      makeNumericLearnerParam(id = "alpha", default = 0.95, lower = 0, upper = 1),
      makeUntypedLearnerParam(id = "lambda", default = NULL)
    ),
    properties = c("numerics", "factors", "prob"),
    name = "Cross-validated Sparse Group Lasso",
    short.name = "cvSGL",
    note = ""
  )
}

trainLearner.surv.cvSGL = function(.learner, .task, .subset, ...) {
  d = getTaskData(.task, subset = .subset, target.extra = TRUE)
  info = mlr:::getFixDataInfo(d$data, factors.to.dummies = TRUE, ordered.to.int = TRUE)
  
  # adjustment for all-zero-vars ---------------------------------------------------
  if (any(apply(d$data, 2, function(x) sum(x) == 0))) {
    d$data = d$data[, -which(apply(d$data, 2, function(x) sum(x == 0)) == nrow(d$data))]
  }
  
  # block adjustment ---------------------------------------------------------------
  blks = get_blocks2(d$data) 
  
  bl_sizes = vapply(blks, length, FUN.VALUE = integer(1L)) 
  sgl_blocks = unlist(mapply(rep.int, times =  bl_sizes, x = seq_along(blks)))
  ind = which(names(list(...)) == "index")
  
  # prepare args -------------------------------------------------------------------
  args = c(list(data = list(x = as.matrix(d$data), 
                            time = d$target$time, 
                            status = d$target$status),
                type = "cox",
                index = sgl_blocks), 
           list(...)[-ind])
  
  # SGL model --------------------------------------------------------------------------
  mod = do.call(SGL::cvSGL, args)
  
  # coefficient extraction and check----------------------------------------------------
  index_nz <- FALSE
  count <- which.min(mod$lldiff)
  while (!any(index_nz)) {
    if (count > length(mod$lldiff)) {
      stop("Null model fitted.")
    }
    
    coeffs = mod$fit$beta[, count]
    index_nz = get_nz(coeffs, check_all_zero = FALSE)
    
    count = count + 1
  }
  
  # coxph model for predicted probs ------------------------------------------------
  dat =  cbind(d$target, d$data)
  cph = compute_coxph_mod(coeffs, index_nz, data = dat) 
  model = mlr:::attachTrainingTime(mod, .task, dat)
  model$cph = cph
  model$index_nz = index_nz
  model
}

predictLearner.surv.cvSGL = function(.learner, .model, .newdata, ...) {
  # Adjust new data matrix for deleted all-zero vars.
  i_selected = names(.newdata) %in% attributes(lmod$fit$X)$dimnames[[2]]
  newX = .newdata[, i_selected]
  
  # compute lin preds 
  coeffs = lmod$cph$coefficients
  vars = attributes(coeffs)$names
  tempX = as.matrix(newX[vars])
  eta = tempX %*% coeffs
  preds = exp(eta) 
  
  if (.learner$predict.type == "response") {
    preds
  } else {
    train.times = sort(unique(c(0, .model$learner.model$times)))
    probs = pec::predictSurvProb(.model$learner.model$cph,
                                 newdata = newX[, .model$learner.model$index_nz, drop = FALSE],
                                 times = train.times)
    list(preds = preds,
         probs = probs,
         train.times = train.times)
  }
}

#  1.9 surv.grridge ----------------------------------------------------------
makeRLearner.surv.grridge = function() {
  makeRLearnerSurv(
    cl = "surv.grridge",
    package = "GRridge",
    par.set = makeParamSet(
      makeUntypedLearnerParam(id = "partitions"),                               # normal blocks (or more general list of list)
      makeUntypedLearnerParam(id = "unpenal"),                                  # needs to be formula
      makeDiscreteLearnerParam(id = "method", 
                               default = "exactstable", 
                               values = c("exactstable", "stable", "exact")),
      makeIntegerLearnerParam(id = "niter", 
                              default = 10, 
                              lower = 1),
      makeLogicalLearnerParam(id = "monotone", 
                              default = NULL, 
                              special.vals = list(NULL)),
      makeNumericLearnerParam(id = "oplt", 
                              lower = 0, 
                              upper = Inf, 
                              default = NULL, 
                              special.vals = list(NULL)),
      makeIntegerLearnerParam(id = "innfold", 
                              lower = 0, 
                              default = NULL, 
                              special.vals = list(NULL)),
      makeLogicalLearnerParam(id = "selectionEN", 
                              default = FALSE, 
                              special.vals = list(NULL)),                       # needs to be true for bench experiment
      makeIntegerVectorLearnerParam(id = "maxsel", 
                                    default = c(25,100),
                                    special.vals = list(NULL)),                 # if set to NULL the max. no of selc. features equals nrow
      makeUntypedLearnerParam(id = "dataunpen"),
      makeLogicalLearnerParam(id = "standardizeX", 
                              default = TRUE, 
                              special.vals = list(NULL)),
      makeDiscreteLearnerParam(id = "savepredobj", 
                               default = "all", 
                               values = c("all", "last", "none"))
    ),
    properties = c("missings", "numerics", "factors", "weights", "prob", "rcens"),
    name = "Group-regularized (logistic) ridge regression",
    short.name = "GRridge",
    note = ""
  )
}

trainLearner.surv.grridge = function (.learner, .task, .subset, .weights = NULL, ...) 
{
  d = getTaskData(.task, subset = .subset, target.extra = TRUE, recode.target = "surv")
  info = mlr:::getFixDataInfo(d$data, factors.to.dummies = TRUE, ordered.to.int = TRUE)
  
  # adjustment for all-zero-vars ---------------------------------------------------
  if (any(apply(d$data, 2, function(x) sum(x) == 0))) {
    all_zeros <- which(apply(d$data, 2, function(x) sum(x == 0)) == nrow(d$data))
    d$data = d$data[, -all_zeros]
  }
  
  args = c(list(highdimdata = t(d$data), 
                response = d$target),
           list(...))
  
  if (is.null(args$maxsel))
    args$maxsel <- ncol(d$data)
  
  # 1. First CV the lambdas separately for the groups of different variabels.
  # 2. The square-root of the inverses of these 5 penalties would then be the 
  #    initial group weights, wt in the code below
  # 3. Then apply these as initial weight for GRridge
  # 4. GRidge will then moderate the penalty weights, thereby appropriately 
  #    accounting for possible collinearities across the 5 groups of variables
  blocks <- get_blocks2(d$data)
  lambdas <- vector("integer", length = length(blocks))
  args_optL2 <- list(response = args$response,
                     fold = 5L)
  scaled_dat <- as.data.frame(scale(d$data))
  
  for (opt in seq_along(blocks)) {
    args_optL2$penalized <- scaled_dat[, blocks[[opt]]]
    lambdas[opt] <- do.call(penalized::optL2, args = args_optL2)$lambda
  }
  ngr <- vapply(blocks, length, FUN.VALUE = integer(1))
  wt <- function(i) rep(1/sqrt(lambdas[i]), ngr[i])
  weights0 <- unlist(sapply(seq_along(blocks), wt))
  weights <- weights0/mean(weights0)
  args$highdimdata <- t(d$data * weights)
  args$partitions <- blocks
  args$method <- "exactstable"
  
  model <- do.call(GRridge::grridge, args)
  
  # model coefficients                                                          
  coeffs <- model$betas
  
  # setting NA-coeffs to zero for stability. Will not be used in model
  if (any(is.na(coeffs))) {
    coeffs[is.na(coeffs)] = 0
  }
  
  # routine only works for one maxsel safely! Adjust!
  index_nz <- names(d$data) %in% names(model$resEN[[1]]$betasEN)                # creating logicals instead of using indices since compute_coxph_mod needs logicals
  
  # check if all coeffs are zero and stop make failure model
  if (!any(index_nz)) {
    stop("Null model fitted.")
  }
  
  dat =  cbind(getTaskTargets(.task), d$data)
  cph <- compute_coxph_mod(coeffs, index_nz, dat)                               # kills R if to many non-zero coefficients
  model$cph = cph
  model$index_nz <- index_nz
  model = mlr:::attachTrainingTime(model, .task, dat) # check for problems
  model$weights <- weights
  if (exists("all_zeros")) {
    model$all_zeros <- all_zeros
  } else {
    model$all_zeros <- NULL
  }  
  class(model) <- "grridge"
  model
}

predictLearner.surv.grridge = function (.learner, .model, .newdata, ...) {
  # Adjust new data matrix for deleted all-zero vars and scale by test data
  # newdata differs in ncol to design matrix X since X had been adjusted for all zeros
  newX <- .newdata[, -lmod$all_zeros]
  
  # weighting since special routine by mvd used
  newX <- newX*lmod$weights
  
  i_selected <- names(.newdata) %in% names(.model$learner.model$resEN[[1]]$betasEN) 
  newX <- .newdata[, i_selected, drop = FALSE]
  
  coeffs <- coef(lmod)
  vars <- attributes(coeffs)$names
  tempX <- as.matrix(newX[vars])
  eta <- tempX %*% coeffs
  preds <- exp(eta)
  
  if (.learner$predict.type == "response") {
    preds
  } else {
    train.times = sort(unique(c(0, .model$learner.model$times)))
    probs = pec::predictSurvProb(.model$learner.model$cph, 
                                 newdata = as.data.frame(newX),
                                 times = train.times)
    
    list(preds = preds, probs = probs, train.times = train.times)
  }
}

#  1.10 surv.tunedrfsrc ---------------------------------------------------------
makeRLearner.surv.randomForestSRC2 = function() {
  makeRLearnerSurv(
    cl = "surv.randomForestSRC2",
    package = c("survival", "randomForestSRC"),
    par.set = makeParamSet(
      makeIntegerLearnerParam(id = "ntree", default = 1000L, lower = 1L),
      makeDiscreteLearnerParam(id = "bootstrap", default = "by.root",
                               values = c("by.root", "by.node", "none")),
      makeIntegerLearnerParam(id = "mtry", lower = 1L),
      makeIntegerLearnerParam(id = "nodesize", lower = 1L, default = 3L),
      makeIntegerLearnerParam(id = "nodedepth", default = -1L),
      makeDiscreteLearnerParam(id = "splitrule", default = "logrank",
                               values = c("logrank", "logrankscore", "random")),
      makeIntegerLearnerParam(id = "nsplit", lower = 0L, default = 0L,
                              requires = quote(splitrule != "random")), # nsplit is ignored and internally set to 1 for splitrule = "random"
      makeLogicalLearnerParam(id = "split.null", default = FALSE),
      makeDiscreteLearnerParam(id = "importance", default = FALSE, tunable = FALSE,
                               values = list(`FALSE` = FALSE, `TRUE` = TRUE, "none", "permute", "random",
                                             "anti", "permute.ensemble", "random.ensemble", "anti.ensemble")),
      makeDiscreteLearnerParam(id = "na.action", default = "na.impute",
                               values = c("na.omit", "na.impute"), when = "both"),
      makeIntegerLearnerParam(id = "nimpute", default = 1L, lower = 1L),
      makeUntypedLearnerParam(id = "ntime"), # can be a single integer with number of time points or a numeric vector of time values
      makeDiscreteLearnerParam(id = "proximity", default = FALSE, tunable = FALSE,
                               values = list("inbag", "oob", "all", `TRUE` = TRUE, `FALSE` = FALSE)),
      makeIntegerLearnerParam(id = "sampsize", lower = 1L,
                              requires = quote(bootstrap == "by.root")),
      makeDiscreteLearnerParam(id = "samptype", default = "swr", values = c("swr", "swor"),
                               requires = quote(bootstrap == "by.root")),
      makeNumericVectorLearnerParam(id = "xvar.wt", lower = 0),
      makeLogicalLearnerParam(id = "forest", default = TRUE, tunable = FALSE),
      makeDiscreteLearnerParam(id = "var.used", default = FALSE, tunable = FALSE,
                               values = list(`FALSE` = FALSE, "all.trees", "by.tree")),
      makeDiscreteLearnerParam(id = "split.depth", default = FALSE, tunable = FALSE,
                               values = list(`FALSE` = FALSE, "all.trees", "by.tree")),
      makeIntegerLearnerParam(id = "seed", upper = 0L, tunable = FALSE),
      makeLogicalLearnerParam(id = "do.trace", default = FALSE, tunable = FALSE, when = "both"), # is currently ignored
      makeLogicalLearnerParam(id = "membership", default = TRUE, tunable = FALSE),
      makeLogicalLearnerParam(id = "statistics", default = FALSE, tunable = FALSE),
      makeLogicalLearnerParam(id = "tree.err", default = FALSE, tunable = FALSE)
    ),
    par.vals = list(na.action = "na.impute"),
    properties = c("missings", "numerics", "factors", "ordered", "weights", "oobpreds", "featimp", "prob"),
    name = "Random Forest",
    short.name = "rfsrc",
    note = '`na.action` has been set to `"na.impute"` by default to allow missing data support.',
    callees = "rfsrc"
  )
}

trainLearner.surv.randomForestSRC2 = function(.learner, .task, .subset, .weights = NULL, ...) {
  data = getTaskData(.task, subset = .subset)
  f = getTaskFormula(.task)
  model = randomForestSRC::tune(f, data = data, case.wt = .weights, doBest = TRUE, ...)
  
  if (.learner$predict.type == "response") {
    model
  } else {
    model = model$rf
    model = mlr:::attachTrainingTime(model, .task, data)
    model
  }
}

predictLearner.surv.randomForestSRC2 = function(.learner, .model, .newdata, ...) {
  preds = predict(.model$learner.model, newdata = .newdata, membership = FALSE, ...)$predicted
  if (.learner$predict.type == "response") {
    preds
  } else {
    train.times = sort(unique(c(0, .model$learner.model$times)))
    probs = pec::predictSurvProb(.model$learner.model, newdata = .newdata, times = train.times)
    list(preds = preds, probs = probs, train.times = train.times)
  }
}

#  1.11 surv.tunedranger -------------------------------------------------------
makeRLearner.surv.tuneMtryFast2 = function() {
  makeRLearnerSurv(
    cl = "surv.tuneMtryFast2",
    package = "tuneRanger",
    par.set = makeParamSet(
      makeIntegerLearnerParam(id = "mtryStart", lower = 1L),
      makeIntegerLearnerParam(id = "num.treesTry", lower = 1L, default = 50L),
      makeNumericLearnerParam(id = "stepFactor", lower = 0L),
      makeNumericLearnerParam(id = "improve", lower = 0L),
      makeLogicalLearnerParam(id = "trace", default = TRUE),
      makeIntegerLearnerParam(id = "min.node.size", lower = 1L),
      makeLogicalLearnerParam(id = "replace", default = TRUE),
      makeNumericLearnerParam(id = "sample.fraction", lower = 0L, upper = 1L),
      makeNumericVectorLearnerParam(id = "split.select.weights", lower = 0, upper = 1),
      makeUntypedLearnerParam(id = "always.split.variables"),
      makeDiscreteLearnerParam("respect.unordered.factors", values = c("ignore", "order", "partition"), default = "ignore"),
      makeDiscreteLearnerParam(id = "importance", values = c("none", "impurity", "permutation"), default = "none", tunable = FALSE),
      makeLogicalLearnerParam(id = "write.forest", default = TRUE, tunable = FALSE),
      makeLogicalLearnerParam(id = "scale.permutation.importance", default = FALSE, requires = quote(importance == "permutation"), tunable = FALSE),
      makeIntegerLearnerParam(id = "num.threads", lower = 1L, when = "both", tunable = FALSE),
      makeLogicalLearnerParam(id = "save.memory", default = FALSE, tunable = FALSE),
      makeLogicalLearnerParam(id = "verbose", default = TRUE, when = "both", tunable = FALSE),
      makeIntegerLearnerParam(id = "seed", when = "both", tunable = FALSE),
      makeDiscreteLearnerParam(id = "splitrule", values = c("gini", "extratrees"), default = "gini"),
      makeIntegerLearnerParam(id = "num.random.splits", lower = 1L, default = 1L, requires = quote(splitrule == "extratrees")),
      makeLogicalLearnerParam(id = "keep.inbag", default = FALSE, tunable = FALSE)
    ),
    properties = c("numerics", "factors", "ordered",  "weights", "prob"),
    name = "tuneMtryFast for ranger",
    short.name = "tuneMtryFast",
    note = ""
  )
}

trainLearner.surv.tuneMtryFast2 = function(.learner, .task, .subset, .weights = NULL, classwt = NULL, cutoff, ...) {
  data = getTaskData(.task, .subset)
  tn = getTaskTargetNames(.task)
  model = tuneRanger::tuneMtryFast(
    formula = NULL,
    data = data,
    dependent.variable.name = tn[1L],
    status.variable.name = tn[2L],
    num.treesTry = 50,
    doBest = TRUE,
    case.weights = .weights,
    ...
  )
  model = mlr:::attachTrainingTime(model, .task, data)
  model
}

predictLearner.surv.tuneMtryFast2 = function(.learner, .model, .newdata, ...) {
  p = predict(object = .model$learner.model, data = .newdata)
  preds = rowMeans(p$chf)
  if (.learner$predict.type == "response") {
    preds
  } else {
    train.times = sort(unique(c(0, .model$learner.model$times)))
    ptemp = p$survival
    pos = prodlim::sindex(jump.times = p$unique.death.times, eval.times = train.times)
    probs = cbind(1, ptemp)[, pos + 1, drop = FALSE]
    list(preds = preds, probs = probs, train.times = train.times)
  }
}

#  Help functions --------------------------------------------------------------
compute_coxph_mod <- function(coeffs, index_nz, data, ...) {
  index_nz_y <- c(TRUE, TRUE, index_nz)
  sub_data <- data[, index_nz_y, drop = FALSE]
  coxph(
    Surv(time, status) ~ .,
    data = sub_data,
    init = coeffs[index_nz],
    iter.max = 0,
    x = TRUE
  )
}

check_coeffs <- function(coeffs) {
  # check for NA coeffs and set to zero 
  if (any(is.na(coeffs))) {
    coeffs[is.na(coeffs)] = 0
  }
  coeffs
}

make_nz <- function(coeffs, check_all_zero = TRUE) {
  index_nz <- coeffs != 0
  
  if (check_all_zero) {
    if (!any(index_nz)) {
      stop("Null model fitted.")
    }
  }
  
  index_nz
}

get_nz <- function(coeffs, check_all_zero = TRUE) {
  coeffs %>%
    check_coeffs() %>%
    make_nz()
}


#  End 1. ---------------------------------------------------------------------

# 2. custom measures ----------------------------------------------------------

get_selected_feats = function(task, model, pred, feats, extra.args) {
  if (inherits(model, "FailureModel")) {
    return(NA_real_)
  }
  
  if (inherits(model$learner.model, "cv.glmnet")) {
    coefs <- coef.cv.glmnet2(model$learner.model)
  } else {
    coefs <- coef(model$learner.model)  
  }
  
  if (any(is.na(coefs))) {
    coefs[is.na(coefs)] <- 0
  }
  
  if (any(coefs == 0)) {
    coefs <- coefs[coefs != 0]
  }
  
  # get selected features
  if (extra.args$group == "default") {
    length(coefs)
  } else {
    sel_gr_feats <- grep(pattern = extra.args$group,
                         x = names(coefs),
                         value = TRUE)
    
    length(sel_gr_feats)
  }
}

#  2.1 group features ----------------------------------------------------------
group_meas <-
  c("featselc_cnv",
    "featselc_rna",
    "featselc_mirna",
    "featselc_clin",
    "featselc_mutation")

properties = c(
  "classif",
  "classif.multi",
  "multilabel",
  "regr",
  "surv",
  "costsens",
  "cluster",
  "req.model",
  "req.pred"
)

for (m in group_meas) {
  grp <- strsplit(m, "_")[[1]][2]
  
  temp <- makeMeasure(
    id = m,
    minimize = TRUE,
    best = 0,
    worst = Inf,
    properties = properties,
    name = paste0("Amount of features of feature group '", grp, "' used for the model"),
    note = "Useful for feature selection. The feature group is specified via extra.arg.
    The names of the features within the group must contain that string. Make sure names of
    features of other groups do not contain this string.",
    fun =  get_selected_feats,
    extra.args = list(group = paste0("_", grp)),
    aggr = test.mean
  )
  
  assign(m, temp)
}

#  2.2 featselc_default -------------------------------------------------------
featselc_default <-
  makeMeasure(
    id = "featselc_default",
    minimize = TRUE,
    best = 0,
    worst = Inf,
    properties = properties,
    name = "Overall amount of selected features.",
    note = "Useful for feature selection.",
    fun = get_selected_feats,
    extra.args = list(group = "default"),
    aggr = test.mean
  )

#  Help functions --------------------------------------------------------------
coef.cv.SGL <- function(model) {
  model$cph$coefficients
}

coef.grridge <- function(model) {
  model$resEN[[1]]$betasEN
}

coef.cv.glmnet2 <- function(model) {
  lambda_min <- model$lambda.min
  lambda_ind <- which(model$lambda == lambda_min)
  index_nz <- model$glmnet.fit$beta[, lambda_ind] != 0
  coefs <- model$glmnet.fit$beta[index_nz, lambda_ind]
}

coef.cvr.ipflasso <- function(model) {
  model$coeff[, model$ind.bestlambda]
}

#  End 2.----------------------------------------------------------------------

# 3. task specific learners ---------------------------------------------------

#  3.1 Main function ----------------------------------------------------------
make_spec_lrns <- function(task, lrnr, args) {
  
  # derive blocks -------------------------------------------------------------
  data <- mlr:::getTaskData(task)
  blocks <- get_blocks2(data)
  clinicals <- blocks$clinical
  
  # create task specific learners ---------------------------------------------
  args <- args[[lrnr]]
  
  # Clinical reference
  if (lrnr == "lrn_clin_ref") {
    args$clinicals <- clinicals
    lrn_clin_ref <- do.call(makeLearner, args = args)
  }
  
  # TS priority lasso
  if (lrnr == "lrn_ts_prior") {
    args$blocks <- blocks
    lrn_ts_prior <- do.call(makeLearner, args = args)
  }
  
  # TS priority lasso favoring
  if (lrnr == "lrn_ts_prior_fav") {
    args$blocks <- blocks
    lrn_ts_prior_fav <- do.call(makeLearner, args = args)
  }
  
  # IPF-Lasso
  if (lrnr == "lrn_tsipf") {
    args$blocks <- blocks
    lrn_tsipf <- do.call(makeLearner, args = args)
  }
  
  # mandatory boostings
  if (lrnr == "lrn_cv_coxboost_unpen") {
    args$unpen.index <- clinicals
    lrn_cv_coxboost_unpen <- do.call(makeLearner, args = args)
  }
  
  # blockForest
  if (lrnr == "lrn_blockForest") {
    args$blocks <- blocks
    lrn_blockForest <- do.call(makeLearner, args = args)
  }
  
  # SGL
  if (lrnr == "lrn_SGL") {
    args$index <- blocks
    lrn_SGL <- do.call(makeLearner, args = args)
  }
  
  # ggridge 
  if (lrnr == "lrn_grridge") {
    args$partitions <- blocks
    lrn_grridge <- do.call(makeLearner, args = args)
  }
  
  list(get(lrnr))
}

#  3.2 Help functions ---------------------------------------------------------

#   3.2.1 Get groups sizes ----------------------------------------------------
get_grp_size = function(group, data, sep = "_") {
  gr_vars <- grep(paste0(sep, group), names(data))
  length(gr_vars)
}

#   3.2.2 Get indices w.r.t. block membership ---------------------------------
get_blocks2 <- function(data, groups = c("clinical", "cnv", "mirna", "mutation", "rna")) {
  blength <- lapply(groups, function(group, data) get_grp_size(group, data), data = data)
  nam <- groups
  
  if (any(blength == 0)) {
    ind <- which(blength == 0)
    blength <- blength[-ind]
    nam <- nam[-ind]
  }
  
  blocks <- rep(seq_along(blength), times = blength)
  blocks <- lapply(seq_along(blength), function(x) which(blocks == x))
  names(blocks) <- nam
  blocks
}

#   3.2.4 Prespecify penalty factors for 2-step ipflasso/priorlasso------------
#    3.2.4.1 Main function ----------------------------------------------------
prespecify_pf <- function(X, 
                          Y, # Needs to be Surv-Object
                          family = "cox", 
                          standardize = TRUE, 
                          nfolds = 5, 
                          type_measure = "deviance",
                          blocks,
                          separate = TRUE, 
                          alpha = 1, 
                          mean_type = "arith") {
  
  l_args <- list(x = as.matrix(X), 
                 y = Y, 
                 family = family, 
                 standardize = standardize, 
                 nfolds = nfolds, 
                 type.measure = type_measure, 
                 alpha = alpha)
  
  # compute mean coeffs 
  if (separate) {
    means <- compute_coeffs_sep(l_args = l_args, blocks = blocks, mean_type = mean_type)
  } else {
    means <- compute_coeffs_comb(l_args = l_args, blocks = blocks, mean_type = mean_type)
  }
  
  # get penalty factors and subset to non zero effect variables
  if (any(means != 0)) {
    means_check <- 1/means  
    
    exc <- create_exclusion_vector(blocks = blocks, means_check = means_check)  
    pf <- means_check[is.finite(means_check)]
    new_blocks <- subset_blocks(blocks = blocks, means_check)
    
    return(list(exc = exc, pf = pf, new_blocks = new_blocks))
  } else {
    warning("All step 1 coefficients are zero")
  }
}

#    3.2.4.2 Help functions ---------------------------------------------------
compute_coeffs_comb <- function(l_args, mean_type, blocks) {
  
  # model 
  mod <- do.call("cv.glmnet", l_args)
  # mod <- cv.glmnet(x = X, y = Y, family = family, standardize = standardize, nfolds = nfolds, type.measure = type.measure, alpha = alpha)
  
  # coeffs 
  coeffs <- mod %>%
    coef(s = "lambda.min") %>%
    abs() %>%
    as.matrix()
  
  # indexing blocks 
  block_length <- vapply(blocks, length, FUN.VALUE = integer(1))
  n_blocks <- length(blocks)
  index = integer(length = 0)
  for (i in 1L:n_blocks) {
    index <- c(index, rep(i,block_length[i]))
  }
  
  # means
  if (mean_type == "arith") {
    means <- as.numeric(tapply(coeffs[,1], index, mean))
  } else {
    means <- as.numeric(tapply(coeffs[,1], index, geometric.mean))
  }
  return(means)
}

compute_coeffs_sep <- function(l_args, blocks, mean_type) {
  
  means <- c(rep(0,length(blocks)))
  names(means) <- names(blocks) 
  
  # compute separate models for every modality
  for (l in 1:length(blocks)) {
    temp_args <- l_args
    temp_args[[1]] <- temp_args[[1]][, blocks[[l]]]
    
    # model
    assign(paste("model",l,sep = ""), do.call("cv.glmnet", temp_args))
    
    # coefficients
    abs_coeffs <- 
      paste("model", l , sep = "") %>%
      get() %>%
      coef(s = "lambda.min") %>%
      abs()
    
    # mean
    if (mean_type == "arith") {
      means[l] <- mean(abs_coeffs)
    } else {
      means[l] <- geometric.mean(abs_coeffs)
    }
  }
  return(means)
}

create_exclusion_vector <- function(blocks, means_check) {
  exc <- NULL
  a = 1    
  
  for (m in 1:length(blocks)) {
    if(is.infinite(means_check[m])) {
      exc <- append(exc, c(a:(a + length(blocks[[m]]) - 1)))
    }
    a = a + length(blocks[[m]])
  }
  return(exc)
}

subset_blocks <- function(blocks, means_check) {  
  for (o in 1:length(blocks)) {
    if(is.infinite(means_check[o])) 
    {blocks[[paste("block", o, sep = "")]] <- NULL}
  }
  
  a = 1
  for (q in 1:length(blocks)) {
    blocks[[q]] <- a:(a + length(blocks[[q]]) - 1)
    a <- a + length(blocks[[q]])
  }
  return(blocks)
}
#  End 3. ---------------------------------------------------------------------

# 4. sample reproduction ------------------------------------------------------
# sample_reproduction <- function(n_lrns = 12L, n_datsets = 1L, seed = 1234L) {
#   set.seed(seed) 
#   lrns <- sort(sample(1:12, n_lrns))
#   datsets <- sort(sample(1:18, n_datsets))
#   
#   list("lrns" = lrns, "datsets" = datsets)
# }
#  End 4. ---------------------------------------------------------------------





