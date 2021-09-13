##############################################
#### Functions for class("mbnma") ####
##############################################




#' Rank parameters from a time-course MBNMA
#'
#' Ranks desired parameters saved from a time-course MBNMA model from "best" to
#' "worst".
#'
#' @inheritParams plot.mb.predict
#' @inheritParams stats::integrate
#' @param params A character vector containing any model parameters monitored
#'   in `mbnma` for which ranking is desired (e.g. `"beta.1"`, `"emax"`).
#'   Parameters must vary by treatment for ranking to be possible. Can include
#'   `"auc"` (see details).
#' @param treats A character vector of treatment/class names (depending on the value of `level`) or
#'   a numeric vector of treatment/class codes (as coded in `mbnma`)
#'   that indicate which treatments/classes to calculate rankings for. If left `NULL``
#'   then rankings will be calculated for all treatments/classes.
#' @param lower_better Indicates whether negative responses are better (`lower_better=TRUE`) or
#'   positive responses are better (`lower_better=FALSE`)
#' @param int.range A numeric vector with two elements that indicates the range
#'   over which to calculate AUC. Takes the form c(lower bound, upper bound). If left
#'   as `NULL` (the default) then the range will be between zero and the maximum follow-up
#'   time in the dataset.
#' @param level A character object to indicate whether the parameters to be ranked are at the treatment
#'   level (`"treatment"`) or class level (`"class"`).
#' @param n.iter The number of iterations for which to calculate AUC (if `"auc"` is included in `params`).
#'   Must be a positive integer. Default is the value used in `mbnma`.
#' @param ... Arguments to be sent to `integrate()`
#'
#' @return A named list whose elements correspond to parameters given in
#'   `params`. Each element contains:
#'   * `summary.rank` A data frame containing
#'   mean, sd, and quantiles for the ranks of each treatment given in `treats`
#'   * `prob.matrix` A matrix of the proportions of MCMC results for which each
#'   treatment in `treats` ranked in which position for the given parameter
#'   * `rank.matrix` A matrix of the ranks of MCMC results for each treatment in
#'   `treats` for the given parameter.
#'
#' @details `"auc"` can be included in `params` to rank treatments based on
#'   Area Under the Curve (AUC). This accounts for the effect of multiple
#'   time-course parameters simultaneously on the treatment response, but will
#'   be impacted by the range of time over which AUC is calculated (`int.range`).
#'   This requires integration over `int.range` and can take some time to run (particularly)
#'   for spline functions as this uses the trapezoid method rather than adaptive quadrature).
#'
#'   As with other post-estimation functions, `rank()` should only be performed on
#'   models which have successfully converged. Note that rankings can be very sensitive to
#'   even small changes in treatment effects and therefore failure to converge in only
#'   one parameter may have substantial impact on rankings.
#'
#' @examples
#' \donttest{
#' # Create an mb.network object from a dataset
#' network <- mb.network(alog_pcfb)
#'
#' # Run an MBNMA model with an Emax time-course
#' emax <- mb.run(network,
#'   fun=temax(pool.emax="rel", method.emax="common",
#'             pool.et50="rel", method.et50="random"),
#'   intercept=FALSE)
#'
#' # Rank treatments by time-course parameter from the model with lower scores being better
#' rank(emax, params=c("emax", "et50"), lower_better=TRUE)
#'
#' # Rank treatments 1-3 by AUC
#' rank(emax, params="auc", treats=c(1:3), lower_better=TRUE,
#'   int.range=c(0,20))
#' }
#'
#' @export
rank.mbnma <- function(x, params="auc", lower_better=FALSE, treats=NULL,
                       int.range=NULL,
                       level="treatment", n.iter=x$BUGSoutput$n.sims,
                       ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mbnma", add=argcheck)
  checkmate::assertCharacter(params, any.missing=FALSE, unique=TRUE, add=argcheck)
  checkmate::assertLogical(lower_better, null.ok=FALSE, len=1, add=argcheck)
  checkmate::assertNumeric(int.range, lower=0, finite=TRUE, any.missing=FALSE, len=2, null.ok=TRUE,
                           sorted=TRUE, add=argcheck)
  checkmate::assertChoice(level, choices=c("treatment", "class"), add=argcheck)
  checkmate::assertInt(n.iter, lower=1, upper=x$BUGSoutput$n.sims, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check level
  if (level=="class" & !("classes" %in% names(x$network))) {
    stop("`level` has been specified as `class` yet `x` is not a class effect model")
  }
  level <- ifelse(level=="treatment", "treatments", "classes")

  # Ensure AUC is the last estimate to be called
  if ("auc" %in% params) {
    params <- c(params[params!="auc"], "auc")

    # if (length(x$model.arg$class.effect)>0) {
    #   stop("AUC cannot currently be calculated for class effect models")
    # }
  }

  # If treats have not been specified then select all of them
  if (is.null(treats)) {
    treats <- x$network[[level]]
  } else if (!is.null(treats)) {
    if (is.character(treats)) {
      if (!all(treats %in% x$network[[level]])) {
        stop("`treats` includes treatments/classes not included in `x`")
      }
    } else if (is.numeric(treats)) {
      if (any(treats > x[["model"]][["data"]]()[["NT"]] | any(treats<1))) {
        stop("If given as numeric treatment/class codes, `treats` must be numbered similarly to treatment/class codes in `x`")
      }
      treats <- x$network[[level]][treats]
    }
  }

  # Provide int.range default values
  if ("auc" %in% params) {

    if (is.null(int.range)) {

      treatsnum <- which(x$network[[level]] %in% treats)
      fupdata <- x$model.arg$jagsdata

      int.max <- max(fupdata$time,
                     na.rm = TRUE)

      int.range <- c(0, int.max)
    }
  }

  # Change beta to d (if present) so that it is identified in mcmc output
  for (i in 1:4) {
    if (paste0("beta.",i) %in% params) {
      if (!paste0("beta.",i) %in% x[["parameters.to.save"]]) {
        params[which(params==paste0("beta.",i))] <- paste0("d.",i)
      }
    }
  }

  rank.result <- list()
  for (i in seq_along(params)) {
    if (params[i] %in% x[["parameters.to.save"]]) {
      param.mod <- x[["BUGSoutput"]][["sims.list"]][[params[i]]]

      # Check that selected parameter is different over multiple treatments
      if (!is.matrix(param.mod) | ncol(param.mod)<=1) {
        msg <- paste0(params[i], " does not vary by treatment and therefore cannot be ranked by treatment")
        stop(msg)
      }

      param.mod <- param.mod[,which(x$network[[level]] %in% treats)]

      rank.mat <- t(apply(param.mod, MARGIN=1, FUN=function(x) {
        order(order(x, decreasing = !lower_better), decreasing=FALSE)
      }))
      colnames(rank.mat) <- treats

      # Ranking probabilityes
      prob.mat <- calcprob(rank.mat, treats=treats)

      # Calculate cumulative ranking probabilities
      cum.mat <- apply(prob.mat, MARGIN=2,
                       FUN=function(col) {cumsum(col)})

      rank.result[[params[i]]] <-
        list("summary"=sumrank(rank.mat),
             "prob.matrix"=prob.mat,
             "rank.matrix"=rank.mat,
             "cum.matrix"=cum.mat,
             "lower_better"=lower_better
             )

    } else if (params[i]=="auc") {

      auc <- rankauc(x, lower_better=lower_better,
                     treats=treats, level=level,
                     int.range=int.range, n.iter=n.iter, ...)

      # Calculate cumulative ranking probs for aux
      auc[["cum.matrix"]] <- apply(auc$prob.matrix, MARGIN=2,
                                   FUN=function(col) {cumsum(col)})
      auc[["lower_better"]] <- lower_better

      rank.result[["auc"]] <- auc
    } else {
      stop(paste0(params[i],
                  " is not a valid paramter saved from the MBNMA model"))
    }
  }

  class(rank.result) <- "mb.rank"
  return(rank.result)

}




#' Print summary MBNMA results to the console
#'
#' @inheritParams predict.mbnma
#' @param ... further arguments passed to `knitr::kable`
#'
#' @export
summary.mbnma <- function(object, ...) {
  checkmate::assertClass(object, "mbnma")

  # State that function does not work if "parameters.to.save" has been specified
  if (!is.null(object$model.arg$parameters.to.save)) {
    stop("Cannot use `summary()` method if `parameters.to.save` have been assigned. Use `print()` instead.")
  }

  # Overall section
  print.overall.str(object)
  cat("\n\n")

  # Print treatment-level section
  print.treat.str(object, ...)

  # Class-effect section
  print.class.str(object, ...)

  # Correlation section
  print.cor.str(object, ...)

  # Model fit statistics section
  cat(print.modfit.str(object))

  # Check for rhat < 1.02
  cat("\n")
  rhat.warning(object)
}





#' Predict responses over time in a given population based on MBNMA time-course
#' models
#'
#' Used to predict responses over time for different treatments or to predict
#' the results of a new study. For MBNMA models that include consistency
#' relative effects on time-course parameters, this is calculated by combining
#' relative treatment effects with a given reference treatment response
#' (specific to the population of interest).
#'
#' @param object An S3 object of `class("mbnma")` generated by running
#'   a time-course MBNMA model
#' @param times A sequence of positive numbers indicating which time points to
#'   predict mean responses for
#' @param E0 An object to indicate the value(s) to use for the response at time = 0
#'   in the prediction. This can take a number of different formats depending
#'   on how it will be used/calculated. The default is 0 but this may lead
#'   to non-sensical predictions.
#'   * `numeric()` A single numeric value representing the deterministic response at time = 0
#'   * `formula()` A formula representing a stochastic distribution for the response
#'   at time = 0. This is specified as a random number generator
#'   (RNG) given as a string, and can take any RNG distribution for which a function exists
#'   in R. For example: `~rnorm(n, 7, 0.5)`.
#' @param treats A character vector of treatment/class names or a numeric vector of treatment/class codes (as coded
#'   in `mbnma`) that indicates which treatments/classes to calculate predictions for. If left as `NULL` then
#'   predictions will be calculated for all treatments/classes. Whether the vector should correspond to treatments or
#'   classes depends on the value of `level`.
#' @param level Can take either `"treatment"` to make predictions for treatments, or `"class"` to make predictions for classes (in
#'   which case `object` must be a class effect model).
#' @param ref.resp An object to indicate the value(s) to use for the reference treatment response in MBNMA models
#'   in which the reference treatment response is not estimated within the model (i.e. those that model any time-
#'   course parameters using `pool="rel"`). This can take a number of different formats depending
#'   on how it will be used/calculated. There are two approaches for this:
#'
#'   1. The reference response can be estimated from a dataset of studies investigating the reference
#'   treatment using meta-analysis. This dataset could be a set of observational
#'   studies that are specific to the population on which to make
#'   predictions, or it could be a subset of the study arms within the MBNMA dataset
#'   that investigate the reference treatment. The data should be provided to `ref.resp` as a
#'   `data.frame()` containing the data in long format (one row per observation). See [ref.synth()]
#'
#'   2. Values for the reference treatment response can be assigned to different time-course parameters
#'   within the model that have been modelled using consistency relative effects (`pool="rel"`).
#'   These are given as a list, in which each named element corresponds to a time-course
#'   parameter modelled in `mbnma`. Their values can be either of the following:
#'   * `numeric()` A numeric value representing the deterministic value of the time-course parameter in
#'   question in individuals given the reference treatment. `0` is used as the default, which assumes no
#'   effect of time on the reference treatment.
#'   * `formula()` A formula representing a stochastic distribution for the value of the time-course
#'   parameter in question. This is specified as a random number generator (RNG) given as a formula,
#'   and can take any RNG distribution for which a function exists in R. For example: `~rnorm(n, -3, 0.2)`.

#' @param synth A character object that can take the value `"common"` or `"random"` that
#'   specifies the the type of pooling to use for synthesis of `ref.resp`. Using `"random"` rather
#'   than `"common"` for `synth` will result in wider 95\\% CrI for predictions.
#' @param ... Arguments to be sent to R2jags for synthesis of the network
#'   reference treatment effect (using [ref.synth()])
#'
#'
#' @return An S3 object of class `mb.predict` that contains the following
#'   elements:
#'   * `summary` A named list of data frames. Each data frame contains
#'   a summary of predicted responses at follow-up times specified in `times`
#'   for each treatment specified in `treats`
#'   * `pred.mat` A named list of
#'   matrices. Each matrix contains the MCMC results of predicted responses at
#'   follow-up times specified in `times` for each treatment specified in
#'   `treats`
#'
#' @details `ref.resp` only needs to be specified if `mbnma` has
#'   been estimated using consistency relative effects (`pool="rel"`) for
#'   any time-course parameters, as these inform the absolute values of the
#'   network reference treatment parameters which can then be added to the
#'   relative effects to calculate specific predictions.
#'
#' @examples
#' \donttest{
#' # Create an mb.network object from a dataset
#' network <- mb.network(osteopain)
#'
#' # Run an MBNMA model with an Emax time-course
#' emax <- mb.run(network,
#'   fun=temax(pool.emax="rel", method.emax="common",
#'     pool.et50="abs", method.et50="common"))
#'
#' # Predict responses using a stochastic baseline (E0) and a distribution for the
#' #network reference treatment
#' preds <- predict(emax, times=c(0:10),
#'   E0=~rnorm(n, 7, 0.5),
#'   ref.resp=list(emax=~rnorm(n, -0.5, 0.05)))
#' summary(preds)
#'
#' # Predict responses using the original dataset to estimate the network reference
#' #treatment response
#' paindata.ref <- osteopain[osteopain$treatname=="Placebo_0",]
#' preds <- predict(emax, times=c(5:15),
#'   E0=10,
#'   ref.resp=paindata.ref)
#' summary(preds)
#'
#' # Repeat the above prediction but using a random effects meta-analysis of the
#' #network reference treatment response
#' preds <- predict(emax, times=c(5:15),
#'   E0=10,
#'   ref.resp=paindata.ref,
#'   synth="random")
#' summary(preds)
#' }
#'
#' @export
predict.mbnma <- function(object, times=seq(0, max(object$model.arg$jagsdata$time, na.rm=TRUE), length.out=30),
                          E0=0,
                          treats = NULL, level="treatment",
                          ref.resp=NULL, synth="common",
                          ...) {
  ######## CHECKS ########

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(object, "mbnma", add=argcheck)
  checkmate::assertNumeric(times, lower=0, finite=TRUE, any.missing=FALSE, unique=TRUE,
                           sorted=TRUE, add=argcheck)
  checkmate::assertChoice(level, choices=c("treatment", "class"), add=argcheck)
  checkmate::assertChoice(synth, choices=c("random", "common"), add=argcheck)
  #checkmate::assertClass(treats, classes=c("numeric", "character"), null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (object$model.arg$link=="rom") {
    stop("'predict()' cannot currently be used with MBNMAs modelled using ratios of means (link='rom')")
  }

  # Check if level="class" that class effect model was fitted
  if (level=="class") {
    if (length(object[["model.arg"]][["class.effect"]])==0) {
      stop(crayon::red(crayon::bold("`level` has been set to `class` but no class effect models were not used")))
    }
    if (!isTRUE(all.equal(
      object$model.arg$fun$params[object$model.arg$fun$apool %in% "rel"],
      names(object$model.arg$class.effect)
    ))) {
      stop(crayon::red(crayon::bold("To predict level='class' all relative effects must be modelled with class effects")))
    }
    level <- "classes"
  } else if (level=="treatment") {
    level <- "treatments"
  }

  # Check whether UME has been used and stop if so
  if (object[["model.arg"]][["UME"]]!=FALSE) {
    stop(crayon::red(crayon::bold("UME model cannot be used for prediction")))
  }

  # Check ref.resp has been specified correctly
  if ("rel" %in% object$model.arg$fun$apool) {
    if (is.null(ref.resp)) {

      # If ref.resp is not given then assign 0 to all rel time-course parameters
      ref.resp <- list()
      rels <- names(object$model.arg$fun$apool)[object$model.arg$fun$apool %in% "rel"]
      for (i in seq_along(rels)) {
        ref.resp[[rels[i]]] <- 0
      }

      # # If ref.resp is not given then assign mbnma MCMC value to all abs time-course parameters
      # abs <- names(object$model.arg$fun$apool)[object$model.arg$fun$apool %in% "abs"]
      # for (i in seq_along(abs)) {
      #   ref.resp[[abs[i]]] <- mbnma$BUGSoutput$sims.list[[abs[i]]]
      # }
    } else {

      # If ref.resp is given ensure it is of the correct class
      if (!(any(class(ref.resp) %in% c("data.frame", "tibble", "list")))) {
        stop(crayon::red(crayon::bold("`object` includes time-course parameters modelled using relative effects (pool=`rel`).
      The reference treatment response for them must be provided to `ref.resp` as a list,
      or estimated from a dataset of reference treatment studies by providing a data frame.")))
      }
    }
  }


  # If treats have not been specified then select all of them
  NT <- ifelse(level=="treatments", object$model.arg$jagsdata$NT, object$model.arg$jagsdata$Nclass)
  if (is.null(treats)) {
    #treats <- c(1:object[["model"]][["data"]]()[["NT"]])
    treats <- object$network[[level]]
  } else if (!is.null(treats)) {
    if (is.numeric(treats)) {
      if (any(treats > NT | any(treats<1))) {
        stop(crayon::red(crayon::bold("If given as numeric treatment/class codes, `treats` must be numbered similarly to treatment/class codes in `object`")))
      }
      treats <- object$network[[level]][treats]
    }
    if (is.character(treats)) {
      if (!all(treats %in% object$network[[level]])) {
        stop(crayon::red(crayon::bold("`treats` includes treatments/classes not included in `object`")))
      }
    }
  }

  #### Check E0 ####
  if (is.null(E0)) {
    stop(crayon::red(crayon::bold("E0 has not been defined")))
  }

  # Check that distribution for E0 is of the correct format
  if (class(E0)=="formula") {
    E0 <- as.character(E0)[2]
    if (grepl("r[A-z]+\\(n,.+\\)", E0)==FALSE) {
      stop(crayon::red(crayon::bold("Stochastic distribution for E0 must be expressed as a string in the form of a supported R distribution (e.g. '~rnorm(n, 5,2)')")))
    }
  } else if (is.numeric(E0)) {
    if (length(E0)!=1) {
      stop(crayon::red(crayon::bold("`E0` can only take a single numeric value if not expressed as a stochastic distribution")))
    }
  } else {
    stop(crayon::red(crayon::bold("'E0' has been incorrectly specified")))
  }


  ###### Extract info from mbnma #######

  n <- object$BUGSoutput$n.sims

  # Initial predict parameters
  # timecourse <- list(init.predict(object)[["timecourse"]])
  # beta.incl <- init.predict(object)[["beta.incl"]]
  # beta.incl <- object$model.arg$fun$params
  # timecourse <- paste0("alpha + ", object$model.arg$fun$jags)

  # Extract parameter values from MBNMA result
  model.vals <- get.model.vals(mbnma=object, E0=E0, level=level)
  timecourse <- model.vals[["timecourse"]]
  time.params <- model.vals[["time.params"]]


  ########## Get reference treatment effect ###########
  #mu.prior <- model.vals[["mu.prior"]]
  mu.params <- time.params[grepl("^mu\\.", time.params)]



  if (!is.null(ref.resp)) {

    # If ref.resp specified as values for each time-course parameter (in a list)
    if (any(class(ref.resp)=="list")) {
      msg <- paste0("Priors required for: ", paste(mu.params, collapse=", "))
      message(msg)

      names(ref.resp) <- paste0("mu.", match(names(ref.resp), object$model.arg$fun$params))

      if (identical(sort(mu.params), sort(names(ref.resp)))==FALSE) {
        msg <- "Named elements of `ref.resp` do not correspond to consistency time-course parameters monitored within the model."
        stop(crayon::bold(crayon::red(msg)))
      } else {
        message(crayon::bold(crayon::green("Success: Elements in prior match consistency time-course treatment effect parameters")))
      }

      # Assign ref.resp to mu values in model
      for (i in seq_along(ref.resp)) {

        if (class(ref.resp[[i]]) %in% c("formula", "character")) {

          if (class(ref.resp[[i]]) %in% "formula") {
            ref.resp[[i]] <- as.character(ref.resp[[i]])[2]
            if (grepl("r[A-z]+\\(n,.+\\)", ref.resp[[i]])==FALSE) {
              stop(crayon::red("Stochastic distribution for ref.resp must be expressed as a formula in the form of a supported R distribution (e.g. ~rnorm(n, 5,2))"))
            }
          }
          assign(mu.params[which(names(ref.resp)[i]==mu.params)],
                 eval(parse(text=ref.resp[[i]])))

        } else if (class(ref.resp[[i]]) %in% "numeric") {
          assign(mu.params[which(names(ref.resp)[i]==mu.params)],
                 rep(ref.resp[[i]], n))
        }
      }
    } else if (any(class(ref.resp) %in% c("data.frame", "tibble"))) {

      ### PLACEBO SYNTHESIS MODEL ###
      args <- list(...)
      synth.result <- do.call(ref.synth, args=c(list(data.ab=ref.resp, mbnma=object, synth=synth), args))

      # synth.result <- ref.synth(data.ab=ref.resp, mbnma=object, synth=synth, ...)

      synth.result <- synth.result$BUGSoutput$median
      synth.result[["deviance"]] <- NULL

      # Assign synth.result to mu values in model
      for (i in seq_along(mu.params)) {
        if (synth=="random") {
          assign(mu.params[i],
                 stats::rnorm(n,
                              synth.result[[mu.params[i]]],
                              synth.result[[paste0("sd.", mu.params[i])]])
          )
        } else if (synth=="common") {
          assign(mu.params[i], rep(synth.result[[mu.params[i]]],
                                   n))
        } else (stop(crayon::red("synth must be either `common` or `random`")))

      }
    }
  }

  # Convert predicted times to splines
  if (any(c("rcs", "bs", "ns", "ls") %in% object$model.arg$fun)) {
    timecourse <- gsub("\\[i\\,", "[", timecourse)
    spline <- genspline(times, spline=object$model.arg$fun$name, knots=object$model.arg$fun$knots, degree=object$model.arg$fun$degree)
  }

  ########## Predict responses ###########

  # Assign E0 to alpha in model
  #alpha <- eval(parse(text=E0)) # TO BE REMOVED
  alpha <- model.vals$alpha


  beta.params <- time.params[grepl("^beta.", time.params)]
  # Assign single beta results to beta values in model
  for (i in seq_along(beta.params)) {
    if (!is.matrix(model.vals[[beta.params[i]]])) {
      assign(beta.params[i], model.vals[[beta.params[i]]])
    } else if (is.matrix(model.vals[[beta.params[i]]])) {
      if (ncol(model.vals[[beta.params[i]]])==1) {
        assign(beta.params[i], model.vals[[beta.params[i]]])
      }
    }
  }

  d.params <- time.params[grepl("^d\\.", time.params)]

  predicts <- list()
  treatsnum <- which(object$network[[level]] %in% treats)
  for (treat in seq_along(treatsnum)) {

    # Assign d results to d values in model
    for (i in seq_along(d.params)) {
      assign(d.params[i], model.vals[[d.params[i]]][,treatsnum[treat]])
    }

    treatpred <- data.frame("pred"=rep(NA,n))
    for (m in seq_along(times)) {
      time <- times[m]

      # Evaluate function
      pred <- eval(parse(text=timecourse))

      if (any(is.na(pred))) {
        pred[is.na(pred)] <- 0
      }

      #treatpred <- cbind(treatpred, pred)
      # if (is.vector(pred)) {
      #   pred <- as.matrix(pred, ncol=1)
      # }

      treatpred[[paste0("time", times[m])]] <- pred

    }

    # predicts[[paste0(treats[treat])]] <- treatpred[,1]
    predicts[[paste0(treats[treat])]] <- treatpred[-1]
  }

  # Generate summary data frame
  sumpred <- list()
  for (i in seq_along(treats)) {
    summary <- data.frame("time"=times)

    summary[["mean"]] <- apply(predicts[[as.character(treats[i])]], MARGIN=2,
                               FUN=function(x) mean(x))
    summary[["sd"]] <- apply(predicts[[as.character(treats[i])]], MARGIN=2,
                             FUN=function(x) stats::sd(x))

    quantiles <- apply(predicts[[as.character(treats[i])]], MARGIN = 2,
                       function(x) stats::quantile(x,
                                                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
    summary <- cbind(summary, t(quantiles))

    sumpred[[as.character(treats[i])]] <- summary
  }

  #predict.result <- list("summary"=sumpred, "pred.mat"=predicts, "treatments"=object$network$treatments, "mbnma"=object)
  predict.result <- list("summary"=sumpred, "pred.mat"=predicts, "network"=object$network,
                         "times"=times, "link"=object$model.arg$link)
  class(predict.result) <- "mb.predict"

  return(predict.result)
}







#' Forest plot for results from time-course MBNMA models
#'
#' Generates a forest plot for time-course parameters of interest from results from time-course MBNMA models.
#' Posterior densities are plotted above each result using `ggdist:stat_:halfeye()`
#'
#' @param x An S3 object of class `"mbnma"` generated by running
#' a time-course MBNMA model
#' @param params A character vector of time-course parameters to plot.
#' Parameters must be given the same name as monitored nodes in `mbnma` and must vary by treatment or class. Can be set to
#' `NULL` to include all available time-course parameters estimated by `mbnma`.
#' @param treat.labs A character vector of treatment labels. If left as `NULL` (the default) then
#' labels will be used as defined in the data.
#' @param class.labs A character vector of class labels if `mbnma` was modelled using class effects
#' If left as `NULL` (the default) then labels will be used as defined in the data.
#' @param ... Arguments to be sent to `ggdist::stat_halfeye()`
#'
#' @return A forest plot of class `c("gg", "ggplot")` that has separate panels for different time-course parameters
#'
#' @examples
#'\donttest{
#' # Create an mb.network object from a dataset
#' alognet <- mb.network(alog_pcfb)
#'
#' # Run an MBNMA model with an Emax time-course
#' emax <- mb.run(alognet,
#'   fun=temax(pool.emax="rel", method.emax="common",
#'     pool.et50="rel", method.et50="common"),
#'   intercept=FALSE)
#'
#' # Generate forest plot
#' plot(emax)
#'
#' # Plot results for only one time-course parameter
#' plot(emax, params="emax")
#' }
#' @export
plot.mbnma <- function(x, params=NULL, treat.labs=NULL, class.labs=NULL, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mbnma", add=argcheck)
  checkmate::assertCharacter(params, null.ok=TRUE, add=argcheck)
  checkmate::assertCharacter(treat.labs, null.ok=TRUE, add=argcheck)
  checkmate::assertCharacter(class.labs, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Declare global variables
  timeparam <- NULL
  Var2 <- NULL
  value <- NULL
  ndistinct <- NULL

  # Change beta to d (if present) so that it is identified in mcmc output
  for (i in 1:4) {
    if (paste0("beta.",i) %in% params) {
      if (!paste0("beta.",i) %in% x[["parameters.to.save"]]) {
        params[which(params==paste0("beta.",i))] <- paste0("d.",i)
      }
    }
  }

  # Check that specified params are monitored in model
  if (!all(params %in% x[["parameters.to.save"]])) {
    stop(paste0("Variable 'params': Must contain elements of set {", paste(x[["parameters.to.save"]], collapse = ", "), "}"))
  }

  # Check that specified params vary over treatment
  for (i in seq_along(params)) {
    if (length(x[["BUGSoutput"]][["median"]][[params[i]]])<=1) {
      stop(paste0(params[i], " does not vary over treatment/class and cannot be plotted"))
    }
  }


  # Add all available params if is.null(params)
  if (is.null(params)) {

    params <- c(x$model.arg$fun$params,
                paste0("sd.", x$model.arg$fun$params),
                toupper(x$model.arg$fun$params),
                paste0("sd.", toupper(x$model.arg$fun$params)),
                paste0("d.", 1:4),
                paste0("sd.beta.", 1:4),
                paste0("D.", 1:4),
                paste0("sd.BETA.", 1:4)
    )

    params <- params[params %in% x$parameters.to.save]

    # Selects only parameters which vary by class or agent (with more than 1 element)
    drop <- vector()
    for (i in seq_along(params)) {
      if (length(x[["BUGSoutput"]][["median"]][[params[i]]])<=1) {
        drop <- append(drop, i)
      }
    }
    if (length(drop)>0) {
      params <- params[-drop]
    }

    if (length(params)==0) {
      stop("No time-course consistency parameters can be identified from the model")
    }
  }

  # # Compile parameter data into one data frame
  # mb.sum <- as.data.frame(x[["BUGSoutput"]][["summary"]])
  # plotdata <- mb.sum[0,]
  # for (i in seq_along(params)) {
  #   paramdata <- mb.sum[grepl(paste0("^", params[i]),rownames(mb.sum)),]
  #   paramdata[["timeparam"]] <- rep(params[i], nrow(paramdata))
  #   plotdata <- rbind(plotdata, paramdata)
  # }
  # plotdata[["param"]] <- as.numeric(gsub("(.+\\[)([0-9]+)(\\])", "\\2", rownames(plotdata)))

  # Compile parameter data into one data frame
  mcmc <- x$BUGSoutput$sims.list
  plotdata <- data.frame(Var2=NA, value=NA, timeparam=NA)
  for (i in seq_along(params)) {
    paramdata <- reshape2::melt(mcmc[[params[i]]])[,2:3]
    paramdata$timeparam <- rep(params[i], nrow(paramdata))
    plotdata <- rbind(plotdata, paramdata)
  }
  plotdata <- plotdata[-1,]

  # Change param labels for treatments
  treatdat <- plotdata[!grepl("[[:upper:]]", plotdata$timeparam),]
  if (!is.null(treat.labs)) {
    treatcodes <- treatdat$Var2
    if (length(treat.labs)!=max(treatcodes)) {
      stop("`treat.labs` length does not equal number of treatments that have been modelled for this time-course parameter")
    } else {
      t.labs <- treat.labs[sort(unique(treatcodes))]
    }
  } else if ("treatments" %in% names(x$network)) {
    t.labs <- x$network[["treatments"]]
  } else {
    t.labs <- sort(unique(treatdat$param))
  }

  # Change param labels for classes
  classdat <- plotdata[grepl("[[:upper:]]", plotdata$timeparam),]
  c.labs <- vector()
  if (nrow(classdat)!=0) {
    if (!is.null(class.labs)) {
      classcodes <- classdat$Var2
      c.labs <- class.labs[classcodes]
    } else if ("classes" %in% names(x$network)) {
      c.labs <- x$network[["classes"]][x$network[["classes"]]!="Placebo"]
    } else {
      c.labs <- sort(unique(classdat$param))
    }
  }

  # Increase param number (Var2) for classes
  ntreat <- ifelse(nrow(treatdat)>0, max(treatdat$Var2), 0)
  plotdata$Var2[grepl("[[:upper:]]", plotdata$timeparam)] <-
    plotdata$Var2[grepl("[[:upper:]]", plotdata$timeparam)] + ntreat

  # Attach labels
  if (nrow(treatdat)>0) {
    all.labs <- c(t.labs, c.labs)
  } else {all.labs <- c.labs}
  plotdata$Var2 <- factor(plotdata$Var2, labels=all.labs)

  if (any(is.na(levels(plotdata$param)))) {
    stop("`treat.labs` or `class.labs` have not been specified correctly")
  }

  # Drop paramers fixed to zero (i.e. network reference)
  plotdata <- plotdata %>% dplyr::group_by(timeparam, Var2) %>%
    dplyr::mutate(ndistinct=dplyr::n_distinct(value)) %>%
    subset(ndistinct>1)

  # Create forest plot
  g <- ggplot2::ggplot(plotdata, ggplot2::aes(x = value, y = Var2)) +
    ggdist::stat_halfeye() +
    ggplot2::facet_wrap(~timeparam, scales="free")

  # Axis labels
  g <- g + ggplot2::xlab("Treatment / Class") +
    ggplot2::ylab("Effect size") +
    theme_mbnma()

  g <- g + do.call(ggdist::stat_halfeye, args = list(...))

  graphics::plot(g)
  return(invisible(g))
}







#' Calculates relative effects/mean differences at a particular time-point
#'
#' Uses mbnma time-course parameter estimates to calculate treatment
#' differences between treatments or classes at a particular time-point.
#' Can be used to compare treatments evaluated in studies at different follow-up times.
#'
#' @param time A numeric value for the time at which to estimate relative effects/mean differences.
#' @param treats A character vector of treatment names for which to calculate relative effects/mean differences.
#' Must be a subset of `mbnma$network$treatments`
#' @param classes A character vector of class names for which to calculate relative effects/mean differences from.
#' Must be a subset of `mbnma$network$classes`. Only works for class effect models.
#' @inheritParams predict.mbnma
#' @inheritParams fitplot
#'
#' @return An object of class `"relative.array"` list containing:
#' * The time-point for which results are estimated
#' * Matrices of posterior means, medians, SDs and upper and lower 95% credible intervals for the
#' differences between each treatment
#' * An array containing MCMC results for the differences between all treatments specified in `treats`
#' or all classes specified in `classes`.
#'
#' Results are reported in tables as the row-defined treatment minus the column-defined treatment.
#'
#' @examples
#' \donttest{
#' # Create an mb.network object from a dataset
#' alognet <- mb.network(alog_pcfb)
#'
#' # Run a quadratic time-course MBNMA using the alogliptin dataset
#' mbnma <- mb.run(alognet,
#'   fun=tpoly(degree=2,
#'   pool.1="rel", method.1="random",
#'   pool.2="rel", method.2="common"
#'   )
#' )
#'
#' # Calculate differences between all treatments at 20 weeks follow-up
#' allres <- get.relative(mbnma, time=20)
#'
#' # Calculate difference between a subset of treatments at 10 weeks follow-up
#' subres <- get.relative(mbnma, time=10,
#'   treats=c("alog_50", "alog_25", "placebo"))
#' }
#' @export
get.relative <- function(mbnma, time=max(mbnma$model.arg$jagsdata$time, na.rm=TRUE),
                         treats=mbnma$network$treatments, classes=NULL) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(time, lower=0, len=1, null.ok=FALSE, add=argcheck)
  checkmate::assertSubset(treats, choices=mbnma$network$treatments, add=argcheck)
  checkmate::assertSubset(classes, choices=mbnma$network$classes, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (!is.null(classes)) {
    treats <- classes
    level <- "class"
  } else {
    level <- "treatment"
  }

  pred <- suppressMessages(predict.mbnma(mbnma, times=time,
                        treats = treats, level = level))


  # Matrix of results
  mat <- do.call(cbind, pred$pred.mat)

  # For lower triangle
  outmat <- array(dim=c(length(treats), length(treats), nrow(mat)))
  for (i in 1:(ncol(mat)-1)) {
    temp <- mat[,-i] - mat[,i]
    #temp <- apply(temp, MARGIN=2, FUN=function(x) {neatCrI(quantile(x, probs=c(0.025, 0.5, 0.975)), digits = 2)})
    outmat[(1+i):dim(outmat)[1],i,] <- t(temp[,i:ncol(temp)])
  }

  # For upper triangle
  for (i in 1:(ncol(mat)-1)) {
    temp <- mat[,i] - mat[,-i]
    outmat[i,(1+i):dim(outmat)[2],] <- t(temp[,i:ncol(temp)])
  }

  dimnames(outmat)[[1]] <- treats
  dimnames(outmat)[[2]] <- treats



  ######### Summary matrixes ######

  xmat <- outmat

  meanmat <- matrix(nrow=nrow(xmat), ncol=ncol(xmat))
  semat <- meanmat
  medmat <- meanmat
  l95mat <- medmat
  u95mat <- medmat

  for (i in 1:nrow(xmat)) {
    for (k in 1:ncol(xmat)) {
      if (!is.na(xmat[i,k,1])) {
        meanmat[i,k] <- mean(xmat[i,k,])
        semat[i,k] <- stats::sd(xmat[i,k,])
        medmat[i,k] <- stats::median(xmat[i,k,])
        l95mat[i,k] <- stats::quantile(xmat[i,k,], probs = 0.025)
        u95mat[i,k] <- stats::quantile(xmat[i,k,], probs = 0.975)
      }
    }
  }

  sumlist <- list("mean"=meanmat, "se"=semat, "median"=medmat, "lower95"=l95mat, "upper95"=u95mat)

  for (i in seq_along(sumlist)) {
    dimnames(sumlist[[i]])[[1]] <- dimnames(xmat)[[1]]
    dimnames(sumlist[[i]])[[2]] <- dimnames(xmat)[[2]]
  }

  out <- list("time"=time, "relarray"=outmat)
  out <- c(out, sumlist)

  class(out) <- "relative.array"

  return(out)
}
