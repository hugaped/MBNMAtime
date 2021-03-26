# Functions for ranking treatment responses
# Author: Hugo Pedder
# Date created: 2019-07-24


#' Set rank as a method
#'
#' @param x An object on which to apply the rank method
#' @param ... Arguments to be passed to methods
#'
#' @export
rank <- function (x, ...) {
  UseMethod("rank", x)
}


#' Rank parameters from a time-course MBNMA
#'
#' Ranks desired parameters saved from a time-course MBNMA model from "best" to
#' "worst".
#'
#' @inheritParams plot.mb.predict
#' @inheritParams stats::integrate
#' @param params A character vector containing any model parameters monitored
#'   in `mbnma` for which ranking is desired (e.g. `"beta.1"`, `"d.emax"`).
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
#'             pool.et50="rel", method.et50="random"))
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

      rank.result[[params[i]]] <-
        list("summary"=sumrank(rank.mat),
             "prob.matrix"=calcprob(rank.mat, treats=treats),
             "rank.matrix"=rank.mat)

    } else if (params[i]=="auc") {
      rank.result[["auc"]] <- rankauc(x, lower_better=lower_better,
                                       treats=treats, level=level,
                                       int.range=int.range, n.iter=n.iter, ...)
    } else {
      stop(paste0(params[i],
                  " is not a valid paramter saved from the MBNMA model"))
    }
  }

  class(rank.result) <- "mb.rank"
  return(rank.result)

}



#' Calculates ranking probabilities for AUC from a time-course MBNMA
#'
#' @inheritParams predict.mbnma
#' @inheritParams rank.mbnma
#' @param mbnma An S3 object of class `"mbnma"` generated by running
#' a time-course MBNMA model
#' @param treats A character vector of treatment/class names (depending on the value of `level`). If left `NULL``
#'   then rankings will be calculated for all treatments/classes. Note that unlike `rank.mbnma()` this argument
#'   cannot take a numeric vector.
#'
#' @inherit rank.mbnma return
#' @inherit rank.mbnma details
#'
rankauc <- function(mbnma, lower_better=FALSE, treats=NULL, level="treatments",
                    int.range=c(0,max(mbnma$network$data.ab$time)),
                    n.iter=mbnma$BUGSoutput$n.sims, subdivisions=100, ...) {

  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertCharacter(treats, null.ok = FALSE)
  checkmate::assertChoice(level, choices = c("treatments", "classes"), add=argcheck)
  checkmate::assertLogical(lower_better, any.missing=FALSE, len=1, add=argcheck)
  checkmate::assertIntegerish(subdivisions, lower=2, add=argcheck)
  checkmate::assertIntegerish(int.range, lower=0, any.missing=FALSE, len=2, sorted=TRUE,
                              add=argcheck)
  #checkmate::assertInt(subdivisions, lower=1, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Initial predict parameters
  #nsims <- mbnma$BUGSoutput$n.sims
  # timecourse <- init.predict(mbnma)[["timecourse"]]
  # beta.incl <- init.predict(mbnma)[["beta.incl"]]

  # Extract parameter values from MBNMA result
  model.vals <- get.model.vals(mbnma, E0=0,
                               level=level)

  # Create vector of parameters in expanded time-course function
  timecourse <- model.vals[["timecourse"]]
  timecourse <- gsub("i\\.", "", timecourse)
  time.params <- model.vals[["time.params"]]

  # Switch spline in timecourse and generate spline matrix
  if (any(c("rcs", "ns", "bs", "ls") %in% mbnma$model.arg$fun$name)) {
    timecourse <- gsub("\\[m\\,", "[,", timecourse)

    seg <- seq(from=int.range[1], to=int.range[2], length.out=subdivisions)
    spline <- genspline(seg,
                        spline=mbnma$model.arg$fun$name, knots=mbnma$model.arg$fun$knots, degree=mbnma$model.arg$fun$degree)
  }

  # Replace mu with 0
  for (i in seq_along(time.params)) {
    if (grepl("^mu", time.params[i])) {
      timecourse <- gsub(time.params[i], 0, timecourse)
      time.params <- time.params[-i]
    }
  }


  # NEW SECTION
  auc.result <- matrix(ncol=length(treats))
  rank.mat <- matrix(ncol=length(treats))
  pb <- utils::txtProgressBar(0, n.iter, style = 3)
  for (mcmc in 1:n.iter) {
    utils::setTxtProgressBar(pb, mcmc)

    # Initialise variables
    auc <- vector()
    rank <- matrix(nrow=length(treats), ncol=length(treats))

    treatsnum <- which(mbnma$network[[level]] %in% treats)
    for (treat in seq_along(treatsnum)) {
      time.mcmc <- timecourse

      #for (i in seq_along(params)) {
      for (i in seq_along(time.params)) {
        # Replace parameter in time-course with value for given treat and MCMC
        #timecourse <- gsub(names(params)[i], params[[i]][mcmc], timecourse)
        temp <- model.vals[[time.params[i]]]
        if (is.vector(temp)) {temp <- matrix(temp, ncol=1)}
        time.mcmc <- gsub(time.params[i],
                          ifelse(is.matrix(temp) & ncol(temp)>1, temp[mcmc,treatsnum[treat]], temp[mcmc]),
                          time.mcmc)
      }

      if (any(c("rcs", "ns", "bs", "ls") %in% mbnma$model.arg$fun$name)) {

        # Using trapezoid method for spline function
        y <- eval(parse(text=time.mcmc))
        auc <- append(auc, sum(diff(seg)*zoo::rollmean(y,2)))

      } else {
        temp <- paste("int.fun <- function(time) {",
                      time.mcmc,
                      "}",
                      sep=" ")

        eval(parse(text=temp))

        integral <- stats::integrate(int.fun,
                                     lower=int.range[1], upper=int.range[2],
                                     subdivisions=subdivisions)

        auc <- append(auc, integral$value)
      }

    }

    # Need to name columns with treats
    auc.result <- rbind(auc.result, auc)

    # Ranking
    #rank.mat <- rbind(rank.mat, order(auc, decreasing = decreasing))
  }
  auc.result <- auc.result[-1,]
  #rank.mat <- rank.mat[-1,]

  # Ranking
  rank.mat <- t(apply(auc.result, MARGIN=1, FUN=function(x) {
    order(order(x, decreasing = !lower_better), decreasing=FALSE)
  }))

  colnames(auc.result) <- treats
  colnames(rank.mat) <- treats

  summary.rank <- sumrank(rank.mat)

  return(list("summary"=summary.rank,
              "prob.matrix"=calcprob(rank.mat, treats=treats),
              "rank.matrix"=rank.mat,
              "auc.int"=auc.result
  ))
}




sumrank <- function(rank.mat) {
  if (is.null(colnames(rank.mat))) {
    colnames(rank.mat) <- c(1:ncol(rank.mat))
  }

  quantiles.rank <- apply(X=rank.mat, MARGIN = 2,
                          function(x) stats::quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
  summary.rank <- data.frame(
    "treatment"=colnames(rank.mat),
    "mean"= apply(X=rank.mat, MARGIN = 2, function(x) {base::mean(x)}),
    "sd"= apply(X=rank.mat, MARGIN = 2, function(x) {stats::sd(x)})
  )
  summary.rank <- cbind(summary.rank, t(quantiles.rank))
  rownames(summary.rank) <- NULL

  return(summary.rank)
}


#' Calculates a matrix of ranking probabilities from a matrix of treatment/agent/class
#' rankings
#'
#' @noRd
calcprob <- function(rank.mat, treats=NULL) {
  NT <- ncol(rank.mat)
  rank.prob <- vector(length=NT)

  for (c in 1:NT) {
    pos.vec <- vector()
    for (r in 1:NT) {
      pos.vec <- append(pos.vec,
                        length(rank.mat[rank.mat[,c]==r,c])/nrow(rank.mat))
    }
    rank.prob <- cbind(rank.prob, pos.vec)
  }
  rank.prob <- rank.prob[,-1]

  if (!is.null(treats)) {
    colnames(rank.prob) <- treats
  }

  return("rank.prob"=rank.prob)
}







#' Rank predictions at a specific time point
#'
#' @param x an object of `class("mb.predict")` that contains predictions from an MBNMA model
#' @param time a number indicating the time point at which predictions should be ranked. It must
#' be one of the time points for which predictions in `x` are available.
#' @param treats A character vector of treatment/class names for which responses have been predicted
#'   in `x` As default, rankings will be calculated for all treatments/classes in `x`.
#' @param lower_better Indicates whether negative responses are better (`lower_better=TRUE`) or
#'
#' @inheritParams rank.mbnma
#'
#' @export
rank.mb.predict <- function(x, time=max(x$summary[[1]]$time), lower_better=FALSE,
                                        treats=names(x$summary)) {

  # Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mb.predict", add=argcheck)
  checkmate::assertNumeric(time, len=1, lower=0, add=argcheck)
  checkmate::assertLogical(lower_better, null.ok=FALSE, len=1, add=argcheck)
  checkmate::assertCharacter(treats, null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (!all(treats %in% names(x$summary))) {
    stop("'treats' includes treatments/classes not included in 'x'")
  }
  if (!time %in% x$summary[[1]]$time) {
    stop("'time' is not a time point given in 'x'")
  }


  #### Compute rankings ####

  # Get time point of interest from x
  index <- which(x$summary[[1]]$time %in% time)
  rank.mat <- lapply(x$pred.mat, FUN=function(k) k[,index])
  rank.mat <- t(do.call(rbind, rank.mat))

  # Compute rankings for each iteration
  rank.mat <- t(apply(rank.mat, MARGIN=1, FUN=function(x) {
    order(order(x, decreasing = !lower_better), decreasing=FALSE)
  }))
  colnames(rank.mat) <- treats

  # Store rankings
  rank.result <- list(temp=
                        list("summary"=sumrank(rank.mat),
                             "prob.matrix"=calcprob(rank.mat, treats=treats),
                             "rank.matrix"=rank.mat)
                        )
  names(rank.result) <- paste0("Predictions at time = ", time)

  class(rank.result) <- "mb.rank"
  return(rank.result)
}
