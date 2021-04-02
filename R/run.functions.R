# Functions for running MBNMA models
# Author: Hugo Pedder
# Date created: 2018-09-10


#' Run MBNMA time-course models
#'
#' Fits a Bayesian time-course model for model-based network meta-analysis
#' (MBNMA) that can account for repeated measures over time within studies by
#' applying a desired time-course function. Follows the methods of \insertCite{pedder2019;textual}{MBNMAtime}.
#'
#' @param network An object of class `"mb.network"`.
#' @param parameters.to.save A character vector containing names of parameters
#'   to monitor in JAGS
#' @param fun An object of class `"timefun"` generated (see Details) using any of
#'   `tloglin()`, `tpoly()`, `texp()`, `temax()`, `tfpoly()`, `tspline()` or `tuser()`
#'
#' @param positive.scale A boolean object that indicates whether all continuous
#'   mean responses (y) are positive and therefore whether the baseline response
#'   should be given a prior that constrains it to be positive (e.g. for scales that cannot be <0).
#' @param intercept A boolean object that indicates whether an intercept is to
#'   be included in the model. Can be used to imply whether mean responses in
#'   data are change from baseline (`FALSE`) or not (setting it to `FALSE`
#'   removes the intercept, `alpha`, from the model).
#' @param link Can take either `"identity"` (the default),
#'   `"log"` (for modelling Ratios of Means \insertCite{friedrich2011}{MBNMAtime}) or
#'   `"smd"` (for modelling Standardised Mean Differences - although this also corresponds to an identity link function).
#'
#' @param rho The correlation coefficient when modelling correlation between time points. The default is a string representing a
#'   prior distribution in JAGS, indicating that it be estimated from the data (e.g. `rho="dunif(0,1)"`). `rho` also be assigned a
#'   numeric value (e.g. `rho=0.7`), which fixes `rho` in the model to this value (e.g. for use in a deterministic sensitivity analysis).
#'   If set to `rho=0` (the default) then this implies modelling no correlation between time points.
#' @param covar A character specifying the covariance structure to use for modelling correlation between time-points. This can
#'   be done by specifying one of the following:
#'   * `"varadj"` - a univariate likelihood with a variance adjustment to assume a constant correlation between subsequent
#'   time points \insertCite{jansen2015}{MBNMAtime}. This is the default.
#'   * `"CS"` - a multivariate normal likelihood with a
#'     \href{https://online.stat.psu.edu/stat502/lesson/10/10.3}{compound symmetry} structure
#'   * `"AR1"` - a multivariate normal likelihood with an
#'     \href{https://online.stat.psu.edu/stat502/lesson/10/10.3}{autoregressive AR1} structure
#'
#' @param var.scale A numeric vector indicating the relative scale of variances between
#' correlated time-course parameters when relative effects are modelled on more than
#' one time-course parameter(see Details). Each element of
#' the vector refers to the relative scale of each of the time-course parameters that is
#' modelled using relative effects.
#'
#' @param class.effect A list of named strings that determines which time-course
#'   parameters to model with a class effect and what that effect should be
#'   (`"common"` or `"random"`). For example: `list(emax="common", et50="random")`.
#' @param UME Can take either `TRUE` or `FALSE` (for an unrelated mean effects
#'   model on all or no time-course parameters respectively) or can be a vector
#'   of parameter name strings to model as UME. For example: `c("beta.1", "beta.2")`.
#'
#' @param pd Can take either:
#'   * `pv` only pV will be reported (as automatically outputted by R2jags).
#'   * `plugin` calculates pD by the plug-in
#'   method \insertCite{spiegelhalter2002}{MBNMAtime}. It is faster, but may output negative
#'   non-sensical values, due to skewed deviances that can arise with non-linear models.
#'   * `pd.kl` (the default) calculates pD by the Kullback–Leibler divergence \insertCite{plummer2008}{MBNMAtime}. This
#'   will require running the model for additional iterations but
#'   will always produce a sensical result.
#'   * `popt` calculates pD using an optimism adjustment which allows for calculation
#'   of the penalized expected deviance \insertCite{plummer2008}{MBNMAtime}
#' @param parallel A boolean value that indicates whether JAGS should be run in
#'   parallel (`TRUE`) or not (`FALSE`). If `TRUE` then the number of cores to
#'   use is automatically calculated. Functions that involve updating the model (e.g. `devplot()`, `fitplot()`)
#'   cannot be used with models implemented in parallel.
#'
#' @param n.iter number of total iterations per chain (including burn in; default: 20000)
#' @param n.thin thinning rate. Must be a positive integer. Set `n.thin > 1`` to save memory
#' and computation time if `n.iter` is large. Default is
#' `max(1, floor(n.chains * (n.iter-n.burnin) / 1000))`` which will only thin if there are at least 2000
#' simulations.
#' @param n.chains number of Markov chains (default: 3)
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the
#' beginning. Default is `n.iter/2``, that is, discarding the first half of the
#' simulations. If `n.burnin` is 0, jags() will run 100 iterations for adaption.
#'
#' @param model.file A JAGS model written as a character object that can be used
#'   to overwrite the JAGS model that is automatically written based on the
#'   specified options. Useful when amending priors using replace.prior()
#' @param ... Arguments to be sent to R2jags.
#'
#' @inheritParams replace.prior
#'
#' @return An object of S3 class `c("mbnma", "rjags")`` containing parameter
#'   results from the model. Can be summarized by `print()` and can check
#'   traceplots using `R2jags::traceplot()` or various functions from the package `mcmcplots`.#'
#'
#'   If there are errors in the JAGS model code then the object will be a list
#'   consisting of two elements - an error message from JAGS that can help with
#'   debugging and `model.arg`, a list of arguments provided to `mb.run()`
#'   which includes `jagscode`, the JAGS code for the model that can help
#'   users identify the source of the error.
#'
#' @section Time-course parameters:
#'   Nodes that are automatically monitored (if present in the model) have the
#'   same name as in the time-course function for named time-course parameters (e.g. `emax`).
#'   However, for named only as `beta.1`, `beta.2`, `beta.3` or `beta.4` parameters
#'   may have an alternative interpretation.
#'
#'   Details of the interpretation and model specification of different parameters can be shown by using the
#'   `summary()` method on an `"mbnma"` object generated by `mb.run()`.
#'
#'   *Parameters modelled using relative effects*
#'   * If pooling is relative (e.g. `pool.1="rel"`) for a given parameter then the named parameter (e.g. `emax`) or a
#'   numbered `d` parameter (e.g. `d.1`) corresponds to the pooled relative effect for a given
#'   treatment compared to the network reference treatment for this time-course parameter.
#'   * `sd.` followed by a named (e.g. `emax`, `beta.1`) is the between-study SD (heterogeneity)
#'   for relative effects, reported if pooling for a time-course parameter is relative (e.g. `pool.1="rel"`) *and* the
#'   method for synthesis is random (e.g. `method.1="random`).
#'   * If class effects are modelled, parameters for classes are represented by the upper case name of the time-course
#'   parameter they correspond to. For example if `class.effect=list(emax="random")`, relative class effects will be
#'   represented by `EMAX`. The SD of the class effect (e.g. `sd.EMAX`, `sd.BETA.1`) is the SD of treatments within a class for the
#'   time-course parameter they correspond to.
#'
#'   *Parameters modelled using absolute effects*
#'   * If pooling is absolute (e.g. `pool.1="abs"`) for a given parameter then the named parameter (e.g. `emax`) or a
#'   numbered `beta` parameter (e.g. `beta.1`) corresponds to the estimated absolute effect for this time-course parameter.
#'   * For an absolute time-course parameter if the corresponding method is common (e.g. `method.1="common"`) the parameter
#'   corresponds to a single common parameter estimated across all studies and treatments. If the corresponding method is
#'   random (e.g. `method.1="random"`) then parameter is a mean effect around which the study-level absolute effects vary
#'   with SD corresponding to `sd.` followed by the named parameter (e.g. `sd.emax`, `sd.beta.1`) .
#'
#'   *Other model parameters*
#'   * `rho` The correlation coefficient for correlation between time-points. Its
#'   interpretation will differ depending on the covariance structure specified in `covar`
#'   * `totresdev` The residual deviance of the model
#'   * `deviance` The deviance of the model
#'
#'
#'
#' @section Time-course function:
#'   Several general time-course functions with up to 4 time-course parameters are provided, but a
#'   user-defined time-course relationship can instead be used. Details can be found in the respective
#'   help files for each function.
#'
#'   Available time-course functions are:
#'   * Log-linear: `tloglin()`
#'   * Polynomial: `tpoly()`
#'   * Exponential: `texp()`
#'   * Emax: `temax()`
#'   * Fractional polynomial: `tfpoly()`
#'   * Splines (various spline types can be used): `tspline()`
#'   * User-defined: `tuser()`
#'
#'
#' @section Correlation between observations:
#'   When modelling correlation between observations using `rho`, values for `rho` must imply a
#'   positive semidefinite covariance matrix.
#'
#' @importFrom Rdpack reprompt
#' @importFrom magrittr "%>%"
#'
#' @references
#'   \insertAllCited
#'
#'
#'
#' @examples
#' \donttest{
#' # Create mb.network object
#' network <- mb.network(osteopain)
#'
#' # Fit a linear time-course MBNMA with:
#' # random relative treatment effects on the slope
#' mb.run(network, fun=tpoly(degree=1, pool.1="rel", method.1="random"))
#'
#' # Fit an emax time-course MBNMA with:
#' # fixed relative treatment effects on emax
#' # a common parameter estimated independently of treatment
#' # a common Hill parameter estimated independently of treatment
#' # a prior for the Hill parameter (normal with mean 0 and precision 0.1)
#' # data reported as change from baseline
#' result <- mb.run(network, fun=temax(pool.emax="rel", method.emax="common",
#'                                     pool.et50="abs", method.et50="common",
#'                                     pool.hill="abs", method.hill="common"),
#'                  priors=list(hill="dnorm(0, 0.1)"),
#'                  intercept=TRUE)
#'
#'
#' # Fit a log-linear MBNMA with:
#' # random relative treatment effects on the rate
#' # an autoregressive AR1 covariance structure
#' # modelled as standardised mean differences
#' result <- mb.run(network, fun=tloglin(pool.rate="rel", method.rate="random"),
#'                  covar="AR1", link="smd")
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(result, c("rate", "sd.rate", "deviance"))
#'
#' # Traceplots
#' mcmcplots::traplot(result)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(result, "rate")
#'
#'
#' ########## Output ###########
#'
#' # Print R2jags output and summary
#' print(result)
#' summary(result)
#'
#' # Plot forest plot of results
#' plot(result)
#'
#'
#' ###### Additional model arguments ######
#'
#' # Use gout dataset
#' goutnet <- mb.network(goutSUA_CFBcomb)
#'
#' # Define a user-defined time-course relationship for use in mb.run
#' timecourse <- ~ exp(beta.1 * time) + (time^beta.2)
#'
#' # Run model with:
#' # user-defined time-course function
#' # random relative effects on beta.1
#' # default common effects on beta.2
#' # default relative pooling on beta.1 and beta.2
#' # common class effect on beta.2
#' mb.run(goutnet, fun=tuser(fun=timecourse, method.1="random"),
#'        class.effect=list(beta.1="common"))
#'
#' # Fit a log-linear MBNMA
#' # with variance adjustment for correlation between time-points
#' result <- mb.run(network, fun=tloglin(),
#'                  rho="dunif(0,1)", covar="varadj")
#' }
#' @export
mb.run <- function(network, fun=tpoly(degree = 1), positive.scale=FALSE, intercept=TRUE,
                      link="identity",
                      parameters.to.save=NULL,
                      rho=0, covar="varadj",
                      var.scale=NULL,
                      class.effect=list(), UME=FALSE,
                      pd="pd.kl", parallel=FALSE,
                      priors=NULL,
                      n.iter=20000, n.chains=3,
                      n.burnin=floor(n.iter/2), n.thin=max(1, floor((n.iter - n.burnin) / 1000)),
                      model.file=NULL, ...
) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(fun, classes = "timefun", add=argcheck)
  checkmate::assertClass(network, "mb.network", add=argcheck)
  checkmate::assertCharacter(model.file, len=1, any.missing=FALSE, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(pd, choices=c("pv", "pd.kl", "plugin", "popt"), null.ok=FALSE, add=argcheck)
  checkmate::assertLogical(parallel, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertList(priors, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Reduce n.burnin by 1 to avoid JAGS error if n.burnin=n.iter
  if (n.iter==n.burnin) {
    n.burnin <- n.burnin - 1
  }


  if (is.null(model.file)) {
    model <- mb.write(fun=fun, link=link,
                      positive.scale=positive.scale, intercept=intercept,
                      rho=rho, covar=covar,
                      class.effect=class.effect, UME=UME,
                      var.scale=var.scale
    )

    if (!is.null(priors)) {
      model <- replace.prior(priors=priors, model=model)
    }

  } else {
    warning("All parameter specifications (time-course, rho, class effects, UME, priors, etc.) are being overwritten by `model.file`")
    model <- model.file
  }

  assigned.parameters.to.save <- parameters.to.save
  if (is.null(parameters.to.save)) {
    parameters.to.save <-
      gen.parameters.to.save(fun=fun, model=model)
  }

  # Add nodes to monitor to calculate plugin pd
  if (pd=="plugin") {
    if (covar!="varadj") {
      stop("pD cannot be calculated via the plugin method if modelling a multivariate normal likelihood - covar!='varadj'")
    }

    pluginvars <- c("theta", "resdev")
    for (param in seq_along(pluginvars)) {
      if (!(pluginvars[param] %in% parameters.to.save)) {
        parameters.to.save <- append(parameters.to.save, pluginvars[param])
      }
    }
    message("The following parameters have been monitored to allow pD plugin calculation: ",
            paste(pluginvars, collapse=", "))
  }

  if (length(class.effect)>0) {
    class <- TRUE
  } else {class <- FALSE}


  #### Run jags model ####

  data.ab <- network[["data.ab"]]
  result.jags <- mb.jags(data.ab, model, fun=fun, link=link,
                       class=class, rho=rho, covar=covar,
                       parameters.to.save=parameters.to.save,
                       n.iter=n.iter, n.chains=n.chains,
                       n.burnin=n.burnin, n.thin=n.thin,
                       ...)
  result <- result.jags[["jagsoutput"]]
  jagsdata <- result.jags[["jagsdata"]]

  if (!("error" %in% names(result))) {
    if (pd == "pd.kl" | pd == "popt") {
      if (pd=="pd.kl") {
        temp <- rjags::dic.samples(result$model, n.iter=n.iter/10, type="pD")
      } else if (pd=="popt") {
        temp <- rjags::dic.samples(result$model, n.iter=n.iter/10, type="popt")
      }
      result$BUGSoutput$pD <- sum(temp$penalty)

    } else if (pd == "plugin") {
      # plugin method
      warning("Plugin method only works for normal likelihood")
      result$BUGSoutput$pD <- pDcalc(obs1=jagsdata[["y"]], obs2=jagsdata[["se"]], fups=jagsdata[["fups"]], narm=jagsdata[["narm"]], NS=jagsdata[["NS"]],
                                     theta.result=result$BUGSoutput$mean$theta, resdev.result=result$BUGSoutput$mean$resdev,
                                     type="time")
    }

    # Recalculate DIC so it is adjusted for choice of pD
    result$BUGSoutput$DIC <- result$BUGSoutput$pD + result$BUGSoutput$median$deviance
  }

  # Add variables for other key model characteristics (for predict and plot functions)
  model.arg <- list("parameters.to.save"=assigned.parameters.to.save,
                    "fun"=fun,
                    "jagscode"=model, "jagsdata"=jagsdata,
                    "link"=link,
                    "positive.scale"=positive.scale, "intercept"=intercept,
                    "rho"=rho, "covar"=covar,
                    "class.effect"=class.effect, "UME"=UME,
                    "var.scale"=var.scale,
                    "parallel"=parallel, "pd"=pd,
                    "priors"=get.prior(model))
  result[["model.arg"]] <- model.arg
  result[["network"]] <- network
  # result[["treatments"]] <- network[["treatments"]]
  # if ("classes" %in% names(network)) {
  #   result[["classes"]] <- network[["classes"]]
  # }
  result[["type"]] <- "time"

  if (!("error" %in% names(result))) {
    class(result) <- c("mbnma", class(result))
  }

  return(result)

}


mb.jags <- function(data.ab, model, fun=NULL, link=NULL,
                       class=FALSE, rho=NULL, covar=NULL,
                       parameters.to.save=parameters.to.save,
                       likelihood=NULL,
                       warn.rhat=FALSE, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertCharacter(model, any.missing=FALSE, add=argcheck)
  checkmate::assertLogical(class, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertCharacter(parameters.to.save, any.missing=FALSE, unique=TRUE,
                  null.ok=TRUE, add=argcheck)
  checkmate::assertClass(fun, "timefun", add=argcheck)
  checkmate::reportAssertions(argcheck)


  if (is.null(likelihood)) {
    # For MBNMAtime
    jagsdata <- getjagsdata(data.ab, class=class, rho=rho, covstruct=covar, fun=fun, link=link) # get data into jags correct format (list("fups", "NT", "NS", "narm", "y", "se", "treat", "time"))
  } else if (is.null(rho) & is.null(covar)) {
    # For MBNMAdose
    # jagsdata <- getjagsdata(data.ab, class=class,
    #                         likelihood=likelihood, link=link) # get data into jags correct format
  }


  # Add variable for maxtime to jagsdata if required
  if (any(grepl("maxtime", model))) {
    jagsdata[["maxtime"]] <- max(data.ab$time)
  }

  # Remove studyID from jagsdata (not used in model)
  tempjags <- jagsdata
  tempjags[["studyID"]] <- NULL

  # Drop time from tempjags in spline models
  if (fun$name %in% c("rcs", "ns", "bs", "ls") & !"AR1" %in% covar) {
    tempjags[["time"]] <- NULL
  }

  # Put data from jagsdata into separate R objects
  for (i in seq_along(tempjags)) {
    ##first extract the object value
    temp <- tempjags[[i]]
    ##now create a new variable with the original name of the list item
    eval(parse(text=paste(names(tempjags)[[i]],"<- temp")))
  }

  # Take names of variables in tempjags for use in rjags
  jagsvars <- list()
  for (i in seq_along(names(tempjags))) {
    jagsvars[[i]] <- names(tempjags)[i]
  }

  # Create a temporary model file
  tmpf=tempfile()
  tmps=file(tmpf,"w")
  cat(paste(model, collapse="\n"),file=tmps)
  close(tmps)

  out <- tryCatch({
    result <- R2jags::jags(data=jagsvars, model.file=tmpf,
                           parameters.to.save=parameters.to.save,
                           ...
    )
  },
  error=function(cond) {
    message(cond)
    return(list("error"=cond))
  }
  )

  # Gives warning if any rhat values > 1.02
  if (warn.rhat==TRUE) {
    if (!("error" %in% names(out))) {
      rhat.warning(out)
    }
  }

  return(list("jagsoutput"=out, "jagsdata"=jagsdata))
}




#' Automatically generate parameters to save for a dose-response MBNMA model
#'
#' @inheritParams mb.run
#' @param model A JAGS model written as a character object
gen.parameters.to.save <- function(fun, model) {
  # model.params is a vector (numeric/character) of the names of the dose-response parameters in the model
  #e.g. c(1, 2, 3) or c("emax", "et50")
  # model is a JAGS model written as a character object

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(fun, classes = "timefun", add=argcheck)
  checkmate::assertCharacter(model, min.len = 10)
  checkmate::reportAssertions(argcheck)


  # Set some automatic parameters based on the model code
  parameters.to.save <- vector()
  for (i in seq_along(fun$params)) {

    # For unnamed parameters
    if (any(grepl(paste0("^d\\.", i), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("d.", i))
    }
    if (any(grepl(paste0("^D\\.", i), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("D.", i))
    }
    if (any(grepl(paste0("^sd\\.", i), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.", i))
    }
    if (any(grepl(paste0("^sd\\.D\\.", i), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.D.", i))
    }

    # For named parameters
    if (any(grepl(fun$params[i], model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, fun$params[i])
    }
    if (any(grepl(toupper(fun$params[i]), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, toupper(fun$params[i]))
    }
    if (any(grepl(paste0("^sd\\.", fun$params[i]), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.", fun$params[i]))
    }
    if (any(grepl(paste0("^sd\\.", toupper(fun$params[i])), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.", toupper(fun$params[i])))
    }

    # Remove if both d and beta are in for any parameter
    if (paste0("d.",i) %in% parameters.to.save & paste0("beta.",i) %in% parameters.to.save) {
      parameters.to.save <- parameters.to.save[!parameters.to.save %in% paste0("beta.",i)]
    }
  }

  # For MBNMAtime
  if (any(grepl("rho", model))==TRUE) {
    parameters.to.save <- append(parameters.to.save, "rho")
  }
  if (any(grepl("totresdev", model))==TRUE) {
    parameters.to.save <- append(parameters.to.save, c("totresdev"))
  }

  return(unique(parameters.to.save))
}







#' Calculate plugin pD from a JAGS model with univariate likelihood for studies
#' with repeated measurements
#'
#' Uses results from MBNMA JAGS models to calculate pD via the
#' plugin method \insertCite{spiegelhalter2002}{MBNMAtime}. Can only be used for models with known
#' standard errors or covariance matrices (typically univariate).
#'
#' @param obs1 A matrix (study x arm) or array (study x arm x time point) containing
#'   observed data for `y` (normal likelihood) or `r` (binomial or Poisson likelihood)
#'   in each arm of each study. This will be the same array
#'   used as data for the JAGS model.
#' @param obs2 A matrix (study x arm) or array (study x arm x time point) containing
#'   observed data for `se` (normal likelihood), `N` (binomial likelihood) or `E` (Poisson likelihood)
#'   in each arm of each study. This will be the same array
#'   used as data for the JAGS model.
#' @param fups A numeric vector of length equal to the number of studies,
#'   containing the number of follow-up mean responses reported in each study. Required for
#'   time-course MBNMA models (if `type="time"`)
#' @param narm A numeric vector of length equal to the number of studies,
#'   containing the number of arms in each study.
#' @param NS A single number equal to the number of studies in the dataset.
#' @param theta.result A matrix (study x arm) or array (study x arm x time point)
#'   containing the posterior mean predicted means/probabilities/rate in each arm of each
#'   study. This will be estimated by the JAGS model.
#' @param resdev.result A matrix (study x arm) or array (study x arm x time point)
#'   containing the posterior mean residual deviance contributions in each arm of each
#'   study. This will be estimated by the JAGS model.
#'
#' @param likelihood A character object of any of the following likelihoods:
#' * `univariate`
#' * `binomial` (does not work with time-course MBNMA models)
#' * `multivar.normal` (does not work with time-course MBNMA models)
#' @param type The type of MBNMA model fitted. Can be either `"time"` or `"dose"`
#'
#' @details Method for calculating pD via the plugin method proposed by
#'   \insertCite{spiegelhalter2002}{MBNMAtime}. Standard errors / covariance matrices must be assumed
#'   to be known. To obtain values for theta.result and resdev.result these
#'   parameters must be monitored when running the JAGS model.
#'
#'   For non-linear time-course MBNMA models residual deviance contributions may be skewed, which
#'   can lead to non-sensical results when calculating pD via the plugin method.
#'   Alternative approaches are to use pV (`pv`) as an approximation \insertCite{plummer2008}{MBNMAtime} or
#'   pD calculated by Kullback–Leibler divergence (`pd.kl`) or using an optimism adjustment (`popt`) \insertCite{plummer2008}{MBNMAtime}.
#'
#' @references TO ADD pV REF
#' @inherit mb.run references
#'
#' @examples
#' \donttest{
#' # Using the alogliptin dataset
#' network <- mb.network(alog_pcfb)
#'
#' # Run Emax model saving predicted means and residual deviance contributions
#' emax <- mb.emax(network, parameters.to.save=c("theta", "resdev"), intercept=FALSE)
#'
#' # Get matrices of observed data
#' jagsdat <- getjagsdata(network$data.ab)
#'
#' # Plugin estimation of pD is problematic with non-linear models as it often leads to
#' #negative values, hence use of pV, pd.kl and popt as other measures for the effective
#' #number of parameters
#' pDcalc(obs1=jagsdat$y, obs2=jagsdat$se,
#'   fups=jagsdat$fups, narm=jagsdat$narm, NS=jagsdat$NS,
#'   theta.result = emax$BUGSoutput$mean$theta,
#'   resdev.result = emax$BUGSoutput$mean$resdev
#'   )
#' }
#' @export
pDcalc <- function(obs1, obs2, fups=NULL, narm, NS, theta.result, resdev.result,
                   likelihood="normal", type="time") {
  # For univariate models only!!

  # likelihood could in theory be c("normal", "multivar.normal", "binomial")
  # theta.result = model$BUGSoutput$mean$theta
  # resdev.result = model$BUGSoutput$mean$totresdev

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertArray(obs1, add=argcheck)
  checkmate::assertArray(obs2, add=argcheck)
  checkmate::assertArray(theta.result, add=argcheck)
  checkmate::assertArray(resdev.result, add=argcheck)
  checkmate::assertNumeric(fups, null.ok=TRUE, add=argcheck)
  checkmate::assertNumeric(narm, add=argcheck)
  checkmate::assertNumeric(NS, add=argcheck)
  checkmate::assertChoice(likelihood, choices=c("normal", "binomial", "poisson"), null.ok=FALSE, add=argcheck)
  checkmate::assertChoice(type, choices=c("dose", "time"), null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (type=="time") {
    if (is.null(fups)) {
      stop("`fups` cannot be NA in pDcalc for time-course MBNMA")
    }
    dev.post <- array(dim=c(NS,max(narm),max(fups)))
    pD <- array(dim=c(NS,max(narm),max(fups)))
  } else if (type=="dose") {
    dev.post <- matrix(nrow=NS, ncol=max(narm))
    pD <- matrix(nrow=NS, ncol=max(narm))
  }

  for (i in 1:NS) {
    for (k in 1:narm[i]) {

      if (type=="time") {
        for (m in 1:fups[i]) {
          # Need to use formula for residual deviance as plugin
          if (likelihood=="normal") {
            dev.post[i,k,m] <- ((obs1[i,k,m] - theta.result[i,k,m])/obs2[i,k,m])^2
            pD[i,k,m] <- resdev.result[i,k,m] - dev.post[i,k,m]
          } else {
            stop("pD cannot be calculated via `plugin` method for time-course MBNMA models without data following a normal likelihood")
          }
        }
      } else if (type=="dose") {
        if (likelihood=="normal") {
          dev.post[i,k] <- ((obs1[i,k] - theta.result[i,k])/obs2[i,k])^2

        } else if (likelihood=="binomial") {
          rhat[i,k] <- theta.result[i,k] * obs2[i,k]
          dev.post[i,k] <- 2*(obs1[i,k] * (log(obs1[i,k])-log(rhat[i,k]))  +
                                (obs2[i,k]-obs1[i,k]) * (log(obs2[i,k]-obs1[i,k]) -
                                                           log(obs2[i,k]-rhat[i,k])))
        } else if (likelihood=="poisson") {
          rhat[i,k] <- theta.result[i,k] * obs2[i,k]
          dev.post[i,k] <- 2*((rhat[i,k]-obs1[i,k]) + (obs1[i,k] * (log(obs1[i,k]/rhat[i,k]))))
        }

        pD[i,k] <- resdev.result[i,k] - dev.post[i,k]

      }

    }
  }

  pD <- sum(pD, na.rm=TRUE)

  return(pD)
}






#' Prepares beta time-course parameters for mb.write()
#'
#' Checks that beta time-course parameters have been specified correctly and converts them
#' to the correct format for `mb.write()` and other subsequent functions.
#'
#' @param beta.1 A two-element list whose elements have the following names:
#' * `pool` Can take either `"rel"`, `"arm"` or `"const"`
#' * `method` Can take either `"common"`, `"random"`, or be assigned a numeric value
compound.beta <- function(beta.1) {
  # Checks
  checkmate::assertList(beta.1, len=2, types=c("character", "numeric"), null.ok = FALSE,
                        unique=TRUE)

  # If the list is not named then name first and second elements pool and method respectively
  if (is.null(names(beta.1))) {
    names(beta.1) <- c("pool", "method")
  }

  if (!identical(sort(names(beta.1)), c("method", "pool"))) {
    stop("Time-course parameters must be each be a list containing named elments `pool` and `method`")
  }

  if (!(beta.1$pool %in% c("rel", "arm", "const"))) {
    stop("`pool` must be specified as either `rel`, `arm`, or `const`")
  }

  if (is.numeric(beta.1$method) & beta.1$pool!="const") {
    stop("A numeric value can only be assigned to a time-course parameter if `pool`=`const`")
  }
  if (is.character(beta.1$method) & !(beta.1$method %in% c("common", "random"))) {
    stop("`method` can either take the values `common` or `random`, or can be assigned a numeric value if `pool`=`const`")
  }

  if (beta.1$pool=="const" & is.numeric(beta.1$method)) {
    beta.out <- beta.1$method
  } else {
    beta.out <- paste(beta.1$pool, beta.1$method, sep=".")
  }
  return(beta.out)
}





check.beta.arg <- function(beta.1) {
  # Checks
  checkmate::assertList(beta.1, len=2, types=c("character", "numeric"),
                        unique=TRUE)

  if (!is.null(beta.1)) {
    if (is.null(names(beta.1))) {
      newbeta <- list()
      for (i in seq_along(beta.1)) {
        if (is.character(beta.1[[i]])) {
          if (beta.1[[i]] %in% c("rel", "arm", "const")) {
            newbeta$pool <- beta.1[[i]]
          } else if (beta.1[[i]] %in% c("common", "random")) {
            newbeta$method <- beta.1[[i]]
          }
        } else if (is.numeric(beta.1[[i]])) {
          newbeta$method <- beta.1[[i]]
        } else {
          stop("Time-course parameter arguments incorrectly specified")
        }
      }
      beta.1 <- newbeta
    }

    if (!all(names(beta.1) %in% c("pool", "method"))) {
      stop("Time-course parameter arguments incorrectly named")
    }
  }

  return(beta.1)
}






#' Update MBNMA to obtain deviance contributions or fitted values
#'
#' @inheritParams devplot
#' @inheritParams R2jags::jags
#' @param param A character object that represents the parameter within the model to monitor when updating. Can
#' currently only be used for monitoring fitted values and deviance contributions and so can take
#' either `"dev"` (for deviance contributions), `"resdev"` (for residual deviance contributions)
#' or `"theta"` (for fitted values).
#'
#' @return A data frame containing posterior means for the specified `param` at each observation, arm and study.
#'
#' @examples
#' \donttest{
#' # Using the alogliptin dataset
#' network <- mb.network(alog_pcfb)
#'
#' # Run Emax model
#' emax <- mb.emax(network, intercept=FALSE)
#'
#' # Update model for 500 iterations to monitor fitted values
#' mb.update(emax, param="theta", n.iter=500)
#'
#' # Update model for 500 iterations to monitor residual deviance contributions
#' mb.update(emax, param="resdev", n.iter=500)
#'
#' # Update model for 500 iterations to monitor deviance contributions
#' mb.update(emax, param="dev", n.iter=500)
#' }
#' @export
mb.update <- function(mbnma, param="theta",
                      n.iter=mbnma$BUGSoutput$n.iter, n.thin=mbnma$BUGSoutput$n.thin) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertChoice(param, choices = c("dev", "resdev", "theta"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  # if (!(grepl(paste0("\n", param), mbnma$model.arg$jagscode)==TRUE |
  #       grepl(paste0(" ", param), mbnma$model.arg$jagscode)==TRUE)) {
  if (!(any(grepl(paste0("^", param), mbnma$model.arg$jagscode)==TRUE) |
        any(grepl(paste0(" ", param), mbnma$model.arg$jagscode)==TRUE))) {
    stop(paste0(param, " not in model code"))
  }

  result <- rjags::jags.samples(mbnma$model, variable.names = param,
                                n.iter=n.iter, n.thin=n.thin)

  # Take means of posteriors and convert to data.frame with indices
  update.mat <- apply(result[[param]], c(1,2,3), function(x) mean(x, na.rm=TRUE))
  update.df <- reshape2::melt(update.mat)
  names(update.df) <- c("study", "arm", "fup", "mean")

  # Remove missing values
  update.df <- update.df[stats::complete.cases(update.df),]

  return(update.df)
}
