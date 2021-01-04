# Functions for running MBNMA models
# Author: Hugo Pedder
# Date created: 2018-09-10


#' Run MBNMA time-course models
#'
#' Fits a Bayesian time-course model for model-based network meta-analysis
#' (MBNMA) that can account for repeated measures over time within studies by
#' applying a desired time-course function. Follows the methods of \insertCite{pedder2019;textual}{MBNMAtime}.
#'
#' @param network An object of class `mb.network`.
#' @param parameters.to.save A character vector containing names of parameters
#'   to monitor in JAGS
#' @param fun is a character specifying a functional form to be assigned to the
#'   time-course. Options are given in `details`.
#' @param user.fun A formula specifying any relationship including `time` and
#'   one/several of: `beta.1`, `beta.2`, `beta.3`, `beta.4`.
#' @param model.file A JAGS model written as a character object that can be used
#'   to overwrite the JAGS model that is automatically written based on the
#'   specified options. Useful when amending priors using replace.prior()
#'
#' @param alpha Refers to the baseline mean response and is a character object
#'   that can take either:
#'   * `"study"` to constrain baseline to be equal for all
#'   arms within a study (`i` index is added))
#'   * `"arm"` to allow baseline to
#'   vary between arms within a study (`i,k` index is added)).
#' @param beta.1 A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' @param beta.2 A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' @param beta.3 A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' @param beta.4 A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' @param knots The number/location of knots if a restricted cubic spline time-course function is fitted (`fun="rcs"`).
#' If a single number is given it indicates the number of knots (they will
#'   be equally spaced across the range of time points). If a numeric vector is given it indicates the location of the knots.
#'   Minimum number of knots is 3.
#' @param positive.scale A boolean object that indicates whether all continuous
#'   mean responses (y) are positive and therefore whether the baseline response
#'   should be given a prior that constrains it to be positive.
#' @param intercept A boolean object that indicates whether an intercept is to
#'   be included in the model. Can be used to imply whether mean responses in
#'   data are change from baseline (`FALSE`) or not (setting it to `FALSE`
#'   removes the intercept, `alpha`, from the model).
#' @param link Can take either `"identity"`, `"log"` (for modelling Ratios of Means \insertCite{friedrichROM}{MBNMAtime}) or
#'   `"smd"` (for modelling Standardised Mean Differences).
#'
#' @param rho The correlation coefficient when modelling correlation between time points. If left
#'   as `NULL` (the default) then this implies modelling no correlation between time points.
#'   Can either be assigned the string `"estimate"` to indicate that rho should be estimated
#'   from the data, or assigned a numeric value, which fixes `rho` in the model to the assigned
#'   value, either for when `rho` is calculated externally or for use in deterministic sensitivity
#'   analyses.
#' @param covar A character specifying the covariance structure to use for the
#'   multivariate normal likelihood. Can currently take either `"CS"` (compound
#'   symmetry) or `"AR1"` (autoregressive AR1).
#'
#' @param var.scale A numeric vector indicating the relative scale of variances between
#' correlated time-course parameters when relative effects are modelled on more than
#' one time-course parameter(see Details LINK). Each element of
#' the vector refers to the relative scale of each of the time-course parameters that is
#' modelled using relative effects.
#'
#' @param class.effect A list of named strings that determines which time-course
#'   parameters to model with a class effect and what that effect should be
#'   (`"common"` or `"random"`). For example: `list("beta.2"="common", "beta.3"="random")`.
#' @param UME Can take either `TRUE` or `FALSE` (for an unrelated mean effects
#'   model on all or no time-course parameters respectively) or can be a vector
#'   of parameter name strings to model as UME. For example: `c("beta.1", "beta.2")`.
#'
#' @param pd Can take either:
#'   * `pv` only pV will be reported (as automatically outputted by R2jags).
#'   * `plugin` calculates pD by the plug-in
#'   method \insertCite{spiegelhalter2002}{MBNMAtime}. It is faster, but may output negative
#'   non-sensical values, due to skewed deviances that can arise with non-linear models.
#'   * `pd.kl` calculates pD by the Kullback–Leibler divergence \insertCite{plummer2008}{MBNMAtime}. This
#'   will require running the model for additional iterations but
#'   will always produce a positive result.
#'   * `popt` calculates pD using an optimism adjustment which allows for calculation
#'   of the penalized expected deviance \insertCite{plummer2008}{MBNMAtime}
#' @param parallel A boolean value that indicates whether JAGS should be run in
#'   parallel (`TRUE`) or not (`FALSE`). If `TRUE` then the number of cores to
#'   use is automatically calculated.
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
#' @param arg.params Contains a list of arguments sent to `mb.run()` by time-course
#' specific wrapper functions
#' @param ... Arguments to be sent to R2jags.
#'
#' @inheritParams replace.prior
#'
#' @return An object of S3 class `c("mbnma", "rjags")`` containing parameter
#'   results from the model. Can be summarized by `print()` and can check
#'   traceplots using `R2jags::traceplot()` or various functions from the package `mcmcplots`.
#'
#'   Nodes that are automatically monitored (if present in the model) have the
#'   following interpretation. They will have an additional suffix that relates
#'   to the name/number of the time-course parameter to which they correspond
#'   (e.g. `d.et50` or `d.1`):
#'   * `d` The pooled relative effect for a given
#'   treatment compared to the network reference treatment for a particular
#'   time-course parameter, reported if `pool="rel"`
#'   * `sd.d` The between-study SD (heterogeneity) for relative effects,
#'   reported if `pool="rel"` and `method="random"`
#'   * `D` The class effect for a given class compared to the
#'   network reference class for a particular time-course parameter
#'   * `sd.D` The standard deviation for the pooled relative effects of treatments within a
#'   given class from a model with a random class effect.
#'   * `beta` If `pool="const"` then only a single node will be present in the
#'   output, which corresponds to the absolute value of a particular time-course
#'   parameter across the network, If `pool="arm"`
#'   then for the relevant time-course parameter there will be one node for
#'   each treatment, which represents the absolute value of the time-course
#'   parameter for each treatment
#'   * `sd.beta` Reported if `method="random"` and `pool` is either `"const"` or `"arm"`.
#'   If `pool="const"` this represents the between-study SD for the absolute value of a particular
#'   time-course parameter exchangeable across the network. If `pool="arm"`
#'   this represents the between-study SD for the absolute value of a
#'   particular time-course parameter exchangeable by treatment
#'   * `rho` The correlation coefficient for correlation between time-points. Its
#'   interpretation will differ depending on the covariance structure used
#'   * `totresdev` The residual deviance of the model
#'   * `deviance` The deviance of the model
#'
#'
#'   If there are errors in the JAGS model code then the object will be a list
#'   consisting of two elements - an error message from JAGS that can help with
#'   debugging and `model.arg`, a list of arguments provided to `mb.run()`
#'   which includes `jagscode`, the JAGS code for the model that can help
#'   users identify the source of the error.
#'
#' @section Time-course parameters:
#' Time-course parameters in the model must be provided as a list with named elements
#' `pool` and `method`.
#'
#' `pool` is used to define the approach used for pooling of a given time-course parameter and
#' can take any of the following values:
#' * `"rel"` indicates that relative effects should be pooled for this time-course parameter.
#' This preserves randomisation
#' within included studies, are likely to vary less between studies (only due to effect modification),
#' and allow for testing of consistency between direct and indirect evidence. Pooling follows the
#' general approach for Network Meta-Analysis proposed by \insertCite{lu2004;textual}{MBNMAtime}.
#' * `"arm"` indicates that study arms should be pooled within each treatment in the network
#' for this time-course parameter.
#' This allows estimation of absolute time-course parameter values, but makes stronger assumptions
#' regarding similarity of studies.
#' * `"const"` indicates that study arms should be pooled across the whole network for this
#' time-course parameter  *independently of assigned treatment*.
#' This implies using a single value across the network for this time-course parameter,
#' and may therefore be making very strong assumptions of similarity.
#'
#' `method` is used to define the model used for meta-analysis for a given time-course parameter
#' and can take any of the following values:
#' * `"common"` implies that all studies estimate the same true effect
#' (akin to a "fixed effects" meta-analysis)
#' * `"random"` implies that all studies estimate a separate true effect, but that each
#' of these true effects vary randomly around a true mean effect. This approach allows
#' for modelling of between-study heterogeneity.
#' * `numeric()` Assigned a numeric value - this can only be used if `pool="const"`. It indicates that
#' this time-course parameter should not be estimated from the data but should be assigned
#' the numeric value determined by the user. This can be useful for fixing specific time-course
#' parameters (e.g. Hill parameters in Emax functions or knot location in piecewise functions).
#'
#' When relative effects are modelled on more than one time-course parameter,
#' correlation between the time-course parameters is automatically
#' estimated using a vague Wishart prior. This prior can be made slightly more informative
#' by specifying the relative scale of variances between the time-course parameters using
#' `var.scale`.
#'
#'
#' @section Time-course function:
#'   Several general time-course functions are provided, but a
#'   user-defined time-course relationship can instead be used.
#'
#'   Built-in time-course functions are:
#'   * `"linear"`: `beta.1` refers to the gradient
#'   * `"quadratic"`: `beta.1` refers to the gradient, `beta.2` refers
#'   to the change in gradient
#'   * `"exponential"`: `beta.1` refers to the rate of gain/decay
#'   * `"emax"` (emax without a Hill parameter): `beta.1` refers to
#'   Emax parameter, `beta.2` refers to ET50 parameter
#'   * `"emax.hill"` (emax with a Hill parameter): `beta.1` refers to Emax parameter, `beta.2` refers
#'   to ET50 parameter, `beta.3` refers to Hill parameter
#'   * `"rcs"` restricted cubic splines with knot number/location defined by `knot`.`beta.1` refers to the
#'   first spline coeffficient, `beta.2` to the second coefficient, etc.
#'   * `"fract.poly.first"` (first-order fractional polynomial): `beta.1` refers to the slope
#'   parameter, `beta.3` refers to the power parameter
#'   * `"fract.poly.second"` (second-order fractional polynomial): `beta.1` refers to the first slope
#'   parameter, `beta.2` refers to the first power parameter, `beta.3` refers to
#'   the first power parameter, `beta.4` refers to the second power parameter
#'   * `"piecelinear"` piecewise linear: `beta.1` refers to the gradient of the
#'   first linear piece, `beta.2` refers to the gradient of the second linear
#'   piece, `beta.3` refers to the knot location (the time at which the two
#'   pieces are joined)
#'   * `"user"` (user-defined function: `user.fun` must be specified in arguments)
#'
#' @section Correlation between observations:
#'   When modelling correlation between observations using `rho`, values for `rho` must imply a
#'   positive semidefinite covariance matrix. If estimating `rho` from the data (by assigning it
#'   `"estimate"`), the default prior distribution (`dunif(-1,1)`) may include values that exclude
#'   a positive semidefinite covariance matrix. This prior can be restricted (e.g. to `dunif(0,1)`)
#'   using the `priors` argument (see \code{\link{get.prior}})
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
#' # Fit a linear time-course MBNMA with random consistency treatment effects on beta.1 (slope)
#' #and equal baselines
#' #in study arms
#' mb.run(network, fun="linear",
#'   alpha="study", beta.1=list(pool="rel", method="random"))
#'
#' # Fit an emax time-course MBNMA with fixed consistency treatment effects on beta.1 (emax),
#' #a common parameter estimated across the network for beta.2 (et50) and a Hill parameter of 0.5
#' #across the network on data reported as change from baseline
#' result <- mb.run(network, fun="emax.hill",
#'   beta.1=list(pool="rel", method="common"),
#'   beta.2=list(pool="const", method="common"),
#'   beta.3=list(pool="const", method=0.5),
#'   intercept=TRUE)
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(result, c("beta.2", "deviance", "totresdev"))
#'
#' # Traceplots
#' mcmcplots::traplot(result)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(result, "d.1")
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
#' # Fit a linear time-course MBNMA that accounts for correlation between time points
#' mb.run(network, fun="linear",
#'   beta.1=list(pool="rel", method="common"),
#'   rho="estimate", covar="CS")
#'
#' # Define a user-defined time-course relationship for use in mb.run
#' time.fun <- ~alpha + exp(beta.1 * time) + (time^beta.2)
#'
#' # Run model using Kullback-Liebler divergence to calculate pD
#' mb.run(network, fun="user", user.fun=time.fun,
#'   beta.1=list(pool="rel", method="random"),
#'   beta.2=list(pool="rel", method="common"),
#'   pd="pd.kl")
#' }
#' @export
mb.run <- function(network, parameters.to.save=NULL,
                      fun="linear", user.fun=NULL,
                      alpha="study", beta.1=list(pool="rel", method="common"),
                      beta.2=NULL, beta.3=NULL, beta.4=NULL,
                      knots=3, link="identity",
                      positive.scale=FALSE, intercept=TRUE, rho=NULL, covar=NULL,
                      var.scale=NULL,
                      class.effect=list(), UME=FALSE,
                      pd="pd.kl", parallel=FALSE,
                      priors=NULL,
                      n.iter=20000, n.chains=3,
                      n.burnin=floor(n.iter/2), n.thin=max(1, floor((n.iter - n.burnin) / 1000)),
                      model.file=NULL,
                      arg.params=NULL, ...
) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mb.network", add=argcheck)
  checkmate::assertCharacter(model.file, len=1, any.missing=FALSE, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(pd, choices=c("pv", "pd.kl", "plugin", "popt"), null.ok=FALSE, add=argcheck)
  checkmate::assertLogical(parallel, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertList(arg.params, unique=TRUE, null.ok=TRUE, add=argcheck)
  checkmate::assertList(priors, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Reduce n.burnin by 1 to avoid JAGS error if n.burnin=n.iter
  if (n.iter==n.burnin) {
    n.burnin <- n.burnin - 1
  }

  # Check betas are specified correctly and prepare format for subsequent functions
  for (i in 1:4) {
    betaname <- paste0("beta.", i)
    if (!is.null(get(betaname))) {
      assign(betaname, check.beta.arg(get(betaname)))
      assign(paste0(betaname, ".str"), compound.beta(get(betaname)))
    } else if (is.null(get(betaname))) {
      assign(paste0(betaname, ".str"), NULL)
    }
  }

  if (!is.null(arg.params)) {
    if (!all((names(arg.params)) %in% c("wrap.params", "run.params"))) {
      stop("arg.params has been incorrectly specified")
    }
    #wrap.params <- c("emax", "et50")
    wrap.params <- arg.params$wrap.params
    #run.params <- c("beta.1", "beta.2")
    run.params <- arg.params$run.params

    fun.params <- list(names(class.effect), UME)

    # Check class effect and UME names match wrapper arguments
    if (!all(names(class.effect) %in% wrap.params)) {
      stop(paste0("The following list element names in `class.effect` do not match time-course parameter names for the given function:\n",
                  paste(names(class.effect)[!(names(class.effect) %in% wrap.params)], collapse=", ")))
    }
    if (is.character(UME) & !all(UME %in% wrap.params)) {
      stop(paste0("The following parameters specified in `UME` do not match time-course parameter names for the given function:\n",
                  paste(UME[!(UME %in% wrap.params)], collapse=", ")))
    }

    for (i in seq_along(fun.params)) {
      for (k in seq_along(fun.params[[i]])) {
        for (m in seq_along(wrap.params)) {
          if (wrap.params[m] %in% fun.params[[i]][k]) {
            fun.params[[i]][k] <- run.params[m]
          }
        }
      }
      # for (k in seq_along(wrap.params)) {
      #   print(wrap.params[k])
      #   if (wrap.params[k] %in% fun.params[[i]]) {
      #     fun.params[[i]] <- run.params[k]
      #   }
      # }
    }

    names(class.effect) <- fun.params[[1]]
    UME <- fun.params[[2]]
  }

  if (is.null(model.file)) {
    model <- mb.write(fun=fun, user.fun=user.fun,
                      alpha=alpha,
                      beta.1=beta.1.str, beta.2=beta.2.str,
                      beta.3=beta.3.str, beta.4=beta.4.str,
                      knots=knots, link=link,
                      positive.scale=positive.scale, intercept=intercept,
                      rho=rho, covar=covar,
                      class.effect=class.effect, UME=UME,
                      var.scale=var.scale
    )

    # Change beta.1 and beta.2 to emax and et50, etc. if necessary
    if (!is.null(arg.params)) {
      code.params <- c("d", "delta", "sd", "beta", "D", "sd.D", "BETA", "sd.BETA")
      for (i in seq_along(wrap.params)) {
        for (k in seq_along(code.params)) {
          model <- gsub(paste(code.params[k], strsplit(run.params[i], split="[.]")[[1]][2], sep="."),
                        paste(code.params[k], wrap.params[i], sep="."), model)
        }
      }

      # This bit is messy...remove if problematic though it will cause issues with parameters.to.save
      param.list <- list("beta.1"=beta.1.str, "beta.2"=beta.2.str, "beta.3"=beta.3.str, "beta.4"=beta.4.str)
      temp <- names(unlist(sapply(param.list, function(x) {
        if(!is.null(x)) {
          if (x %in% c("arm.common", "arm.random", "const.common", "const.random")) {
            x
          }
        }
      })))
      temp <- unlist(sapply(temp, FUN=function(x) strsplit(x, split="[.]")[[1]][2]))
      wrap.params <- append(wrap.params, temp)

    } else {
      wrap.params <- which(sapply(list(beta.1.str, beta.2.str, beta.3.str, beta.4.str),
                                  is.character))
    }

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
      gen.parameters.to.save(model.params=wrap.params, model=model)
  }

  # Add nodes to monitor to calculate plugin pd
  if (pd=="plugin") {
    if (!is.null(rho)) {
      stop("pD cannot be calculated via the plugin method if modelling correlation between time points - is.null(rho)==FALSE")
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

  # Ensure timecor and class are set correctly for getjagsdata
  # if (grepl("rho", model)==TRUE) {
  #   timecor <- TRUE
  # } else {timecor <- FALSE}
  if (length(class.effect)>0) {
    class <- TRUE
  } else {class <- FALSE}


  #### Run jags model ####

  data.ab <- network[["data.ab"]]
  result.jags <- mb.jags(data.ab, model, fun=fun, link=link,
                       class=class, rho=rho, covar=covar, knots=knots,
                       parameters.to.save=parameters.to.save,
                       n.iter=n.iter, n.chains=n.chains,
                       n.burnin=n.burnin, n.thin=n.thin,
                       ...)
  result <- result.jags[["jagsoutput"]]
  jagsdata <- result.jags[["jagsdata"]]

  if (!("error" %in% names(result))) {
    if (pd == "pd.kl" | pd == "popt") {
      if (pd=="pd.kl") {
        temp <- rjags::dic.samples(result$model, n.iter=2000, type="pD")
      } else if (pd=="popt") {
        temp <- rjags::dic.samples(result$model, n.iter=2000, type="popt")
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
                    "fun"=fun, "user.fun"=user.fun,
                    "jagscode"=model, "jagsdata"=jagsdata,
                    "alpha"=alpha,
                    "beta.1"=beta.1, "beta.2"=beta.2,
                    "beta.3"=beta.3, "beta.4"=beta.4,
                    "knots"=knots, "link"=link,
                    "positive.scale"=positive.scale, "intercept"=intercept,
                    "rho"=rho, "covar"=covar,
                    "class.effect"=class.effect, "UME"=UME,
                    "var.scale"=var.scale,
                    "parallel"=parallel, "pd"=pd,
                    "priors"=get.prior(model), "arg.params"=arg.params)
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
                       class=FALSE, rho=NULL, covar=NULL, knots=3,
                       parameters.to.save=parameters.to.save,
                       likelihood=NULL,
                       warn.rhat=FALSE, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertCharacter(model, any.missing=FALSE, len=1, add=argcheck)
  checkmate::assertLogical(class, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertCharacter(parameters.to.save, any.missing=FALSE, unique=TRUE,
                  null.ok=TRUE, add=argcheck)
  checkmate::assertCharacter(fun, any.missing=FALSE, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)


  if (is.null(likelihood)) {
    # For MBNMAtime
    jagsdata <- getjagsdata(data.ab, class=class, rho=rho, covstruct=covar, knots=knots, fun=fun, link=link) # get data into jags correct format (list("fups", "NT", "NS", "narm", "y", "se", "treat", "time"))
  } else if (is.null(rho) & is.null(covar)) {
    # For MBNMAdose
    # jagsdata <- getjagsdata(data.ab, class=class,
    #                         likelihood=likelihood, link=link) # get data into jags correct format
  }


  # Add variable for maxtime to jagsdata if required
  if (grepl("maxtime", model)) {
    jagsdata[["maxtime"]] <- max(data.ab$time)
  } else if (grepl("maxdose", model)) {
    # Only used in MBNMAdose
    #jagsdata[["maxdose"]] <- index.dose(network[["data.ab"]])[["maxdose"]]
  }

  # Remove studyID from jagsdata (not used in model)
  tempjags <- jagsdata
  tempjags[["studyID"]] <- NULL

  # Drop time from tempjags in spline models
  if (fun %in% c("rcs", "ns", "bs")) {
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
  cat(model,file=tmps)
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
#' Identical to `gen.parameters.to.save()` in MBNMAdose
#'
#' @param model.params A character or numeric vector containing the names of the
#' dose-response parameters in the model
#' @param model A JAGS model written as a character object
gen.parameters.to.save <- function(model.params, model) {
  # model.params is a vector (numeric/character) of the names of the dose-response parameters in the model
  #e.g. c(1, 2, 3) or c("emax", "et50")
  # model is a JAGS model written as a character object

  checkmate::assertCharacter(model, len=1)

  model.params <- as.character(model.params)

  # Set some automatic parameters based on the model code
  parameters.to.save <- vector()
  for (i in seq_along(model.params)) {
    # if (grepl(paste0("\\\nd\\.", model.params[i], "\\[(c,)?k\\] ~"), model)==TRUE |
    #     grepl(paste0("\\\nd\\.", model.params[i], "\\[k\\] <- mult\\["), model)==TRUE) {
    if (grepl(paste0("\\\nd\\.", model.params[i]), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("d.", model.params[i]))
    } else if (grepl(paste0("\\\nd\\.", model.params[i]), model)==FALSE) {
      if (grepl(paste0("\\\nbeta\\.", model.params[i], "(\\[k\\])? ~"), model)==TRUE) {
        parameters.to.save <- append(parameters.to.save, paste0("beta.", model.params[i]))
      }
    }
    if (grepl(paste0("\\\nsd\\.", model.params[i], " ~"), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.", model.params[i]))
    }
    if (grepl(paste0("\\\nsd\\.beta.", model.params[i]), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.beta.", model.params[i]))
    }
    if (grepl(paste0("\\\nD\\.", model.params[i], "(\\[k\\])? ~"), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("D.", model.params[i]))
    }
    if (grepl(paste0("\\\nsd\\.D\\.", model.params[i], " ~"), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.D.", model.params[i]))
    }
    if (grepl(paste0("\\\nBETA\\.", model.params[i]), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("BETA.", model.params[i]))
    }
    if (grepl(paste0("\\\nsd\\.BETA\\.", model.params[i]), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.BETA.", model.params[i]))
    }
  }

  for (i in 1:4) {
    if (grepl(paste0("\\\nd\\.", i, " ~"), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("d.", i))
    }
  }

  # For MBNMAtime
  if (grepl("rho", model)==TRUE) {
    parameters.to.save <- append(parameters.to.save, "rho")
  } else {
    parameters.to.save <- append(parameters.to.save, c("totresdev"))
  }

  # For MBNMAdose
  if (grepl("\\\nsd ~", model)==TRUE) {
    parameters.to.save <- append(parameters.to.save, "sd")
  }

  return(unique(parameters.to.save))

}




#' Run MBNMA model with an Emax time-course function and a Hill parameter
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' time-course function. This function accounts for repeated measures over time
#' within studies by applying an Emax time-course function with a Hill
#' parameter. Follows the methods of \insertCite{pedder2019;textual}{MBNMAtime}. This function acts as a
#' wrapper for `mb.run()` that allows for more clearly defined parameter
#' names.
#'
#' @inherit mb.run return references
#' @param network An object of class `mb.network`.
#' @param emax A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' @param et50 A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' @param hill A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' @param ... Arguments to be sent to `mb.run()`
#' @inheritParams mb.run
#'
#' @inheritSection mb.run Time-course parameters
#' @inheritSection mb.run Correlation between observations
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Create mb.network object
#' network <- mb.network(osteopain)
#'
#' # Fit model with random consistency effects on emax, fixed consistency effects on et50, and a
#' #Hill parameter of 0.2 across the network
#' result <- mb.emax.hill(network,
#'   emax=list(pool="rel", method="random"),
#'   et50=list(pool="rel", method="common"),
#'   hill=list(pool="const", method=0.2))
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(result, c("d.et50[2]", "deviance", "sd.emax", "totresdev"))
#'
#' # Traceplots
#' mcmcplots::traplot(result)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(result, "d.emax")
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
#' # Create mb.network object containing class variable
#' network <- mb.network(goutSUA_CFB)
#'
#' # Fit model with class effects and correlation between time points and user-defined priors
#' mb.emax.hill(network, alpha="study",
#'   parameters.to.save=c("d.emax", "d.et50", "beta.hill"),
#'   emax=list(pool="rel", method="random"),
#'   et50=list(pool="rel", method="common"),
#'   hill=list(pool="const", method="common"),
#'   rho="estimate", covar="AR1",
#'   priors=list("rho"="dunif(0,1)"),
#'   class.effect=list("et50"="random")
#'   )
#' }
#' @export
mb.emax.hill <- function(network, fun="emax.hill",
                            emax=list(pool="rel", method="common"),
                            et50=list(pool="rel", method="common"),
                            hill=list(pool="const", method="common"),
                            alpha="study",
                            positive.scale=FALSE, intercept=TRUE, rho=NULL, covar=NULL,
                            var.scale=NULL,
                            class.effect=list(), UME=FALSE,
                            pd="pv", parallel=TRUE,
                            priors=NULL,
                            ...)
{
  arg.params <- list(
    wrap.params=c("emax", "et50", "hill"),
    run.params=c("beta.1", "beta.2", "beta.3")
  )

  # index <- which(sapply(list(emax, et50, hill), is.character))
  #
  # wrap.params=c("emax", "et50", "hill")
  # run.params=c("beta.1", "beta.2", "beta.3")

  # arg.params <- list(
  #   wrap.params=wrap.params[index],
  #   run.params=run.params[index]
  # )

  result <- mb.run(network=network, fun="emax.hill", user.fun=NULL, model.file=NULL,
                      beta.1=emax, beta.2=et50, beta.3=hill, beta.4=NULL,
                      alpha=alpha,
                      positive.scale=positive.scale, intercept=intercept, rho=rho, covar=covar,
                      var.scale=var.scale,
                      class.effect=class.effect, UME=UME,
                      pd=pd, parallel=parallel,
                      priors=priors,
                      arg.params=arg.params, ...
                      )

  return(result)
}



#' Run MBNMA model with an Emax time-course function (without a Hill parameter)
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' time-course function. This function accounts for repeated measures over time
#' within studies by applying an Emax time-course function. Follows the methods
#' of \insertCite{pedder2019;textual}{MBNMAtime}. This function acts as a wrapper for `mb.run()` that
#' allows for more clearly defined parameter names.
#'
#' @inheritParams mb.emax.hill
#' @inherit mb.run return references
#'
#' @inheritSection mb.run Time-course parameters
#' @inheritSection mb.run Correlation between observations
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Create mb.network object
#' network <- mb.network(osteopain)
#'
#' result <- mb.emax(network,
#'   emax=list(pool="rel", method="random"),
#'   et50=list(pool="const", method="common"))
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(result, parms=c("deviance", "sd.emax", "d.emax[2]", "d.emax[3]"))
#'
#' # Traceplots
#' mcmcplots::traplot(result)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(result, "d.emax")
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
#' # Fit model with correlation between time points
#' mb.emax(network, alpha="study",
#'   parameters.to.save=c("d.emax", "d.et50", "rho"),
#'   emax=list(pool="rel", method="random"),
#'   et50=list(pool="rel", method="common"),
#'   rho="estimate", covar="CS"
#'   )
#' }
#' @export
mb.emax <- function(network,
                       emax=list(pool="rel", method="common"),
                       et50=list(pool="rel", method="common"),
                       alpha="study",
                       positive.scale=FALSE, intercept=TRUE, rho=NULL, covar=NULL,
                       var.scale=NULL,
                       class.effect=list(), UME=FALSE,
                       pd="pv", parallel=TRUE,
                       priors=NULL,
                       ...)
{

  #index <- which(sapply(list(emax, et50), is.character))

  # wrap.params=c("emax", "et50")
  # run.params=c("beta.1", "beta.2")

  arg.params <- list(
    wrap.params=c("emax", "et50"),
    run.params=c("beta.1", "beta.2")
  )

  # arg.params <- list(
  #   wrap.params=wrap.params[index],
  #   run.params=run.params[index]
  # )


  result <- mb.run(network=network, fun="emax", user.fun=NULL, model.file=NULL,
                      beta.1=emax, beta.2=et50, beta.3=NULL, beta.4=NULL,
                      alpha=alpha,
                      positive.scale=positive.scale, intercept=intercept, rho=rho, covar=covar,
                      var.scale=var.scale,
                      class.effect=class.effect, UME=UME,
                      pd=pd, parallel=parallel,
                      priors=priors,
                      arg.params=arg.params, ...
  )

  return(result)
}




#' Run MBNMA model with an exponential time-course function
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' time-course function. This function accounts for repeated measures over time
#' within studies by applying an exponential time-course function. Follows the
#' methods of \insertCite{pedder2019;textual}{MBNMAtime}. This function acts as a wrapper for `mb.run()`
#' that allows for more clearly defined parameter names.
#'
#' @inheritParams mb.emax.hill
#' @inherit mb.run return references
#' @param lambda A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#'
#' @inheritSection mb.run Time-course parameters
#' @inheritSection mb.run Correlation between observations
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Create mb.network object
#' network <- mb.network(osteopain)
#'
#' # Fit exponential time-course with random consistency treatment effects
#' result <- mb.exponential(network,
#'   lambda=list(pool="rel", method="random"))
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(result, parms=c("sd.lambda", "d.lambda[2]"))
#'
#' # Traceplots
#' mcmcplots::traplot(result)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(result, "d.lambda")
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
#' # Fit model with unrelated mean effects that saves residual deviance contributions
#' mb.exponential(network, alpha="study",
#'   parameters.to.save=c("d.lambda", "resdev"),
#'   lambda=list(pool="rel", method="random"),
#'   UME=TRUE
#'   )
#' }
#' @export
mb.exponential <- function(network, lambda=list(pool="rel", method="common"),
                              alpha="study",
                              positive.scale=FALSE, intercept=TRUE, rho=NULL, covar=NULL,
                              var.scale=NULL,
                              class.effect=list(), UME=FALSE,
                              pd="pv", parallel=TRUE,
                              priors=NULL,
                              ...)
{
  arg.params <- list(
    wrap.params=c("lambda"),
    run.params=c("beta.1")
  )

  # index <- which(sapply(list(lambda), is.character))
  #
  # wrap.params=c("lambda")
  # run.params=c("beta.1")
  #
  # arg.params <- list(
  #   wrap.params=wrap.params[index],
  #   run.params=run.params[index]
  # )


  result <- mb.run(network=network, fun="exponential", user.fun=NULL, model.file=NULL,
                      beta.1=lambda, beta.2=NULL, beta.3=NULL, beta.4=NULL,
                      alpha=alpha,
                      positive.scale=positive.scale, intercept=intercept, rho=rho, covar=covar,
                      var.scale=var.scale,
                      class.effect=class.effect, UME=UME,
                      pd=pd, parallel=parallel,
                      priors=priors,
                      arg.params=arg.params, ...
  )

  return(result)
}





#' Run MBNMA model with a first-order fractional polynomial time-course function
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' time-course function. This function accounts for repeated measures over time
#' within studies by applying a first-order fractional polynomial time-course
#' function \insertCite{jansen2015}{MBNMAtime}. Follows the methods of \insertCite{pedder2019;textual}{MBNMAtime}. This function acts as a
#' wrapper for `mb.run()` that allows for more clearly defined parameter
#' names.
#'
#' @inheritParams mb.emax.hill
#' @inherit mb.run return references
#' @param slope list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' @param power A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' For this parameter, `pool` must be set to `"arm"` or `"const"` (i.e. it cannot
#' be `"rel"`).
#'
#' @inheritSection mb.run Time-course parameters
#' @inheritSection mb.run Correlation between observations
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \donttest{
#' # Create mb.network object
#' network <- mb.network(osteopain)
#'
#' # Fit 1st order fractional polynomial time-course with random consistency treatment effects
#' #on the slope and a common parameter for power estimated across the network
#' result <- mb.fract.first(network,
#'   slope=list(pool="rel", method="random"),
#'   power=list(pool="const", method="common"))
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(result, parms=c("sd.slope", "deviance", "beta.power"))
#'
#' # Traceplots
#' mcmcplots::traplot(result)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(result, "d.slope")
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
#' # Fit model with unrelated mean effects and correlation between time points
#' mb.fract.first(network, alpha="study",
#'   slope=list(pool="rel", method="random"),
#'   power=list(pool="const", method="common"),
#'   rho=0.5, covar="AR1",
#'   UME="slope"
#'   )
#' }
#' @export
mb.fract.first <- function(network, slope=list(pool="rel", method="common"),
                              power=list(pool="const", method="common"),
                              alpha="study",
                              positive.scale=FALSE, intercept=TRUE, rho=NULL, covar=NULL,
                              var.scale=NULL,
                              class.effect=list(), UME=FALSE,
                              pd="pv", parallel=TRUE,
                              priors=NULL,
                              ...)
{
  # Checks- TO ADD

  arg.params <- list(
    wrap.params=c("slope", "power"),
    run.params=c("beta.1", "beta.3")
  )

  # index <- which(sapply(list(slope, power), is.character))
  #
  # wrap.params=c("slope", "power")
  # run.params=c("beta.1", "beta.3")
  #
  # arg.params <- list(
  #   wrap.params=wrap.params[index],
  #   run.params=run.params[index]
  # )

  result <- mb.run(network=network, fun="fract.poly.first", user.fun=NULL, model.file=NULL,
                      beta.1=slope, beta.2=NULL, beta.3=power, beta.4=NULL,
                      alpha=alpha,
                      positive.scale=positive.scale, intercept=intercept, rho=rho, covar=covar,
                      var.scale=var.scale,
                      class.effect=class.effect, UME=UME,
                      pd=pd, parallel=parallel,
                      priors=priors,
                      arg.params=arg.params, ...
  )

  return(result)
}









#' Run MBNMA model with a second-order fractional polynomial time-course
#' function
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' time-course function. This function accounts for repeated measures over time
#' within studies by applying a second-order fractional polynomial time-course
#' function \insertCite{jansen2015}{MBNMAtime}. Follows the methods of
#' \insertCite{pedder2019;textual}{MBNMAtime}. This function acts as a wrapper for
#' `mb.run()` that allows for more clearly defined parameter names.
#'
#' @inheritParams mb.emax.hill
#' @inherit mb.run return references
#' @param slope.1 A list with named elements `pool` and `method` that refers to
#'   time-course parameter(s) specified within the time-course function (see
#'   details).
#' @param slope.2 A list with named elements `pool` and `method` that refers to
#'   time-course parameter(s) specified within the time-course function (see
#'   details).
#' @param power.1 A list with named elements `pool` and `method` that refers to
#'   time-course parameter(s) specified within the time-course function (see
#'   details). For this parameter, `pool` must be set to `"arm"` or `"const"`
#'   (i.e. it cannot be `"rel"`).
#' @param power.2 A list with named elements `pool` and `method` that refers to
#'   time-course parameter(s) specified within the time-course function (see
#'   details). For this parameter, `pool` must be set to `"arm"` or `"const"`
#'   (i.e. it cannot be `"rel"`).
#'
#' @inheritSection mb.run Time-course parameters
#' @inheritSection mb.run Correlation between observations
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \donttest{
#' # Create mb.network object
#' network <- mb.network(osteopain)
#'
#' # Fit 2nd order fractional polynomial time-course with fixed consistency treatment effects on
#' #beta.1 and beta.2, absolute time-course parameters estimated by treatment for power.1,
#' #and an exchangeable parameter for power estimated across the network
#' result <- mb.fract.second(network,
#'   slope.1=list(pool="rel", method="common"),
#'   slope.2=list(pool="rel", method="common"),
#'   power.1=list(pool="arm", method="common"),
#'   power.2=list(pool="arm", method="random"))
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(result, parms=c("beta.power.1", "beta.power.2"))
#'
#' # Traceplots
#' mcmcplots::traplot(result)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(result, "beta.power.1")
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
#' # Fit model with correlation between time points
#' mb.fract.second(network, alpha="study",
#'   slope.1=list(pool="rel", method="common"),
#'   slope.2=list(pool="rel", method="common"),
#'   power.1=list(pool="const", method="common"),
#'   power.2=list(pool="const", method="common"),
#'   rho=0.5, covar="AR1"
#'   )
#' }
#' @export
mb.fract.second <- function(network, slope.1=list(pool="rel", method="common"),
                               slope.2=list(pool="rel", method="common"),
                               power.1=list(pool="const", method="common"),
                               power.2=list(pool="const", method="common"),
                               alpha="study",
                               positive.scale=FALSE, intercept=TRUE, rho=NULL, covar=NULL,
                               var.scale=NULL,
                               class.effect=list(), UME=FALSE,
                               pd="pv", parallel=TRUE,
                               priors=NULL,
                               ...)
{
  # Checks - TO ADD

  arg.params <- list(
    wrap.params=c("slope.1", "slope.2", "power.1", "power.2"),
    run.params=c("beta.1", "beta.2", "beta.3", "beta.4")
  )

  # index <- which(sapply(list(power.1, power.2), is.character))
  #
  # wrap.params=c("power.1", "power.2")
  # run.params=c("beta.3", "beta.4")
  #
  # arg.params <- list(
  #   wrap.params=wrap.params[index],
  #   run.params=run.params[index]
  # )


  result <- mb.run(network=network, fun="fract.poly.second", user.fun=NULL, model.file=NULL,
                      beta.1=slope.1, beta.2=slope.2, beta.3=power.1, beta.4=power.2,
                      alpha=alpha,
                      positive.scale=positive.scale, intercept=intercept, rho=rho, covar=covar,
                      var.scale=var.scale,
                      class.effect=class.effect, UME=UME,
                      pd=pd, parallel=parallel,
                      priors=priors,
                      arg.params=arg.params, ...
  )

  return(result)
}





#' Run MBNMA model with a linear time-course function
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' time-course function. This function accounts for repeated measures over time
#' within studies by applying a linear time-course function. Follows the methods
#' of \insertCite{pedder2019;textual}{MBNMAtime}. This function acts as a wrapper for `mb.run()` that
#' allows for more clearly defined parameter names.
#'
#' @inheritParams mb.emax.hill
#' @inherit mb.run return references
#' @param slope A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#'
#' @inheritSection mb.run Time-course parameters
#' @inheritSection mb.run Correlation between observations
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Create mb.network object
#' network <- mb.network(osteopain)
#'
#' # Fit linear time-course with random consistency treatment effects on the slope
#' result <- mb.linear(network,
#'   slope=list(pool="rel", method="random"))
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(result, parms="sd.slope")
#'
#' # Traceplots
#' mcmcplots::traplot(result)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(result, "d.slope")
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
#' # Fit model with unrelated mean effects and correlation between time points
#' mb.linear(network, alpha="study",
#'   slope=list(pool="rel", method="random"),
#'   rho=0.5, covar="AR1",
#'   UME=TRUE
#'   )
#' }
#' @export
mb.linear <- function(network, slope=list(pool="rel", method="common"),
                         alpha="study",
                         positive.scale=FALSE, intercept=TRUE, rho=NULL, covar=NULL,
                         var.scale=NULL,
                         class.effect=list(), UME=FALSE,
                         pd="pv", parallel=TRUE,
                         priors=NULL,
                         ...)
{

  arg.params <- list(
    wrap.params=c("slope"),
    run.params=c("beta.1")
  )

  # index <- which(sapply(list(slope), is.character))
  #
  # wrap.params=c("slope")
  # run.params=c("beta.1")
  #
  # arg.params <- list(
  #   wrap.params=wrap.params[index],
  #   run.params=run.params[index]
  # )

  result <- mb.run(network=network, fun="linear", user.fun=NULL, model.file=NULL,
                      beta.1=slope, beta.2=NULL, beta.3=NULL, beta.4=NULL,
                      alpha=alpha,
                      positive.scale=positive.scale, intercept=intercept, rho=rho, covar=covar,
                      var.scale=var.scale,
                      class.effect=class.effect, UME=UME,
                      pd=pd, parallel=parallel,
                      priors=priors,
                      arg.params=arg.params, ...
  )

  return(result)
}






#' Run MBNMA model with a quadratic time-course function
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' time-course function. This function accounts for repeated measures over time
#' within studies by applying a quadratic time-course function. Follows the
#' methods of \insertCite{pedder2019;textual}{MBNMAtime}. This function acts as a wrapper for `mb.run()`
#' that allows for more clearly defined parameter names.
#'
#' @inheritParams mb.emax.hill
#' @inherit mb.run return references
#' @param beta.1 A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' @param beta.2 A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#'
#' @inheritSection mb.run Time-course parameters
#' @inheritSection mb.run Correlation between observations
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Create mb.network object
#' network <- mb.network(osteopain)
#'
#' # Fit quadratic time-course with fixed consistency treatment effects on beta.1 and
#' #random consistency treatment effects on beta.2
#' result <- mb.quadratic(network,
#'   beta.1=list(pool="rel", method="common"),
#'   beta.2=list(pool="rel", method="random"))
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(result, parms=c("sd.2", "d.1[3]", "d.2[3]", "totresdev"))
#'
#' # Traceplots
#' mcmcplots::traplot(result)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(result, "d.2")
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
#' # Fit model with unrelated mean effects on beta.1
#' mb.quadratic(network, alpha="study",
#'   beta.1=list(pool="rel", method="random"),
#'   beta.2=list(pool="const", method="common"),
#'   UME="beta.1"
#'   )
#' }
#' @export
mb.quadratic <- function(network, beta.1=list(pool="rel", method="common"),
                            beta.2=list(pool="rel", method="common"),
                            alpha="study",
                            positive.scale=FALSE, intercept=TRUE, rho=NULL, covar=NULL,
                            var.scale=NULL,
                            class.effect=list(), UME=FALSE,
                            pd="pv", parallel=TRUE,
                            priors=NULL,
                            ...)
{
  result <- mb.run(network=network, fun="quadratic", user.fun=NULL, model.file=NULL,
                      beta.1=beta.1, beta.2=beta.2, beta.3=NULL, beta.4=NULL,
                      alpha=alpha,
                      positive.scale=positive.scale, intercept=intercept, rho=rho, covar=covar,
                      var.scale=var.scale,
                      class.effect=class.effect, UME=UME,
                      pd=pd, parallel=parallel,
                      priors=priors,
                      ...
  )

  return(result)
}






#' Run MBNMA model with a piecewise linear time-course function
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' time-course function. This function accounts for repeated measures over time
#' within studies by applying a piecewise linear time-course function. Follows
#' the methods of \insertCite{pedder2019;textual}{MBNMAtime}. This function acts as a wrapper for
#' `mb.run()` that allows for more clearly defined parameter names.
#'
#' @inheritParams mb.emax.hill
#' @inherit mb.run return references
#' @param slope.1 A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' @param slope.2 A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' @param knot A list with named elements `pool` and `method` that refers to
#' time-course parameter(s) specified within the time-course function (see details).
#' For this parameter, `pool` must be set to `"arm"` or `"const"` (i.e. it cannot
#' be `"rel"`).
#'
#' @inheritSection mb.run Time-course parameters
#' @inheritSection mb.run Correlation between observations
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Create mb.network object
#' network <- mb.network(osteopain)
#'
#' # Fit piecewise linear time-course with random consistency treatment effects on slope.1,
#' #absolute time-course parameters by treatment on slope.2, and an exchangeable parameter
#' #for knot estimated across the network
#' result <- mb.piecelinear(network,
#'   slope.1=list(pool="rel", method="random"),
#'   slope.2=list(pool="arm", method="common"),
#'   knot=list(pool="const", method="common"))
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(result, parms=c("beta.slope.2", "beta.knot"))
#'
#' # Traceplots
#' mcmcplots::traplot(result)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(result, "beta.slope.2")
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
#' # Fit model with unrelated mean effects on beta.1 and correlation between time points
#' mb.piecelinear(network, alpha="study",
#'   slope.1=list(pool="rel", method="random"),
#'   slope.2=list(pool="rel", method="common"),
#'   knot=list(pool="const", method=1),
#'   rho="estimate", covar="CS",
#'   UME="slope.1"
#'   )
#' }
#' @export
mb.piecelinear <- function(network, slope.1=list(pool="rel", method="common"), slope.2=list(pool="rel", method="common"), knot=list(pool="const", method="common"),
                              alpha="study",
                              positive.scale=FALSE, intercept=TRUE, rho=NULL, covar=NULL,
                              var.scale=NULL,
                              class.effect=list(), UME=FALSE,
                              pd="pv", parallel=TRUE,
                              priors=NULL,
                              ...)
{

  arg.params <- list(
    wrap.params=c("slope.1", "slope.2", "knot"),
    run.params=c("beta.1", "beta.2", "beta.3")
  )

  # index <- which(sapply(list(slope.1, slope.2, knot), is.character))
  #
  # wrap.params=c("slope.1", "slope.2", "knot")
  # run.params=c("beta.1", "beta.2", "beta.3")
  #
  # arg.params <- list(
  #   wrap.params=wrap.params[index],
  #   run.params=run.params[index]
  # )

  result <- mb.run(network=network, fun="piecelinear", user.fun=NULL, model.file=NULL,
                      beta.1=slope.1, beta.2=slope.2, beta.3=knot, beta.4=NULL,
                      alpha=alpha,
                      positive.scale=positive.scale, intercept=intercept, rho=rho, covar=covar,
                      var.scale=var.scale,
                      class.effect=class.effect, UME=UME,
                      pd=pd, parallel=parallel,
                      priors=priors,
                      arg.params=arg.params, ...
  )

  return(result)
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
#'   Alternative approaches are to use pV (`pv`) as an approximation (Plummer REF) or
#'   pD calculated by Kullback–Leibler divergence (`pd.kl`) or using an optimism adjustment (`popt`) (REF).
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
#' emax <- mb.emax(network, parameters.to.save=c("theta", "resdev"))
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
#' emax <- mb.emax(network)
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

  if (!(grepl(paste0("\n", param), mbnma$model.arg$jagscode)==TRUE |
        grepl(paste0(" ", param), mbnma$model.arg$jagscode)==TRUE)) {
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
