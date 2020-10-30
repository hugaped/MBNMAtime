# Functions for writing MBNMA models
# Author: Hugo Pedder
# Date created: 2018-09-10


#' Write MBNMA time-course models JAGS code
#'
#' Writes JAGS code for a Bayesian time-course model for model-based network
#' meta-analysis (MBNMA).
#'
#' @inheritParams mb.run
#'
#' @return A single long character string containing the JAGS model generated
#'   based on the arguments passed to the function.
#'
#' @inherit mb.run details
#'
#' @examples
#' # Write a linear time-course MBNMA with random treatment effects on beta.1 and equal baselines
#' #in study arms
#' model <- mb.write(fun="linear", alpha="study", beta.1="rel.random")
#' cat(model) # Concatenates model representations making code more easily readable
#'
#' # Write an emax time-course MBNMA with a Hill parameter of 0.5 with no intercept
#' model <- mb.write(fun="emax.hill",
#'   beta.1="rel.common", beta.2="const.common", beta.3=0.5,
#'   intercept=TRUE)
#' cat(model) # Concatenates model representations making code more easily readable
#'
#' # Write an exponential time-course MBNMA that accounts for correlation between time points
#' model <- mb.write(fun="exponential",
#'   alpha="arm", beta.1="rel.common",
#'   rho="estimate", covar="AR1")
#' cat(model)
#'
#' # Define a user-defined time-course relationship for the MBNMA JAGS model
#' time.fun <- ~ alpha + (exp(beta.1 * time) / (beta.2 * time))
#' model <- mb.write(fun="user", user.fun=time.fun,
#'   beta.1="rel.random", beta.2="rel.common")
#' cat(model)
#' @export
mb.write <- function(fun="linear", user.fun=NULL, alpha="study", beta.1="rel.common", beta.2=NULL, beta.3=NULL, beta.4=NULL,
                        positive.scale=TRUE, intercept=TRUE, rho=NULL, covar=NULL, knots=3,
                        var.scale=NULL,
                        class.effect=list(), UME=FALSE) {


  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertFormula(user.fun, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(alpha, choices=c("arm", "study"), null.ok=FALSE, add=argcheck)
  checkmate::assertLogical(positive.scale, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertLogical(intercept, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  #checkmate::assertLogical(timecor, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertChoice(covar, choices=c("CS", "AR1"), null.ok=TRUE, add=argcheck)
  checkmate::assertList(class.effect, unique=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)


  ####### VECTORS #######

  parameters <- c("alpha", "beta.1", "beta.2", "beta.3", "beta.4")

  # Change UME to relate to parameters in model
  if (length(UME)==1) {
    if (UME==TRUE) {
      UME <- vector()
      treat.params <- parameters[2:5]
      for (i in seq_along(treat.params)) {
        if (!is.null(get(treat.params[i]))) {
          if (get(treat.params[i]) %in% c("rel.common","rel.random")) {
            UME <- append(UME, treat.params[i])
          }
        }
      }
    }
  }

  write.check(fun=fun, user.fun=user.fun, alpha=alpha,
                                      beta.1=beta.1, beta.2=beta.2, beta.3=beta.3, beta.4=beta.4,
                                      positive.scale=positive.scale, intercept=intercept,
                                      rho=rho, covar=covar, var.scale=var.scale, knots=knots,
                                      class.effect=class.effect, UME=UME)

  model <- write.model()

  timecourse <- time.fun(fun=fun, user.fun=user.fun, knots=knots,
                         alpha=alpha, beta.1=beta.1, beta.2=beta.2,
                         beta.3=beta.3, beta.4=beta.4)[["jagscode"]]

  alphacode <- write.alpha(model, timecourse,
                       intercept=intercept, positive.scale=positive.scale)

  timecourse <- alphacode[["timecourse"]]
  model <- alphacode[["model"]]

  model <- write.likelihood(model, timecourse, rho, covar)

  #model <- write.fract.poly(model, timecourse)

  model <- write.beta(model, timecourse,
                      beta.1=beta.1, beta.2=beta.2, beta.3=beta.3, beta.4=beta.4,
                      UME=UME, class.effect=class.effect)

  # Ammend fract.poly for s.beta changes
  #model <- write.fract.poly(model, timecourse)

  #model <- write.piecelinear(model, beta.1, beta.3)

  model <- write.piece.fract(model, fun, beta.3)

  model <- write.cor(model, var.scale=var.scale, class.effect = class.effect)

  model <- write.remove.loops(model)

  return(model)
}










#' Write time-course function component of JAGS code for MBNMA time-course
#' models
#'
#' Writes a single line of JAGS code representing the time-course function
#' component of an MBNMA time-course model, outputted as a single character
#' string.
#'
#' @inheritParams mb.run
#'
#' @return A single character string containing JAGS model representing the
#'   time-course function component of an MBNMA time-course model, generated
#'   based on the arguments passed to the function.
#'
#' @inherit mb.run details
#'
time.fun <- function(fun="linear", user.fun=NULL, alpha="arm", knots=3,
                     beta.1="rel.common", beta.2=NULL, beta.3=NULL, beta.4=NULL) {
  # fun can be any of the following functions: linear, exponential, piecewise linear, Emax, Emax with Hill parameter, fractional polynomial (Jansen 2015)
  # user.fun can be any user function combining alpha, beta.1 / beta.2 / beta.3, gamma.1 and gamma.2, and time
  # alpha can be:
  # "study" to constrain baseline to be equal for all arms within a study ([i] index is added)
  # "arm" to allow baseline to vary between arms within a study ([i,k] index is added)
  # beta.1, beta.2, beta.3 and beta.4 can be:
  # "rel.common" to allow it to be split into reference and relative fixed treatment effects (for estimating consistency relationships using indirect evidence)
  # "rel.random" to allow it to be split into reference and relative random treatment effects (for estimating consistency relationships using indirect evidence)
  # "common" to estimate a single parameter across the network that is independent of treatment
  # Given a numeric value to allow user to provide data for a single parameter across the network that is independent of treatment

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertFormula(user.fun, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (fun=="linear") {
    timecourse <- "alpha + (beta.1 * time)"
  } else if (fun=="exponential") {
    #timecourse <- "alpha + exp(beta.1 * time)"
    timecourse <- "alpha + (beta.1 * (1 - exp(- time)))"
  } else if (fun=="emax") {
    timecourse <- "alpha + ((beta.1 * time) / (exp(beta.2) + time))"
    message("ET50 parameters (beta.2) are on exponential scale to ensure they take positive values on the natural scale")
  } else if (fun=="emax.hill") {
    timecourse <- "alpha + ((beta.1 * (time^exp(beta.3))) / ((exp(beta.2)^exp(beta.3)) + (time^exp(beta.3))))"
    message("ET50 parameters (beta.2) are on exponential scale to ensure they take positive values on the natural scale")
    message("Hill parameters (beta.3) are on exponential scale to ensure they take positive values on the natural scale")
  } else if (fun=="quadratic") {
    timecourse <- "alpha + (beta.1 * time) + (beta.2 * (time^2))"
  } else if (fun=="fract.poly.first") {
    #timecourse <- "alpha + (beta.1 * time.fp.1)"
    timecourse <- "alpha + (beta.1 * ifelse(time>0, ifelse(beta.3==0, log(time), time^beta.3), 0) )"
  } else if (fun=="fract.poly.second") {
    #timecourse <- "alpha + (beta.1 * time.fp.1) + (beta.2 * time.fp.2)"
    timecourse <- "alpha + (beta.1 * ifelse(time>0, ifelse(beta.3==0, log(time), time^beta.3), 0)) + (beta.2 * ifelse(beta.4==beta.3, ifelse(time>0, ifelse(beta.4==0, log(time)^2, (time^beta.4) * log(time)), 0), ifelse(time>0, ifelse(beta.4==0, log(time), time^beta.4), 0)))"
  } else if (fun=="piecelinear") {
    #timecourse <- "((time < beta.3) * alpha) + ((time < beta.3) * (beta.1 * time)) + ((time >= beta.3) * alpha.knot) + ((time >= beta.3) * (beta.2 * time))"
    timecourse <- "((time < beta.3) * alpha) + ((time < beta.3) * (beta.1 * time)) + ((time >= beta.3) * (alpha + (beta.1 * beta.3))) + ((time >= beta.3) * (beta.2 * time))"
  } else if (any(c("rcs", "ns", "bs") %in% fun)) {
    knotnum <- ifelse(length(knots)>1, length(knots), knots)
    timetemp <- "alpha + "
    for (knot in 1:(knotnum-1)) {
      timetemp <- paste0(timetemp, "(beta.", knot, " * spline[i,m,", knot, "])")
      if (knot<knotnum-1) {
        timetemp <- paste0(timetemp, " + ")
      }
    }
    timecourse <- timetemp
  }

  if (fun=="user") {
    user.str <- as.character(user.fun[2])
    timecourse <- user.str
  }


  # Checks that beta parameters not required are NULL for non fract.poly
  if (!(fun=="fract.poly.first" | fun=="fract.poly.second")) {
    betaparams <- c("beta.1", "beta.2", "beta.3", "beta.4")
    for (i in seq_along(betaparams)) {
      if (grepl(betaparams[i], timecourse)==TRUE) {
        if (is.null(get(betaparams[i]))) {
          stop(paste0(betaparams[i], " has not been specified and is required for specified timecourse"))
        }
      } else if (grepl(betaparams[i], timecourse)==FALSE) {
        if (!is.null(get(betaparams[i]))) {
          stop(paste0(betaparams[i], " has been specified for the selected time-course when it is not necessary"))
        }
      }
    }
  }

  # Checks that beta.3 is not fixed or random for piecelinear
  if (fun=="piecelinear") {
    if (!is.numeric(beta.3) & (beta.3=="rel.common" | beta.3=="rel.random")) {
      stop("`beta.3` is a knot parameter and cannot have consistency equations placed on it. It must either be a numeric value or assigned `common` if it is to be estimated from the data.")
    }
    if (is.numeric(beta.3) & beta.3<=0) {
      stop("`beta.3` is a parameter indicating the location of a piecewise knot and must take a positive value.")
    }
  }


  jagscode <- timecourse

  jagscode <- gsub("time", "time[i,m]", jagscode, perl=T)

  if (alpha=="study") {
    jagscode <- gsub("alpha(?!\\.)", "alpha[i]", jagscode, perl=T)
    #jagscode <- gsub("alpha\\.1", "alpha.1[i]", jagscode, perl=T) # For piecewise linear
  } else if (alpha=="arm") {
    jagscode <- gsub("alpha(?!\\.)", "alpha[i,k]", jagscode, perl=T)
    #jagscode <- gsub("alpha\\.1", "alpha.1[i,k]", jagscode, perl=T) # For piecewise linear
  }

  for (i in 1:4) {
    if (!is.null(get(paste0("beta.", i)))) {
      if (grepl("rel.common", get(paste0("beta.", i)))==TRUE |
                grepl("rel.random", get(paste0("beta.", i)))==TRUE) {
        jagscode <- gsub(paste0("beta\\.", i), paste0("beta.", i, "[i,k]"),
                         jagscode, perl=T)
      }
    }
  }


  return(list("time.function"=fun, "relationship"=timecourse, "jagscode"=jagscode))
}













#' Checks validity of arguments for mb.write
#'
#' @inheritParams mb.run
#'
#' @return A boolean object that indicates whether the arguments imply modelling
#'   correlation between time points.
#'
#' @details Used to check if the arguments given to mb.write are valid. The
#'   function will return informative errors if arguments are mispecified and
#'   will return an object that indicates whether the arguments imply modelling a
#'   correlation between time points if it passes.
#'
write.check <- function(fun="linear", user.fun=NULL, alpha="arm", beta.1="rel.common", beta.2=NULL, beta.3=NULL, beta.4=NULL,
                        positive.scale=TRUE, intercept=TRUE, rho=NULL, covar=NULL, knots=3,
                        var.scale=NULL,
                        class.effect=list(), UME=FALSE) {
  parameters <- c("alpha", "beta.1", "beta.2", "beta.3", "beta.4")

  # Run argument checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(alpha, choices=c("arm", "study"), null.ok=FALSE, add=argcheck)
  checkmate::assertList(class.effect, unique=FALSE, add=argcheck)
  checkmate::assertNumeric(var.scale, null.ok = TRUE)
  checkmate::assertNumeric(knots, null.ok=FALSE, add=argcheck)

  timefuns <- c("none", "linear", "quadratic", "exponential", "emax", "emax.hill",
                "fract.poly.first", "fract.poly.second", "piecelinear", "user", "rcs", "bs", "ns")
  checkmate::assertChoice(fun, choices=timefuns,
                          null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check knots
  knoterr <- "Minimum number of `knots` for fun=`rcs` is 3"
  if (length(knots)==1) {
    if (knots<3) {
      stop(knoterr)
    }
  } else if (length(knots)>1 & length(knots)<3) {
    stop(knoterr)
  }
  if (length(knots)>1) {
    if (!(all(knots<=1 & all(knots>=0)))) {
      stop("`knots` specified as quantiles must be between 0 and 1")
    }
  }


  # Checks that beta parameters are single values
  for (i in 1:4) {
    betaparam <- get(paste0("beta.", i))
    if (!is.null(betaparam)) {
      if (length(betaparam)>1) {
        stop("Time-course parameters cannot be assigned a vector. They can be only `rel.common`, `rel.random`, `abs.i`, `arm.common`, `arm.random`, `const.common`, `const.random`, or given a numeric value.")
      }

      if (is.character(betaparam) & !(betaparam %in% c("rel.common", "rel.random", "abs.i", "arm.common", "arm.random", "const.common", "const.random"))) {
        msg <- paste0("A time-course parameter has currently been given the value: `", betaparam,
                      "`. Time-course parameters can only be `rel.common`, `rel.random`, `abs.i`, `arm.common`, `arm.random`, `const.common`, `const.random`, or given a numeric value.")
        stop(msg)
      }
    }
  }

  # Checks for parameter classifications
  if (fun=="emax.hill" & is.null(beta.3)) {
    stop("Hill parameter (beta.3) for emax.hill has not been specified.")
  }

  if (fun=="fract.poly.first" | fun=="fract.poly.second") {
    if (!is.null(beta.3)) {
      if (!(beta.3 %in% c("arm.common", "arm.random", "const.common", "const.random") |
            is.numeric(beta.3))) {
        stop("For fractional polynomials, power.1 (beta.3) and power.2 (beta.4) represent powers of time and must be estimated as a single absolute value across the whole network or given as data (`pool=`const`)")
      }
    }
  }

  if (fun=="fract.poly.second") {
    if (!is.null(beta.4)) {
      if (!(beta.4 %in% c("arm.common", "arm.random", "const.common", "const.random") |
            is.numeric(beta.4))) {
        stop("For fractional polynomials, power.1 (beta.3) and power.2 (beta.4) represent powers of time and must be estimated as a single absolute value across the whole network or given as data (`pool=`const`)")
      }
    }
  }

  if (fun=="piecelinear") {
    if (!is.null(beta.3)) {
      if (beta.3 %in% c("rel.common", "rel.random")) {
        stop("Relative effects cannot be used for pooling for knot (beta.3) in piecewise linear function - set `pool` to `arm` or `const` for this parameter")
      }
    } else if (is.null(beta.3)) {
      stop("knot (beta.3) in piecewise linear function cannot be NULL - `pool` must be set to `arm` or `const` for this parameter")
    }
  }

  #### Check user.fun ####
  if (fun!="user" & !is.null(user.fun)) {
    warning(paste0("user.fun is only applied if fun=`user`. Time-course function used for this model will be ", fun))
  }

  if (fun=="user") {
    if (is.null(user.fun)) {
      stop("user.fun must contain a function that includes a combination of alpha and beta parameters with time")
    }

    user.str <- as.character(user.fun[2])
    if (grepl("beta.2", user.str)==TRUE & grepl("beta.1", user.str)==FALSE) {
      stop("user.fun cannot contain beta.2 if beta.1 is not present")
    } else if (grepl("beta.3", user.str)==TRUE & grepl("beta.2", user.str)==FALSE | grepl("beta.1", user.str)==FALSE) {
      stop("user.fun cannot contain beta.3 if beta.2 and beta.1 are not present")
    }

    for (i in 1:4) {
      if (grepl(paste0("beta.",i), user.str)==TRUE) {
        if(is.null(get(paste0("beta.",i)))) {
          msg <- paste0("beta.",i, " has been specified in `user.fun` time-course function yet no arguments have been given for it")
          stop(msg)
        }
      }
    }
  }

  if (!is.null(rho) & is.null(covar)) {
    stop("Arguments for `rho` imply modelling correlation betwen time points, which requires a covariance structure to be specified using `covar`.")
  }

  # if ((!is.null(rho) & rho!="dunif(0,1)") & timecor==FALSE) {
  #   warning("Arguments for rho imply modelling a correlation between time points. timecor has been set to TRUE to allow for this")
  #   timecor <- TRUE
  # }

  if (!is.null(rho)) {
    if (is.character(rho) & rho!="estimate") {
      stop("`rho` must either be assigned the value `estimate` to be estimated from the data or must be assigned a single numeric value")
    }
    if (is.numeric(rho) & (rho < -1 | rho > 1)) {
      stop("Numeric values for `rho` cannot be outside the bounds of c(-1, 1)")
    }
    if (is.null(covar)) {
      stop("`rho` has been assigned a value. Must therefore also define `covar`")
    }
  }

  if (length(class.effect)>0) {
    for (i in seq_along(class.effect)) {

      if (is.null(names(class.effect)[i])) {
        stop(paste0("Elements of class.effect must be named with the names of any of the following parameters: ", paste(parameters[2:5], collapse=", ")))
      }
      if (!(names(class.effect)[i] %in% parameters[2:5])) {
        stop(paste0("Can only specify class effects for parameters which have been assigned consistency relationships and given the value `common` or `random`"))
      }
      if (!(get(names(class.effect)[i]) %in% c("rel.common", "rel.random", "arm.common", "arm.random"))) {
        stop("Can only specify class effects for parameters which have been assigned consistency relationships and given the value `common` or `random`")
      }

      for (k in seq_along(class.effect[[i]])) {
        if (length(class.effect[[i]])!=1) {
          stop("Each class effect can only take a single value which must be either `common` or `random`")
        }
        if (!(class.effect[[i]]=="common" | class.effect[[i]]=="random")) {
          stop("Each class effect can only take a single value which must be either `common` or `random`")
        }
      }
    }

    # Check that UME and class effects are not on same parameter
    for (i in seq_along(names(class.effect))) {
      if (names(class.effect)[i] %in% UME) {
        stop(paste0("UME / node-splitting model cannot be applied to parameters with class effects. ", names(class.effect)[i], " has been assigned as a class effect and as a parameter for node-splitting / UME"))
      }
    }
  }

  if (!all(UME==FALSE)) {
    for (i in seq_along(UME)) {
      if (!(get(UME[i]) %in% c("rel.common", "rel.random"))) {
        stop("Can only specify UME for parameters which have been assigned consistency relationships and given the value `rel.common` or `rel.random`")
      }
    }
  }

}




#' Write the basic JAGS model code for MBNMA to which other lines of model
#' code can be added
#'
#' @return A character object of JAGS model code
#'
write.model <- function() {
  model <-
"
model{ 			# Begin Model Code

for(i in 1:NS){ # Run through all NS trials

for (k in 1:narm[i]){ # Run through all arms within a study

for (m in 1:fups[i]) {	# Run through all observations within a study

resdev[i,k,m] <- pow((y[i,k,m] - theta[i,k,m]),2) * prec[i,k,m] # residual deviance for normal likelihood
dev[i,k,m] <- -2* (log(pow((prec[i,k,m]/(2*3.14159)),0.5) * exp(-0.5*(pow((y[i,k,m]-theta[i,k,m]),2)*prec[i,k,m])))) # deviance for normal likelihood

}

resarmdev[i,k] <- sum(resdev[i,k,1:fups[i]])
armdev[i,k] <- sum(dev[i,k,1:fups[i]])

}

resstudydev[i] <- sum(resarmdev[i, 1:narm[i]])
studydev[i] <- sum(armdev[i, 1:narm[i]])

for(k in 2:narm[i]){ # Treatment effects
}
}

for (k in 2:NT){ # Priors on relative treatment effects
}

for (k in 1:NT){ # Priors on absolute treatment effects
}

for (k in 2:Nclass){ # Priors on relative class effects
}

for (k in 1:Nclass){ # Priors on absolute class effects
}

for (c in 1:(NT-1)) {
for (k in (c+1):NT) { # UME priors
}
}

totresdev <- sum(resstudydev[])
totdev <- sum(studydev[])

# Model ends
}
"

  return(model)
}



#' Writes insert points for RegEx in MBNMA JAGS code
#'
#' @return A list with named elements containing character strings that match
#'   points in MBNMA JAGS code. These points can therefore be used to insert
#'   other lines of JAGS code into the correct section within the overall model
#'   code.
#'
write.inserts <- function() {
  insert.start <- "(+.# Begin Model Code\n)(+.)"
  insert.study <- "(+.# Run through all NS trials\n)(+.)"
  insert.arm <- "(+.# Run through all arms within a study\n)(+.)"
  insert.obs <- "(+.# Run through all observations within a study\n)(+.)"
  insert.te <- "(+.# Treatment effects\n)(+.)"
  insert.te.priors <- "(+.# Priors on relative treatment effects\n)(+.)"
  insert.end <- "(.+)(\n# Model ends)"
  insert.class.priors <- "(+.# Priors on relative class effects\n)(+.)"
  insert.absclass.priors <- "(+.# Priors on absolute class effects\n)(+.)"
  insert.ume.priors <- "(+.# UME priors\n)(+.)"
  insert.te.abs <- "(+.# Priors on absolute treatment effects\n)(+.)"

  return(inserts <- list("insert.start"=insert.start,
                         "insert.study"=insert.study,
                         "insert.arm"=insert.arm,
                         "insert.obs"=insert.obs,
                         "insert.te"=insert.te,
                         "insert.te.priors"=insert.te.priors,
                         "insert.end"=insert.end,
                         "insert.class.priors"=insert.class.priors,
                         "insert.absclass.priors"=insert.absclass.priors,
                         "insert.ume.priors"=insert.ume.priors,
                         "insert.te.abs"=insert.te.abs
  ))
}



#' Adds sections of JAGS code for an MBNMA model that correspond to alpha
#' parameters
#'
#' @param timecourse A character object that contains JAGS code for the
#'   time-course component of the model
#' @param model A character string representing the MBNMA model in JAGS code
#' @inheritParams mb.run
#'
#' @return A list of named elements: `model` is a character object of JAGS MBNMA
#'   model code that includes alpha parameter components of the model
#'   `timecourse` is a character object that contains JAGS code for the
#'   time-course component of the model, for which alpha will be indexed
#'   correctly
#'
write.alpha <- function(model, timecourse,
                        intercept, positive.scale) {

alpha.study.prior <- "
alpha[i] ~ dnorm(0,0.0001)
"
alpha.arm.prior <- "
alpha[i,k] ~ dnorm(0,0.0001)
"

  inserts <- write.inserts()

  if (intercept==FALSE) {
    timecourse <- gsub("(alpha\\[i,k\\] \\+)(+.)", "\\2", timecourse)	 # Removes alpha from timecourse
    timecourse <- gsub("(alpha\\[i\\] \\+)(+.)", "\\2", timecourse)	 # Removes alpha from timecourse
  }

  if (positive.scale==TRUE & grepl("alpha\\[i(,k)?\\]", timecourse)==TRUE) {
    timecourse <- gsub("(alpha\\[i(,k)?\\])", paste0("exp(", "\\1", ")"), timecourse)
  }

  if (grepl("alpha\\[i\\]", timecourse)==TRUE) {
    model <- gsub(inserts[["insert.study"]], paste0("\\1", alpha.study.prior, "\\2"), model)
  } else if (grepl("alpha\\[i,k\\]", timecourse)==TRUE) {
    model <- gsub(inserts[["insert.arm"]], paste0("\\1", alpha.arm.prior, "\\2"), model)
  }

  return(list("model"=model, "timecourse"=timecourse))
}





#' Adds sections of JAGS code for an MBNMA model that correspond to the
#' likelihood
#'
#' @inheritParams write.beta
#' @inheritParams write.alpha
#' @inheritParams mb.run
#'
#' @return A character object of JAGS MBNMA model code that includes likelihood
#'   components of the model
#'
write.likelihood <- function(model, timecourse, rho=NULL, covar=NULL) {

# Linear predictor
predictor <- "theta[i,k,m] <- "

norm.like <- "
y[i,k,m] ~ dnorm(theta[i,k,m], prec[i,k,m])
prec[i,k,m] <- pow(se[i,k,m], -2)
"

mnorm.like <- "
y[i,k,1:fups[i]] ~ dmnorm.vcov(theta[i,k,1:fups[i]], cov.mat[i,k,1:fups[i],1:fups[i]])
"

if (!is.null(rho)) {
  if (rho=="estimate") {
    rho.prior <- paste0("\n", "rho ~ dunif(-1,1)", "\n")
  } else if (is.numeric(rho)) {
    rho.prior <- paste0("\n", "rho <- ", rho, "\n")
  }
}

# if (is.character(rho)) {
#   rho.prior <- paste0("\n", "rho ~ ", rho, "\n")
# } else if (is.numeric(rho)) {
#   rho.prior <- paste0("\n", "rho <- ", rho, "\n")
# }

covar.cs <- "
# Generates covariance matrix upper and lower triangles
for (c in 1:(fups[i]-1)) {
for (r in (c+1):fups[i]) {
cov.mat[i,k,r,c] <- se[i,k,c]*se[i,k,r]*rho   # Lower triangle
cov.mat[i,k,c,r] <- se[i,k,c]*se[i,k,r]*rho   # Upper triangle
}
}

# Generates covariance matrix diagonals
for (m in 1:fups[i]) {
cov.mat[i,k,m,m] <- pow(se[i,k,m],2)
}
"

covar.ar1 <- "
# Generates covariance matrix upper and lower triangles
		for (c in 1:(fups[i]-1)) {
for (r in (c+1):fups[i]) {
cov.mat[i,k,r,c] <- se[i,k,c]*se[i,k,r]*cor[i,r,c]   # Lower triangle
cov.mat[i,k,c,r] <- se[i,k,c]*se[i,k,r]*cor[i,c,r]   # Upper triangle
}
}

# Generates covariance matrix diagonals
for (m in 1:fups[i]) {
cov.mat[i,k,m,m] <- pow(se[i,k,m],2)
}
"

timepoint.corr <- "
# Generates separate correlation coefficients for each two time points
	for (c in 1:(fups[i]-1)) {
for (r in (c+1):fups[i]) {
cor[i,r,c] <- pow(rho, (abs(time[i,r] - time[i,c])) / timedif.0[i])
cor[i,c,r] <- pow(rho, (abs(time[i,r] - time[i,c])) / timedif.0[i])
}
}
"

  inserts <- write.inserts()

  if (is.null(rho)) {
    model <- gsub(inserts[["insert.obs"]], paste0("\\1", norm.like, "\\2"), model)
  } else if (!is.null(rho)) {
    model <- gsub(inserts[["insert.arm"]], paste0("\\1", mnorm.like, "\\2"), model)
    model <- gsub(inserts[["insert.end"]], paste0("\\1", rho.prior, "\\2"), model)
    model <- gsub("(?<=\n){1}(.+?)(dev)(.+?)(?=\n){1}", "", model, perl=T) # Remove devs
    # if (rho=="estimate") {
    #   model <- gsub("(?<=\n){1}(.+?)(dev)(.+?)(?=\n){1}", "", model, perl=T) # Remove devs
    # }
    if (covar=="CS") {
      model <- gsub(inserts[["insert.arm"]], paste0("\\1", covar.cs, "\\2"), model)
    } else if (covar=="AR1") {
      model <- gsub(inserts[["insert.arm"]], paste0("\\1", covar.ar1, "\\2"), model)
      model <- gsub(inserts[["insert.study"]], paste0("\\1", timepoint.corr, "\\2"), model)
    }
  }

  # Add linear predictor
  model <- gsub(inserts[["insert.obs"]], paste0("\\1", paste0(predictor, timecourse), "\\2"), model)

  return(model)
}





#' Adds additional sections of JAGS code for an MBNMA model required for models
#' with a fractional polynomial time-course
#'
#' @inheritParams write.beta
#' @inheritParams write.alpha
#'
#' @details FUNCTION IS DEPRECATED
#'
#' @return A character object of JAGS MBNMA model code that includes fractional
#'   polynomial components of the model
#'
write.fract.poly <- function(model, timecourse) {

  if (grepl("time\\.fp\\.1", model)==TRUE) {

    # Fractional polynomial segments
    fract.poly.1 <- "
time.fp.1[i,1] <- 0
for (m in 2:fups[i]) {
time.fp.1[i,m] <- (equals(beta.3,0)*log(time[i,m]) + (1-equals(beta.3,0))*pow(time[i,m],beta.3))}
"
    fract.poly.2 <- "
time.fp.1[i,1] <- 0
time.fp.2[i,1] <- 0
for (m in 2:fups[i]) {
time.fp.1[i,m] <- (equals(beta.3,0)*log(time[i,m]) + (1-equals(beta.3,0))*pow(time[i,m],beta.3))
time.fp.2[i,m] <- ((1-equals(beta.4,beta.3))*(equals(beta.4,0)*log(time[i,m]) + (1-equals(beta.4,0))*pow(time[i,m],beta.4)  ) + equals(beta.4, beta.4)*(equals(beta.4,0)*log(time[i,m])*log(time[i,m]) + (1-equals(beta.4,0))*pow(time[i,m],beta.4) *log(time[i,m])))}
"

    inserts <- write.inserts()

    # Remove existing segments in order to replace
    model <- gsub(paste0("(.+)", "(\\ntime\\.fp\\.1\\[i\\,1\\]", ".+\\)\\})", "(.+)"), "\\1\\3", model)

    # # Replace beta.i with s.beta.i if absolute effects are used
    # for (i in 1:4) {
    #   if (grepl(paste0("beta.",i,"\\[treat"), model)) {
    #     beta.str <- paste0("beta.",i,"[k]")
    #   } else {beta.str <- paste0("beta.",i)}
    #
    #   fract.poly.1 <- gsub(paste0("(.+)(beta.",i,")(.+)"), paste0("\\1", beta.str, "\\3"), fract.poly.1)
    #   fract.poly.2 <- gsub(paste0("(.+)(beta.",i,")(.+)"), paste0("\\1", beta.str, "\\3"), fract.poly.2)
    # }
    #
    # # Add additional loop over arms if fract.poly contains s.beta
    # fract.poly.1 <- gsub("(.+\\{)(.+)(\\})", paste0("\\1", "\nfor (k in 1:NT){", "\\2", "}\n", "\\3"), fract.poly.1)
    # fract.poly.2 <- gsub("(.+\\{)(.+)(\\})", paste0("\\1", "\nfor (k in 1:NT){", "\\2", "}\n", "\\3"), fract.poly.2)

    # Add time powers for fractional polynomials (Jansen 2015)
    if (grepl("time\\.fp\\.1", timecourse)==TRUE & grepl("time\\.fp\\.2", timecourse)==FALSE) {
      model <- gsub(inserts[["insert.study"]], paste0("\\1", fract.poly.1, "\\2"), model)
    } else if (grepl("time\\.fp\\.2", timecourse)==TRUE) {
      model <- gsub(inserts[["insert.study"]], paste0("\\1", fract.poly.2, "\\2"), model)
    }

  }

  return(model)
}



#' Adds additional sections of JAGS code for an MBNMA model required for models
#' with a piecewise linear or fractional polynomial time-course
#'
#' @inheritParams write.beta
#' @inheritParams mb.run
#'
#' @return A character object of JAGS MBNMA model code that includes piecewise
#'   linear components of the model
write.piece.fract <- function(model, fun, beta.3) {

  if (fun=="piecelinear") {
    if (beta.3 %in% c("arm.common", "arm.random", "const.common", "const.random")) {
      model <- gsub("(.+)(beta\\.3 ~ dnorm\\(0,0\\.0001\\))(.+)",
                    paste0("\\1", "beta.3 ~ dunif(0, maxtime)", "\\3"),
                    model)

      model <- gsub("(.+)(beta\\.3\\[k\\] ~ dnorm\\(0,0\\.001\\))(.+)",
                    paste0("\\1", "beta.3[k] ~ dunif(0, maxtime)", "\\3"),
                    model)
    }
  }

  return(model)
}




#' Adds additional sections of JAGS code for an MBNMA model required for models
#' with a piecewise linear time-course
#'
#' @inheritParams write.beta
#' @inheritParams mb.run
#'
#' @details FUNCTION DEPRECATED
#'
#' @return A character object of JAGS MBNMA model code that includes piecewise
#'   linear components of the model
#'
write.piecelinear <- function(model, beta.1, beta.3) {
  if (grepl("alpha\\.knot", model)) {

    inserts <- write.inserts()

    if (grepl("alpha\\[i\\]", model)) {
      alpha.str <- "alpha[i]"
    } else if (grepl("alpha\\[i,k\\]", model)) {
      alpha.str <- "alpha[i,k]"
    }

    # Ensure alpha is on same scale (e.g. exponential) in knot as in overall model
    if (grepl("exp\\(alpha", model)) {
      alpha.str <- paste0("exp(", alpha.str, ")")
    }

    alpha.knot <- paste0("\nalpha.knot[i,k] <- ", alpha.str, " + (beta.1*beta.3)\n" )

    for (i in 1:4) {
      if (grepl(paste0("s.beta.",i), model)) {
        beta.str <- paste0("s.beta.",i,"[i,k]")
      } else {beta.str <- paste0("beta.",i)}

      alpha.knot <- gsub(paste0("(.+)(beta.",i,")(.+)"), paste0("\\1", beta.str, "\\3"), alpha.knot)
    }



#     if (grepl("alpha\\[i\\]", model)) {
# alpha.knot <- "
# alpha.knot[i,k] <- alpha[i] + (beta.1*beta.3)
# "
#     } else if (grepl("alpha\\[i,k\\]", model)) {
# alpha.knot <- "
# alpha.knot[i,k] <- alpha[i,k] + (beta.1*beta.3)
# "
#     }

    if (beta.1 %in% c("rel.common", "rel.random")) {
      alpha.knot <- gsub("(.+beta\\.1)(+.)", paste0("\\1", "[i,k]", "\\2"), alpha.knot)
    }

    model <- gsub("alpha.knot", "alpha.knot[i,k]", model)

    model <- gsub(inserts[["insert.arm"]], paste0("\\1", alpha.knot, "\\2"), model)

    if (beta.3 %in% c("abs.i", "arm.common", "arm.random", "const.common", "const.random")) {
      model <- gsub("(.+)(beta\\.3 ~ dnorm\\(0,0\\.0001\\))(.+)",
                    paste0("\\1", "beta.3 ~ dunif(0, maxtime)", "\\3"),
                    model)

      model <- gsub("(.+)(beta\\.3\\[k\\] ~ dnorm\\(0,0\\.001\\))(.+)",
                    paste0("\\1", "beta.3[k] ~ dunif(0, maxtime)", "\\3"),
                    model)
    }
  }
  return(model)
}




#' Generate objects required for write.beta and write.beta.ref
#' @noRd
write.beta.vars <- function(beta.1, beta.2, beta.3, beta.4) {


  for (i in 1:4) {
    # Single parameters
    if (!is.null(get(paste("beta", i, sep=".")))) {
      if (get(paste("beta", i, sep=".")) == "abs.i") {
        assign(paste("beta", i, "prior", sep="."),
               paste0("\n", "beta.", i, "[i,k] ~ dnorm(0,0.0001)", "\n")
        )

        assign(paste("beta", i, "abs", sep="."),
               paste0("\n", "s.beta.", i, "[i,k] <- beta.", i, "[i,k]", "\n")
        )
      } else if (get(paste("beta", i, sep=".")) == "arm.common") {
        assign(paste("beta", i, "prior", sep="."),
               paste0("\n", "beta.", i, "[k] ~ dnorm(0,0.001)", "\n")
        )

        assign(paste("beta", i, "abs", sep="."),
               paste0("\n", "s.beta.", i, "[i,k] <- beta.", i, "[treat[i,k]]", "\n")
        )
      } else if (get(paste("beta", i, sep=".")) == "arm.random") {
        assign(paste("beta", i, "prior", sep="."),
               paste0("\n", "beta.", i, "[k] ~ dnorm(0,0.001)", "\n")
        )

        assign(paste("beta", i, "abs", sep="."),
               paste0("\n", "s.beta.", i, "[i,k] ~ dnorm(beta.", i, "[treat[i,k]], prec.beta.", i, ")\n")
        )

        # Changed whilst fixing summary(MBNMA)
        assign(paste("beta", i, "abs.sd", sep="."),
               paste0("\n", "prec.beta.", i, " <- pow(sd.beta.", i, ", -2)", "\n", "sd.beta.", i, " ~ dnorm(0,0.0025) T(0,)", "\n")
        )
      } else if (get(paste("beta", i, sep=".")) %in% c("const.common", "const.random")) {
        assign(paste("beta", i, "prior", sep="."),
               paste0("\n", "beta.", i, " ~ dnorm(0,0.0001)", "\n")
        )

        if (get(paste("beta", i, sep=".")) == "const.common") {
          assign(paste("beta", i, "abs", sep="."),
                 paste0("\n", "s.beta.", i, "[i,k] <- beta.", i, "\n")
          )
        } else if (get(paste("beta", i, sep=".")) == "const.random") {
          assign(paste("beta", i, "abs", sep="."),
                 paste0("\n", "s.beta.", i, "[i,k] ~ dnorm(beta.", i, ", prec.beta.", i, ")", "\n")
          )
          assign(paste("beta", i, "abs.sd", sep="."),
                 paste0("\n", "prec.beta.", i, " <- pow(sd.beta.", i, ", -2)", "\n", "sd.beta.", i, " ~ dnorm(0,0.0025) T(0,)", "\n")
          )
        }
      } else if (is.numeric(get(paste("beta", i, sep=".")))) {
        assign(paste("beta", i, "prior", sep="."),
               paste0("\n", "beta.", i, " <- ", get(paste("beta", i, sep=".")), "\n")
        )
      }
    }


    # Split beta parameters
    assign(paste("beta", i, "rel", sep="."),
           paste0("\n", "beta.", i, "[i,k] <- mu.", i, "[i] + delta.", i, "[i,k]", "\n")
    )

    # Reference treatment prior
    assign(paste("mu", i, "prior", sep="."),
           paste0("\n", "mu.", i, "[i] ~ dnorm(0,0.0001)", "\n")
    )

    # Reference relative effects set to zero
    assign(paste("delta", i, "ref", sep="."),
           paste0("\n", "delta.", i, "[i,1] <- 0", "\n")
    )

    # Fixed treatment effects
    assign(paste("delta", i, "fe", sep="."),
           paste0("\n", "delta.", i, "[i,k] <- d.", i, "[treat[i,k]] - d.", i, "[treat[i,1]]", "\n")
    )

    # Random treatment effects with multi-arm correction
    assign(paste("delta", i, "re", sep="."),
           paste0("\n", "delta.", i, "[i,k] ~ dnorm(md.", i, "[i,k], taud.", i, "[i,k])", "\n",
                  "md.", i, "[i,k] <- d.", i, "[treat[i,k]] - d.", i, "[treat[i,1]] + sw.", i, "[i,k]", "\n",
                  "taud.", i, "[i,k] <- tau.", i, " * 2*(k-1)/k", "\n",
                  "w.", i, "[i,k] <- (delta.", i, "[i,k] - d.", i, "[treat[i,k]] + d.", i, "[treat[i,1]])", "\n",
                  "sw.", i, "[i,k] <- sum(w.", i, "[i,1:(k-1)])/(k-1)", "\n"
           )
    )

    # Multi-arm correction for random treatment effects
    assign(paste("multiarm", i, sep="."),
           paste0("\n", "w.", i, "[i,1] <- 0", "\n")
    )

    # Fixed treatment effects - UME
    assign(paste("delta", i, "fe", "ume", sep="."),
           paste0("\n", "delta.", i, "[i,k] <- d.", i, "[treat[i,1],treat[i,k]]", "\n")
    )

    # Random treatment effects  - UME
    assign(paste("delta", i, "re", "ume", sep="."),
           paste0("\n", "delta.", i, "[i,k] ~ dnorm(d.", i, "[treat[i,1],treat[i,k]], tau.", i, ")", "\n")
    )

    # Placebo relative effect
    assign(paste("d.zero", i, sep="."),
           paste0("\n", "d.", i, "[1] <- 0", "\n")
    )

    # Prior active treatment relative effect
    assign(paste("d.prior", i, sep="."),
           paste0("\n", "d.", i, "[k] ~ dnorm(0,0.001)", "\n")
    )

    # Placebo relative class effect
    assign(paste("D.zero", i, sep="."),
           paste0("\n", "D.", i, "[1] <- 0", "\n")
    )

    # UME prior active treatment relative effect
    assign(paste("ume.prior", i, sep="."),
           paste0("\n", "d.", i, "[c,k] ~ dnorm(0,0.001)", "\n")
    )

    # Placebo relative effect for UME
    assign(paste("ume.zero", i, sep="."),
           paste0("\n", "d.", i, "[1,1] <- 0", "\n")
    )

    # Class effect on active treatment relative effect - common
    assign(paste("d.class.fe", i, sep="."),
           paste0("\n", "d.", i, "[k]  <- D.", i, "[class[k]]", "\n")
    )
    assign(paste("d.class.re", i, sep="."),
           paste0("\n", "d.", i, "[k]  ~ dnorm(D.", i, "[class[k]], tau.D.", i, ")\n")
    )

    # Prior on active treatment class effect
    assign(paste("D.prior", i, sep="."),
           paste0("\n", "D.", i, "[k] ~ dnorm(0,0.001)", "\n")
    )

    # SD prior
    assign(paste("sd.prior", i, sep="."),
           paste0("\n", "sd.", i, " ~ dnorm(0,0.0025) T(0,)", "\n",
                  "tau.", i, " <- pow(sd.", i, ", -2)", "\n"
           )
    )

    # SD prior D class effects
    assign(paste("sd.D.prior", i, sep="."),
           paste0("\n", "sd.D.", i, " ~ dnorm(0,0.0025) T(0,)", "\n",
                  "tau.D.", i, " <- pow(sd.D.", i, ", -2)", "\n"
           )
    )



    #### BETA effects ####

    # Class effect on active treatment relative effect - common
    assign(paste("beta.class.fe", i, sep="."),
           paste0("\n", "beta.", i, "[k]  <- BETA.", i, "[class[k]]", "\n")
    )
    assign(paste("beta.class.re", i, sep="."),
           paste0("\n", "beta.", i, "[k]  ~ dnorm(BETA.", i, "[class[k]], tau.BETA.", i, ")\n")
    )

    # Prior on active treatment class effect
    assign(paste("BETA.prior", i, sep="."),
           paste0("\n", "BETA.", i, "[k] ~ dnorm(0,0.001)", "\n")
    )

    # SD prior BETA class effects
    assign(paste("sd.BETA.prior", i, sep="."),
           paste0("\n", "sd.BETA.", i, " ~ dnorm(0,0.0025) T(0,)", "\n",
                  "tau.BETA.", i, " <- pow(sd.BETA.", i, ", -2)", "\n"
           )
    )


    #### For write.beta.ref ####

    # Beta parameters for reference treatment only
    assign(paste("beta", i, "ref", sep="."),
           paste0("\n", "beta.", i, "[i,k] <- mu.", i, "[i]", "\n")
    )

    # Reference treatment fixed
    assign(paste("mu", i, "fixed", sep="."),
           paste0("\n", "mu.", i, "[i] <- m.mu.", i, "\n")
    )

    # Reference treatment random
    assign(paste("mu", i, "random", sep="."),
           paste0("\n", "mu.", i, "[i] ~ dnorm(m.mu.", i, ", tau.mu.", i, ")",  "\n")
    )

    # Reference treatment common prior
    assign(paste("m.mu", i, sep="."),
           paste0("\n", "m.mu.", i, " ~ dnorm(0,0.0001)",  "\n")
    )

    # Reference treatment sd prior
    assign(paste("sd.mu", i, sep="."),
           paste0("\n", "sd.mu.", i, " ~ dnorm(0,0.0025) T(0,)",  "\n",
                  "tau.mu.", i, " <- pow(sd.mu.", i, ", -2)", "\n"
           )
    )
  }
  varnames <- ls()

  vars <- list()
  for (i in seq_along(varnames)) {
    vars[[varnames[i]]] <- get(varnames[i])
  }
  return(vars)
}







#' Adds sections of JAGS code for an MBNMA model that correspond to beta
#' parameters
#'
#' @param model A character object of JAGS MBNMA model code
#' @param timecourse A character object representing the time-course used in the MBNMA model
#' @inheritParams mb.run
#'
#' @return A character object of JAGS MBNMA model code that includes beta
#'   parameter components of the model
#'
write.beta <- function(model, timecourse,
                       beta.1, beta.2=NULL, beta.3=NULL, beta.4=NULL,
                       UME, class.effect
                       ) {

  inserts <- write.inserts()

  # Assign treatment effect segments
  vars <- write.beta.vars(beta.1, beta.2, beta.3, beta.4)

  # Add to model
  for (i in 1:4) {
    # If a beta parameter has relative effects for indirect evidence calculation
    if (!is.null(get(paste0("beta.", i)))) {

      # RELATIVE TIME-COURSE PARAMETERS
      if ((get(paste0("beta.", i))=="rel.common" | get(paste0("beta.", i))=="rel.random") &
          grepl(paste0("beta.", i , "\\[i,k\\]"), timecourse)==TRUE) {

        model <- gsub(inserts[["insert.arm"]], paste0("\\1", vars[[paste("beta", i, "rel", sep=".")]], "\\2"), model)
        model <- gsub(inserts[["insert.study"]], paste0("\\1", vars[[paste("mu", i, "prior", sep=".")]], "\\2"), model)
        model <- gsub(inserts[["insert.study"]], paste0("\\1", vars[[paste("delta", i, "ref", sep=".")]], "\\2"), model)

        # For treatment effects across network
        if (!((paste0("beta.", i) %in% UME))) {
          model <- gsub(inserts[["insert.start"]], paste0("\\1", vars[[paste("d.zero", i, sep=".")]], "\\2"), model)
        }

        if (!  ((paste0("beta.", i) %in% names(class.effect)) | ((paste0("beta.", i) %in% UME)))) {
          model <- gsub(inserts[["insert.te.priors"]], paste0("\\1", vars[[paste("d.prior", i, sep=".")]], "\\2"), model)

        # CLASS EFFECTS
        } else if ((paste0("beta.", i) %in% names(class.effect))) {

          if (class.effect[[paste0("beta.", i)]] == "common") {
            model <- gsub(inserts[["insert.te.priors"]], paste0("\\1", vars[[paste("d.class.fe", i, sep=".")]], "\\2"), model)
            model <- gsub(inserts[["insert.start"]], paste0("\\1", vars[[paste("D.zero", i, sep=".")]], "\\2"), model)
            model <- gsub(inserts[["insert.class.priors"]], paste0("\\1", vars[[paste("D.prior", i, sep=".")]], "\\2"), model)
          } else if (class.effect[[paste0("beta.", i)]] == "random") {
            model <- gsub(inserts[["insert.te.priors"]], paste0("\\1", vars[[paste("d.class.re", i, sep=".")]], "\\2"), model)
            model <- gsub(inserts[["insert.start"]], paste0("\\1", vars[[paste("D.zero", i, sep=".")]], "\\2"), model)
            model <- gsub(inserts[["insert.class.priors"]], paste0("\\1", vars[[paste("D.prior", i, sep=".")]], "\\2"), model)
            model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("sd.D.prior", i, sep=".")]], "\\2"), model)
          }
        }

        # UME
        if (get(paste0("beta.", i))=="rel.common") {
          if ((paste0("beta.", i) %in% UME)) {
            model <- gsub(inserts[["insert.te"]], paste0("\\1", vars[[paste("delta", i, "fe", "ume", sep=".")]], "\\2"), model)
            model <- gsub(inserts[["insert.ume.priors"]], paste0("\\1", vars[[paste("ume", "prior", i, sep=".")]], "\\2"), model)
            model <- gsub(inserts[["insert.start"]], paste0("\\1", vars[[paste("ume", "zero", i, sep=".")]], "\\2"), model)
          } else {
            model <- gsub(inserts[["insert.te"]], paste0("\\1", vars[[paste("delta", i, "fe", sep=".")]], "\\2"), model)
          }
        } else if (get(paste0("beta.", i))=="rel.random") {
          model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("sd.prior", i, sep=".")]], "\\2"), model)
          if ((paste0("beta.", i) %in% UME)) {
            model <- gsub(inserts[["insert.te"]], paste0("\\1", vars[[paste("delta", i, "re", "ume", sep=".")]], "\\2"), model)
            model <- gsub(inserts[["insert.ume.priors"]], paste0("\\1", vars[[paste("ume", "prior", i, sep=".")]], "\\2"), model)
          } else {
            model <- gsub(inserts[["insert.te"]], paste0("\\1", vars[[paste("delta", i, "re", sep=".")]], "\\2"), model)
            model <- gsub(inserts[["insert.study"]], paste0("\\1", vars[[paste("multiarm", i, sep=".")]], "\\2"), model)
          }
        }
      }

      # ABSOLUTE TIME-COURSE PARAMETERS
      if (grepl(paste0("beta.", i), model)==TRUE) {
        # If beta is estimmated as an absolute value
        if (get(paste0("beta.", i)) %in% c("arm.common", "arm.random", "const.common", "const.random")) {
          model <- gsub(paste0("(+.)", paste0("(beta\\.", i, ")")),
                        paste0("\\1", "s.", "\\2", "[i,k]"),
                        model)
        }

        if (get(paste0("beta.", i)) %in% c("arm.common", "arm.random")) {
          model <- gsub(inserts[["insert.arm"]], paste0("\\1", vars[[paste("beta", i, "abs", sep=".")]], "\\2"), model)

          if (!(paste0("beta.", i) %in% names(class.effect))) {
            model <- gsub(inserts[["insert.te.abs"]], paste0("\\1", vars[[paste("beta", i, "prior", sep=".")]], "\\2"), model)
          }

          if (get(paste0("beta.", i)) == "arm.random") {
            model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("beta", i, "abs.sd", sep=".")]], "\\2"), model)
          }

        } else if (get(paste0("beta.", i)) %in% c("const.common", "const.random")) {
          model <- gsub(inserts[["insert.arm"]], paste0("\\1", vars[[paste("beta", i, "abs", sep=".")]], "\\2"), model)
          model <- gsub(inserts[["insert.start"]], paste0("\\1", vars[[paste("beta", i, "prior", sep=".")]], "\\2"), model)

          if (get(paste0("beta.", i)) == "const.random") {
            model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("beta", i, "abs.sd", sep=".")]], "\\2"), model)
          }
        }

        # If beta is given a numeric value
        if (is.numeric(get(paste0("beta.", i)))) {
          model <- gsub(inserts[["insert.start"]], paste0("\\1", vars[[paste("beta", i, "prior", sep=".")]], "\\2"), model)
        }

        # CLASS EFFECTS
        if (get(paste0("beta.",i)) %in% c("arm.common", "arm.random")) {
          if ((paste0("beta.", i) %in% names(class.effect))) {
            if (class.effect[[paste0("beta.", i)]] == "common") {
              model <- gsub(inserts[["insert.te.abs"]], paste0("\\1", vars[[paste("beta.class.fe", i, sep=".")]], "\\2"), model)
              model <- gsub(inserts[["insert.absclass.priors"]], paste0("\\1", vars[[paste("BETA.prior", i, sep=".")]], "\\2"), model)
            } else if (class.effect[[paste0("beta.", i)]] == "random") {
              model <- gsub(inserts[["insert.te.abs"]], paste0("\\1", vars[[paste("beta.class.re", i, sep=".")]], "\\2"), model)
              model <- gsub(inserts[["insert.absclass.priors"]], paste0("\\1", vars[[paste("BETA.prior", i, sep=".")]], "\\2"), model)
              model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("sd.BETA.prior", i, sep=".")]], "\\2"), model)
            }
          }
        }

      }

    }
  }

  return(model)
}







#' Adds sections of JAGS code for an MBNMA reference synthesis model that
#' correspond to beta parameters
#' @noRd
write.beta.ref <- function(model, timecourse,
                       beta.1=beta.1, beta.2=beta.2, beta.3=beta.3, beta.4=beta.4,
                       mu.synth="random"
) {

  inserts <- write.inserts()

  # # Switch arm beta for const beta (the two are equivalent in )
  # for (i in 1:4) {
  #   if (!is.null(get(paste0("beta.", i)))) {
  #     if (get(paste0("beta.", i)) == "arm.common") {
  #       assign(paste0("beta.", i), "const.common")
  #     } else if (get(paste0("beta.", i)) == "arm.random") {
  #       assign(paste0("beta.", i), "const.random")
  #     }
  #   }
  # }

  # Assign treatment effect segments
  vars <- write.beta.vars(beta.1, beta.2, beta.3, beta.4)

  # Add to model
  for (i in 1:4) {
    # If a beta parameter has relative effects for indirect evidence calculation
    if (!is.null(get(paste0("beta.", i)))) {
      if ((get(paste0("beta.", i))=="rel.common" | get(paste0("beta.", i))=="rel.random") &
          grepl(paste0("beta.", i , "\\[i,k\\]"), timecourse)==TRUE) {

        model <- gsub(inserts[["insert.arm"]], paste0("\\1", vars[[paste("beta", i, "ref", sep=".")]], "\\2"), model)

        if (mu.synth=="fixed") {
          model <- gsub(inserts[["insert.study"]], paste0("\\1", vars[[paste("mu", i, "fixed", sep=".")]], "\\2"), model)
          model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("m.mu", i, sep=".")]], "\\2"), model)
        } else if (mu.synth=="random") {
          model <- gsub(inserts[["insert.study"]], paste0("\\1", vars[[paste("mu", i, "random", sep=".")]], "\\2"), model)
          model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("m.mu", i, sep=".")]], "\\2"), model)
          model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("sd.mu", i, sep=".")]], "\\2"), model)
        }
      }


      # ABSOLUTE TIME-COURSE PARAMETERS
      if (grepl(paste0("beta.", i), model)==TRUE) {
        # If beta is estimmated as an absolute value
        if (get(paste0("beta.", i)) %in% c("abs.i", "arm.common", "arm.random", "const.common", "const.random")) {
          model <- gsub(paste0("(+.)", paste0("(beta\\.", i, ")")),
                        paste0("\\1", "s.", "\\2", "[i,k]"),
                        model)
        }

        # if (get(paste0("beta.", i)) == "arm.common") {
        #   model <- gsub(inserts[["insert.te.abs"]], paste0("\\1", vars[[paste("beta", i, "prior", sep=".")]], "\\2"), model)
        #   model <- gsub(inserts[["insert.arm"]], paste0("\\1", vars[[paste("beta", i, "abs", sep=".")]], "\\2"), model)
        # } else if (get(paste0("beta.", i)) == "arm.random") {
        #   model <- gsub(inserts[["insert.te.abs"]], paste0("\\1", vars[[paste("beta", i, "prior", sep=".")]], "\\2"), model)
        #   model <- gsub(inserts[["insert.arm"]], paste0("\\1", vars[[paste("beta", i, "abs", sep=".")]], "\\2"), model)
        #   model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("beta", i, "abs.sd", sep=".")]], "\\2"), model)
        # } else if (get(paste0("beta.", i)) == "const.common") {

        if (get(paste0("beta.", i)) == "const.common") {
          model <- gsub(inserts[["insert.arm"]], paste0("\\1", vars[[paste("beta", i, "abs", sep=".")]], "\\2"), model)
          model <- gsub(inserts[["insert.start"]], paste0("\\1", vars[[paste("beta", i, "prior", sep=".")]], "\\2"), model)
        } else if (get(paste0("beta.", i)) == "const.random") {
          model <- gsub(inserts[["insert.arm"]], paste0("\\1", vars[[paste("beta", i, "abs", sep=".")]], "\\2"), model)
          model <- gsub(inserts[["insert.start"]], paste0("\\1", vars[[paste("beta", i, "prior", sep=".")]], "\\2"), model)
          model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("beta", i, "abs.sd", sep=".")]], "\\2"), model)
        }

        # If beta is given a numeric value
        if (is.numeric(get(paste0("beta.", i)))) {
          model <- gsub(inserts[["insert.start"]], paste0("\\1", vars[[paste("beta", i, "prior", sep=".")]], "\\2"), model)
        }
      }
    }
  }

  return(model)
}




#' Adds correlation between time-course relative effects
#'
#' This uses a Wishart prior as default for modelling the correlation
#'
#' @inheritParams mb.run
#' @inheritParams write.beta
write.cor <- function(model, var.scale=NULL, class.effect=list()) {

  if (length(class.effect)>0) {
    message("Class effects cannot be modelled with correlation between time-course relative effects - correlation will be ignored")
  } else {

    sufparams <- vector()
    for (i in 1:4) {
      if (grepl(paste0("\nd\\.", i, "\\[k\\] ~"), model)) {
        sufparams <- append(sufparams, i)
      }
    }

    # Check the number of parameters modelled using relative effects
    mat.size <- length(sufparams)
    if (mat.size>=2) {
      model <- write.cov.mat(model, sufparams=sufparams,
                             cor="estimate", cor.prior="wishart",
                             var.scale=var.scale)
    }
  }
  return(model)
}






#' Adds covariance matrix for correlation between relative effects
#'
#' @param sufparams A numeric vector of dose-response/time-course parameter suffixes. It
#'  should be the same length as the number of relative effects (i.e. the covariance
#'  matrix size).
#' @noRd
write.cov.mat <- function(model, sufparams, cor="estimate", cor.prior="wishart",
                          var.scale=NULL) {

  inserts <- write.inserts()

  jagswish <- "
for (r in 1:mat.size) {
d.prior[r] <- 0
}

inv.R ~ dwish(Omega[,], mat.size)

for (r in 1:(mat.size-1)) {  # Covariance matrix upper/lower triangles
for (c in (r+1):mat.size) {
Omega[r,c] <- 0   # Lower triangle
Omega[c,r] <- 0   # Upper triangle
}
}
"

  jagsrho <- "
for (r in 1:mat.size) {
d.prior[r] <- 0
R[r,r] <- 1000    # Covariance matrix diagonals
}

for (r in 1:(mat.size-1)) {  # Covariance matrix upper/lower triangles
for (c in (r+1):mat.size) {
R[r,c] <- 1000*rho[1]   # Lower triangle
R[c,r] <- 1000*rho[1]   # Upper triangle
}
}
"

  mat.size <- length(sufparams)
  for (i in seq_along(sufparams)) {
    # Change d.1[k] ~ dnorm(0,0.001)  to   d.1[k] <- d.mult[1,k]
    model <- gsub(paste0("d\\.", sufparams[i], "\\[k\\] ~ [a-z]+\\([0-9]+(\\.[0-9]+)?,[0-9]+(\\.?[0-9]+)?\\)\\\n"),
                  paste0("d.", sufparams[i], "[k] <- mult[", i, ",k]\n"),
                  model
    )
  }

  if (cor.prior=="wishart") {
    addcode <- jagswish
    model <- gsub(inserts[["insert.te.priors"]],
                  paste0("\\1mult[1:", mat.size, ",k] ~ dmnorm(d.prior[], inv.R[1:", mat.size, ", 1:", mat.size, "])\\2"),
                  model
    )

    # Check that var.scale has correct length and add omega to code
    if (is.null(var.scale)) {
      var.scale <- rep(1,mat.size)
    } else if (length(var.scale)!=mat.size) {
      stop(paste0("`var.scale` must be a numeric vector whose length is the size of the covariance matrix for dose-response parameters.\nCovariance matrix size = ", mat.size))
    }
    for (i in seq_along(var.scale)) {
      model <- gsub(inserts[["insert.end"]],
                    paste0("\\1Omega[", i, ",", i, "] <- ", var.scale[i], "\n\\2"),
                    model)
    }

  } else if (cor.prior=="rho") {
    addcode <- jagsrho
    model <- gsub(inserts[["insert.te.priors"]],
                  paste0("\\1mult[1:", mat.size, ",k] ~ dmnorm.vcov(d.prior[], R[1:", mat.size, ", 1:", mat.size, "])\\2"),
                  model
    )

    if (cor=="estimate") {
      addcode <- paste(addcode, "
                       for (m in 1:(mat.size-1)) {
                       rho[m] ~ dunif(-1,1)
                       }
                       ")
    } else if (is.numeric(cor)) {
      # Add values for rho assigned by user
      if (length(cor)!=(mat.size)-1) {
        stop("Length of numeric vector assigned to `cor` must equal the size of the correlation matrix - 1")
      }
      for (m in seq_along(cor)) {
        model <- gsub(inserts[["insert.end"]],
                      paste0("\\1rho[", m, "] <- ", cor[m]))
      }
  }
}

  addcode <- gsub("mat\\.size", mat.size, addcode)
  model <- gsub(inserts[["insert.end"]], paste0("\\1", addcode, "\\2"), model)
  return(model)
}






#' Removes any loops from MBNMA model JAGS code that do not contain any
#' expressions
#'
#' @inheritParams write.beta
#'
#' @return A character object of JAGS MBNMA model code that has had empty loops
#'   removed from it
#'
write.remove.loops <- function(model) {
  # Remove empty loops
  empty.loops <- list(
    "for \\(k in 2\\:NT\\)\\{ \\# Priors on relative treatment effects\\\n}", # ume.loop
    "for \\(k in 1\\:NT\\)\\{ \\# Priors on absolute treatment effects\\\n}", # absolute.te.loop
    "for\\(k in 2\\:narm\\[i\\]\\)\\{ \\# Treatment effects\\\n\\}", # treatment effects loop
    "for \\(k in 2\\:Nclass\\)\\{ \\# Priors on relative class effects\\\n\\}", # rel class loop
    "for \\(k in 1\\:Nclass\\)\\{ \\# Priors on absolute class effects\\\n\\}", # arm class loop
    "for \\(c in 1\\:\\(NT-1\\)\\) \\{\\\nfor \\(k in \\(c\\+1\\)\\:NT\\) \\{ \\# UME priors\\\n\\}\\\n\\}" # UME loop
  )

  for (i in seq_along(empty.loops)) {
    if (grepl(empty.loops[[i]], model)==TRUE) {
      model <- gsub(paste0("(.+)(", empty.loops[[i]], ")(.+)"),
                    paste0("\\1", "\\3"), model)
    }
  }

  return(model)
}





#' Write MBNMA time-course models JAGS code for synthesis of studies
#' investigating reference treatment
#'
#' Writes JAGS code for a Bayesian time-course model for model-based network
#' meta-analysis (MBNMA) that pools reference treatment effects from different
#' studies. This model only pools single study arms and therefore does not pool
#' relative effects.
#'
#' @param mu.synth A string that takes the value `fixed` or `random`, indicating
#'   the type of synthesis model to use
#' @inheritParams mb.run
#'
#' @return A character object of JAGS MBNMA model code that includes beta
#'   parameter components of the model
#'
#' @examples
#' # Write an exponential time-course MBNMA synthesis model
#' model <- write.ref.synth(fun="exponential",
#'   alpha="arm", beta.1="rel.common", mu.synth="fixed")
#' cat(model) # Concatenates model representations making code more easily readable
#'
#' @export
write.ref.synth <- function(fun="linear", user.fun=NULL, alpha="arm", beta.1="rel.common", beta.2=NULL, beta.3=NULL, beta.4=NULL,
                        positive.scale=TRUE, intercept=TRUE, rho=NULL, covar=NULL,
                        mu.synth="random",
                        class.effect=list(), UME=FALSE, priors=NULL) {

  ####### VECTORS #######

  parameters <- c("alpha", "beta.1", "beta.2", "beta.3", "beta.4")

  # Change UME to relate to parameters in model
  if (UME!=FALSE) {
    stop("UME models cannot be used for synthesis of reference treatment")
  } else if (length(class.effect)!=0) {
    warning("Synthesis of reference treatment effect is independent of class effects specified in the MBNMA model")
    class.effect <- list()
  }

  write.check(fun=fun, user.fun=user.fun, alpha=alpha,
              beta.1=beta.1, beta.2=beta.2, beta.3=beta.3, beta.4=beta.4,
              positive.scale=positive.scale, intercept=intercept,
              rho=rho, covar=covar,
              class.effect=class.effect, UME=UME)

  model <- write.model()

  suppressMessages(
    timecourse <- time.fun(fun=fun, user.fun=user.fun,
                           alpha=alpha, beta.1=beta.1, beta.2=beta.2,
                           beta.3=beta.3, beta.4=beta.4)[["jagscode"]]
  )

  alphacode <- write.alpha(model, timecourse,
                           intercept=intercept, positive.scale=positive.scale)
  timecourse <- alphacode[["timecourse"]]
  model <- alphacode[["model"]]

  model <- write.likelihood(model, timecourse, rho, covar)

  #model <- write.fract.poly(model, timecourse)

  model <- write.beta.ref(model, timecourse,
                          beta.1=beta.1, beta.2=beta.2, beta.3=beta.3, beta.4=beta.4,
                          mu.synth=mu.synth
                          )

  # Ammend fract.poly for s.beta changes
  #model <- write.fract.poly(model, timecourse)

  #model <- write.piecelinear(model, beta.1, beta.3)

  model <- write.piece.fract(model, fun, beta.3)

  model <- write.remove.loops(model)

  # Add user-specific priors (but only those for reference treatment)
  if (!is.null(priors)) {
    ref.priors <- get.prior(model)
    for (i in seq_along(ref.priors)) {
      if (names(ref.priors)[i] %in% names(priors)) {
        ref.priors[[i]] <- priors[[names(ref.priors)[i]]]
      }
    }
    model <- replace.prior(ref.priors, model=model)
  }

  return(model)
}








#' Get current priors from JAGS model code
#'
#' Identical to `get.prior()` in MBNMAdose.
#' This function takes JAGS model presented as a string and identifies what
#' prior values have been used for calculation.
#'
#' @inheritParams write.beta
#'
#' @return A character vector, each element of which is a line of JAGS code
#'   corresponding to a prior in the JAGS code.
#'
#' @details Even if an MBNMA model that has not initialised successfully and
#'   results have not been calculated, the JAGS model for it is saved in
#'   `MBNMA$model.arg$jagscode` and therefore priors can still be obtained.
#'   This allows for priors to be changed even in failing models, which may help
#'   solve issues with initialisation.
#'
#' @examples
#' \donttest{
#' # Create mb.network object using an MBNMAtime dataset
#' network <- mb.network(osteopain)
#'
#' # Create mb.network object using an MBNMAdose dataset
#'
#' # Run linear MBNMA
#' result <- mb.linear(network,
#'   slope=list(pool="rel", method="random"))
#'
#' # Obtain model prior values
#' get.prior(result$model.arg$jagscode)
#'
#' # ...also equivalent to
#' print(result$model.arg$priors)
#' }
#'
#' @export
get.prior <- function(model) {

  # Run Checks
  checkmate::assertCharacter(model, len=1)

  model <- strsplit(model, split="\n")[[1]]

  priorcode <- model[c(grep("^.+~ [A-z]+\\([-?0-9]", model),
                       grep("^.+~ [A-z]+\\(Omega", model))]

  priorlist <- strsplit(priorcode, split=" +?~ +?")
  priors <- list()
  for (i in seq_along(priorlist)) {
    priorname <- unlist(strsplit(priorlist[[i]][1], split="\\["))[1]
    priors[[priorname]] <- priorlist[[i]][2]
  }

  return(priors)
}



#' Replace original priors in an MBNMA model with new priors
#'
#' Identical to `replace.prior()` in MBNMAdose.
#'
#' This function takes new priors, as specified by the user, and adds them to
#' the JAGS code from an MBNMA model. New priors replace old priors in the JAGS
#' model.
#'
#' @inheritParams get.prior
#' @param mbnma An S3 object of class `c("mbnma", "rjags")` generated by running a
#'   time-course MBNMA model.
#' @param priors A named list of parameter values (without indices) and
#'   replacement prior distribution values given as strings
#'   **using distributions as specified in JAGS syntax**.
#'
#' @details Values in `priors` can include any JAGS functions/distributions
#'   (e.g. censoring/truncation).
#'
#' @return A character object of JAGS MBNMA model code that includes the new
#'   priors in place of original priors
#'
replace.prior <- function(priors, model=NULL, mbnma=NULL) {
  # priors is a named list of parameter values (without indices) and replacement
  #prior values given as strings USING DISTRIBUTIONS AS SPECIFIED IN JAGS SYNTAX (i.e.
  #dnorm() is specified using mean and precision rather than mean and SD.
  #It can include JAGS functions (e.g. censoring/truncation)
  #e.g. for a half-normal SD prior list("sd.et50"="dnorm(0,0.5) T(0,)")

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", null.ok=TRUE, add=argcheck)
  checkmate::assertCharacter(model, len=1, null.ok=TRUE, add=argcheck)
  checkmate::assertList(priors, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (!is.null(mbnma) & !is.null(model)) {
    stop("Must provide EITHER an existing MBNMA model (using `mbnma`) OR MBNMA JAGS code (using `model`)")
  }

  if (!is.null(mbnma)) {
    model <- strsplit(mbnma$model.arg$jagscode, split="\n")[[1]]
  } else if (!is.null(model)) {
    model <- strsplit(model, split="\n")[[1]]
  } else {
    stop("Must provide EITHER an existing MBNMA model (using `mbnma`) OR MBNMA JAGS code (using `model`)")
  }

  for (i in seq_along(priors)) {
    # Checks
    if (length(grep(paste0("^( +)?", names(priors)[i]), model))==0) {
      stop("Prior named ", names(priors)[i], " not found in the model code. Check priors currently present in model code using get.prior()")
    # } else if (length(grep(paste0("^( +)?", names(priors)[i]), model))>1) {
    #   stop("Prior named ", names(priors)[i], " has matched on multiple instances in the model code. Check priors currently present in model code using get.prior()")
    }

    #line <- grep(paste0("^( +)?", names(priors)[i]), model)
    line <- grep(paste0("^( +)?", names(priors)[i], ".+~"), model)
    state <- model[line]
    model[line] <- gsub("(^.+~ )(.+$)", paste0("\\1", priors[[i]]), state)
  }

  # Cut irrelevant section from JAGS code
  start <- grep("^model\\{", model)
  end <- grep("# Model ends", model) + 1

  model <- paste(model[start:end], collapse="\n")

  return(model)
}
