# Functions for writing MBNMA models in JAGS
# Author: Hugo Pedder
# Date created: 2020-01-04



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
#' # Write a linear time-course MBNMA:
#' # random treatment effects on beta.1
#' # equal baselines in study arms
#' model <- mb.write(fun=tpoly(degree=1, pool.1="rel", method.1="random"))
#'
#' # Write an emax time-course MBNMA with:
#' # a Hill parameter
#' # no intercept
#' model <- mb.write(fun=temax(pool.emax="rel", method.emax="common",
#'     pool.et50="abs", method.et50="common", pool.hill="abs", method.hill="common"),
#'   intercept=TRUE)
#'
#' # Write a log-linear time-course MBNMA with:
#' # AR1 correlation between time points
#' model <- mb.write(fun=tloglin(),
#'   rho="dunif(0,1)", covar="AR1")
#'
#' # Define a user-defined time-course relationship for the MBNMA JAGS model
#' userfun <- ~ (exp(beta.1 * time) / (beta.2 * time))
#' model <- mb.write(fun=tuser(fun=userfun,
#'     pool.1="rel", method.1="random",
#'     pool.2="rel", method.2="common"))
#'
#' @export
mb.write <- function(fun=tpoly(degree = 1), link="identity", positive.scale=TRUE, intercept=NULL,
                     rho=0, covar="varadj", omega=NULL, corparam=TRUE, sdscale=FALSE,
                     class.effect=list(), UME=FALSE) {


  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(fun, classes = "timefun", add=argcheck)
  checkmate::assertChoice(link, choices=c("identity", "log", "smd"), add=argcheck)
  checkmate::assertLogical(positive.scale, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertLogical(intercept, len=1, null.ok=TRUE, any.missing=FALSE, add=argcheck)
  checkmate::assertChoice(covar, choices=c("varadj", "CS", "AR1"), null.ok=FALSE, add=argcheck)
  checkmate::assertList(class.effect, unique=FALSE, add=argcheck)
  checkmate::assertLogical(corparam, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)


  ####### VECTORS #######

  parameters <- c("alpha", "beta.1", "beta.2", "beta.3", "beta.4")

  # Change UME to relate to parameters in model
  if (length(UME)==1) {
    if (UME==TRUE) {
      UME <- fun$params[which(fun$apool=="rel")]
    } else if (UME==FALSE) {
      UME <- c()
    }
  }

  write.check(fun=fun, positive.scale=positive.scale, intercept=intercept, link=link,
              rho=rho, covar=covar, omega=omega, sdscale=sdscale,
              class.effect=class.effect, UME=UME)

  model <- write.model()

  alphacode <- write.timecourse(model=model, fun=fun, intercept=intercept, positive.scale=positive.scale)

  timecourse <- alphacode[["timecourse"]]
  model <- alphacode[["model"]]

  model <- write.likelihood(model=model, timecourse=timecourse, rho=rho, covar=covar,
                            link=link, fun=fun, sdscale=sdscale)

  model <- write.beta(model=model, timecourse=timecourse, fun=fun,
                      UME=UME, class.effect=class.effect)


  if (corparam==TRUE) {
    model <- write.cor(model=model, fun=fun, omega=omega, class.effect = class.effect)
  }

  model <- add.funparams(model=model, fun=fun)

  model <- remove.loops(model)

  return(paste(model, sep="\n"))
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
write.check <- function(fun=tpoly(degree=1), positive.scale=TRUE, intercept=NULL, rho=0, covar=NULL,
                        omega=NULL, link="identity", sdscale=FALSE,
                        class.effect=list(), UME=c()) {

  # Run argument checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(link, choices=c("identity", "log", "smd"), null.ok=FALSE, add=argcheck)
  checkmate::assertList(class.effect, unique=FALSE, add=argcheck)
  checkmate::assertMatrix(omega, null.ok = TRUE, add=argcheck)
  checkmate::assertClass(fun, "timefun", add = argcheck)
  checkmate::assertLogical(sdscale, null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)


  if (!is.null(rho)) {
    # if (is.character(rho) & rho!="estimate") {
    #   stop("`rho` must either be assigned the value `estimate` to be estimated from the data or must be assigned a single numeric value")
    # }
    if (is.numeric(rho) & (rho < -1 | rho > 1)) {
      stop("Numeric values for `rho` cannot be outside the bounds of c(-1, 1)")
    }
    # if (is.null(covar)) {
    #   stop("`rho` has been assigned a value. Must therefore also define `covar`")
    # }
  } else {
    stop("`rho` cannot be NULL. To model no correlation between time-points set rho=0")
  }


  if (!all(names(class.effect) %in% fun$params)) {
    stop(paste0("The following list element names in `class.effect` do not match time-course parameter names in 'fun':\n",
                paste(names(class.effect)[!(names(class.effect) %in% fun$params)], collapse=", ")))
  }
  if (!all(UME %in% fun$params)) {
    stop(paste0("The following parameters specified in 'UME' do not match time-course parameter names in 'fun':\n",
                paste(UME[!(UME %in% fun$params)], collapse=", ")))
  }

  if (any(UME %in% names(class.effect))) {
    stop(paste0("UME / node-splitting model cannot be applied to parameters with class effects.\n",
                paste(UME[which(UME %in% names(class.effect))], collapse=", "), " have been assigned UME and class effects"))
  }

  if (any(fun$apool[which(fun$params %in% UME)] == "abs")) {
    stop("'UME' can only be specified for time-course parameters in 'fun' that have been modelled as pool='rel'")
  }

  if (any(fun$apool[which(fun$params %in% names(class.effect))] == "abs")) {
    stop("Class effects can only be specified for time-course parameters in 'fun' that have been modelled as pool='rel'")
  }

  # Check omega is symmetric positive definite matrix with correct dimensions
  if (!is.null(omega)) {
    err <- FALSE

    nrel <- sum(fun$apool %in% "rel" & !names(fun$apool) %in% names(class.effect))
    if (!all(dim(omega)==nrel)) {
      err <- TRUE
    }
    if (!isSymmetric(omega)) {
      err <- TRUE
    }
    if (any(eigen(omega)$values <= 0)) {
      err <- TRUE
    }
    if (err==TRUE) {
      stop("omega must be a symmetric positive definite matrix with dimensions equal to the number of\ntime-course parameters modelled as pool='rel'")
    }
  }

}






#' Write the basic JAGS model code for MBNMA to which other lines of model
#' code can be added
#'
#' @return A character vector of JAGS model code
#'
#'
write.model <- function() {
  model <- c(
    start="model{ 			# Begin Model Code",
    study="for(i in 1:NS){ # Run through all NS trials",
    arm="for (k in 1:narm[i]){ # Run through all arms within a study",
    fup="for (m in 1:fups[i]) {	# Run through all observations within a study",
    "resdev[i,k,m] <- pow((y[i,k,m] - theta[i,k,m]),2) * prec[i,k,m] # residual deviance for normal likelihood",
    "dev[i,k,m] <- -2* (log(pow((prec[i,k,m]/(2*3.14159)),0.5) * exp(-0.5*(pow((y[i,k,m]-theta[i,k,m]),2)*prec[i,k,m])))) # deviance for normal likelihood",
    "}",
    "",
    "resarmdev[i,k] <- sum(resdev[i,k,1:fups[i]])",
    "armdev[i,k] <- sum(dev[i,k,1:fups[i]])",
    "}",
    "",
    "resstudydev[i] <- sum(resarmdev[i, 1:narm[i]])",
    "studydev[i] <- sum(armdev[i, 1:narm[i]])",
    "",
    te="for(k in 2:narm[i]){ # Treatment effects",
    "}",
    "}",
    "",
    trt.prior="for (k in 2:NT){ # Priors on relative treatment effects",
    "}",
    "",
    class.prior="for (k in 2:Nclass){ # Priors on relative class effects",
    "}",
    "",
    "for (c in 1:(NT-1)) {",
    ume.prior="for (k in (c+1):NT) { # UME priors",
    "}",
    "}",
    "",
    "totresdev <- sum(resstudydev[])",
    "totdev <- sum(studydev[])",
    end="",
    "# Model ends",
    "}"
  )
}






#' Adds sections of JAGS code for an MBNMA model that correspond to alpha
#' parameters
#'
#' @param timecourse A character object that contains JAGS code for the
#'   time-course component of the model
#' @param model A character string representing the MBNMA model in JAGS code
#' @inheritParams mb.run
#'
#' @return A list of named elements: `model` is a character vector of JAGS MBNMA
#'   model code that includes alpha parameter components of the model
#'   `timecourse` is a character object that contains JAGS code for the
#'   time-course component of the model, for which alpha will be indexed
#'   correctly
#'
write.timecourse <- function(model, fun,
                        intercept, positive.scale) {

  timecourse <- fun$jags
  priors <- default.priors(fun)
  if (!is.null(intercept)) {
    if (intercept==TRUE) {

      # Insert prior for alpha
      model <- model.insert(model, pos=which(names(model)=="study"), priors[["alpha"]])

      if (positive.scale==TRUE) {
        timecourse <- paste0("exp(alpha[i]) + ", timecourse)
      } else {
        timecourse <- paste0("alpha[i] + ", timecourse)
      }
    }
  } else if (is.null(intercept)) {
    # Insert prior for alpha
    model <- model.insert(model, pos=which(names(model)=="study"), priors[["alpha"]])

    if (positive.scale==TRUE) {
      timecourse <- paste0("exp(ifelse(intercept[i]==1, alpha[i], 0)) + ", timecourse)
    } else {
      timecourse <- paste0("ifelse(intercept[i]==1, alpha[i], 0) + ", timecourse)
    }
  }

  # timecourse <- paste0("theta[i,k,m] + ", timecourse)
  #
  # # Insert timecourse into model
  # model <- model.insert(model, pos=which(names(model)=="fup"), timecourse)

  return(list("model"=model, "timecourse"=timecourse))
}






#' Insert element into model vector at desired location
#'
#' @noRd
model.insert <- function(a, pos, x){
  # dots <- list(...)
  # stopifnot(length(dots)==length(pos))
  # result <- vector("list",2*length(pos)+1)
  # result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
  # result[c(FALSE,TRUE)] <- dots
  # unlist(result)

  if (pos>length(a)) {
    stop("'pos' cannot be greater than length(a)")
  }
  start <- a[1:pos]
  end <- a[(pos+1):length(a)]
  return(c(start, x, end))
}






#' Adds sections of JAGS code for an MBNMA model that correspond to the
#' likelihood
#'
#' @inheritParams write.beta
#' @inheritParams write.timecourse
#' @inheritParams mb.run
#'
#' @return A character vector of JAGS MBNMA model code that includes likelihood
#'   components of the model
#'
write.likelihood <- function(model, timecourse, rho=0, covar="varadj",
                             link="identity", sdscale=FALSE, fun) {

  # Likelihoods
  norm.like <- c(
    "y[i,k,m] ~ dnorm(theta[i,k,m], prec[i,k,m])",
    "prec[i,k,m] <- pow(se[i,k,m], -2)"
  )

  mnorm.like <- c(
    "y[i,k,1:fups[i]] ~ dmnorm.vcov(theta[i,k,1:fups[i]], cov.mat[i,k,1:fups[i],1:fups[i]])"
  )

  # Linear predictor
  if (link=="identity") {
    predictor <- timecourse <- paste0("theta[i,k,m] <- ", timecourse)
  } else if (link=="log") {
    predictor <- timecourse <- paste0("log(theta[i,k,m]) <- ", timecourse)
  } else if (link=="smd") {
    predictor <- c(
      "phi[i,k,m] <- theta[i,k,m] * pool.sd[i]",
      paste0("theta[i,k,m] <- ", timecourse)
    )
    norm.like <- gsub("theta", "phi", norm.like)
    mnorm.like <- gsub("theta", "phi", mnorm.like)

    model <- gsub("theta", "phi", model) # Swap in resdev
  }


  # Write rho prior and multivariate code sections
  if (!is.null(rho)) {
    if (is.character(rho)) {
      rho.prior <- default.priors(fun)[["rho"]]
    } else if (is.numeric(rho)) {
      rho.prior <- paste0("rho <- ", rho)
    }
  }

  covar.cs <- c(
    "# Generates covariance matrix upper and lower triangles",
    "for (c in 1:(fups[i]-1)) {",
    "for (r in (c+1):fups[i]) {",
    "cov.mat[i,k,r,c] <- se[i,k,c]*se[i,k,r]*rho   # Lower triangle",
    "cov.mat[i,k,c,r] <- se[i,k,c]*se[i,k,r]*rho   # Upper triangle",
    "}",
    "}",
    "",
    "# Generates covariance matrix diagonals",
    "for (m in 1:fups[i]) {",
    "cov.mat[i,k,m,m] <- pow(se[i,k,m],2)",
    "}"
  )

  covar.ar1 <- c(
    "# Generates covariance matrix upper and lower triangles",
    "for (c in 1:(fups[i]-1)) {",
    "for (r in (c+1):fups[i]) {",
    "cov.mat[i,k,r,c] <- se[i,k,c]*se[i,k,r]*cor[i,r,c]   # Lower triangle",
    "cov.mat[i,k,c,r] <- se[i,k,c]*se[i,k,r]*cor[i,c,r]   # Upper triangle",
    "}",
    "}",
    "",
    "# Generates covariance matrix diagonals",
    "for (m in 1:fups[i]) {",
    "cov.mat[i,k,m,m] <- pow(se[i,k,m],2)",
    "}"
  )

  timepoint.corr <- c(
    "# Generates separate correlation coefficients for each two time points",
    "for (c in 1:(fups[i]-1)) {",
    "for (r in (c+1):fups[i]) {",
    "cor[i,r,c] <- pow(rho, (abs(time[i,r] - time[i,c])) / timedif.0[i])",
    "cor[i,c,r] <- pow(rho, (abs(time[i,r] - time[i,c])) / timedif.0[i])",
    "}",
    "}"
  )

  # Insert likelihood into model
  if (is.null(rho)) {
    model <- model.insert(model, pos=which(names(model)=="fup"), x=norm.like)
  } else if (!is.null(rho)) {
    model <- model.insert(model, pos=which(names(model)=="end"), x=rho.prior)

    if (covar %in% c("CS", "AR1")) {
      model <- model.insert(model, pos=which(names(model)=="arm"), x=mnorm.like)

      # Add covariance matrices
      if (covar=="CS") {
        model <- model.insert(model, pos=which(names(model)=="arm"), x=covar.cs)
      } else if (covar=="AR1") {
        model <- model.insert(model, pos=which(names(model)=="arm"), x=covar.ar1)
        model <- model.insert(model, pos=which(names(model)=="study"), x=timepoint.corr)
      }

      # Remove residual deviance calculations
      # Drop resdev, dev and totresdev
      model <- subset(model, !grepl("dev", model))
    }
    if (covar %in% c("varadj")) {
      # norm.like[1] <- gsub("(prec\\[i\\,k\\,m\\])", "\\1*(1-rho2)", norm.like[1])
      norm.like[1] <- gsub("prec\\[i\\,k\\,m\\] \\<\\- pow\\(se\\[i\\,k\\,m\\]",
                           "\\1*(1-rho2)", norm.like[1])

      model <- model.insert(model, pos=which(names(model)=="fup"), x=norm.like)

      model <- model.insert(model, pos=which(names(model)=="start"), x="rho2 <- rho*rho")
    }

  }

  # Add linear predictor
  model <- model.insert(model, pos=which(names(model)=="fup"), x=predictor)

  if (link=="smd" & sdscale==FALSE) {

    smd.sub <- c(
      "sd[i,k] <- se[i,k,1] * pow(n[i,k],0.5)",
      "nvar[i,k] <- (n[i,k]-1) * pow(sd[i,k],2)"
    )
    model <- model.insert(model, pos=which(names(model)=="arm"), x=smd.sub)

    pool.sd <- c(
      "df[i] <- sum(n[i,1:narm[i]]) - narm[i]",
      "pool.var[i] <- sum(nvar[i,1:narm[i]])/df[i]",
      "pool.sd[i] <- pow(pool.var[i], 0.5)"
    )
    model <- model.insert(model, pos=which(names(model)=="study"), x=pool.sd)
  }

  return(model)
}






#' Adds sections of JAGS code for an MBNMA model that correspond to beta
#' parameters
#'
#' @param model A character object of JAGS MBNMA model code
#' @param timecourse A character object representing the time-course used in the MBNMA model
#' @inheritParams mb.run
#'
#' @return A character vector of JAGS MBNMA model code that includes beta
#'   parameter components of the model
#'
write.beta <- function(model, timecourse, fun, UME, class.effect) {

  priors <- default.priors(fun)

  for (i in seq_along(fun$apool)) {
    if (fun$apool[i]=="rel") {

      # Split beta into relative effects
      model <- model.insert(model, pos=which(names(model)=="arm"),
                            x=paste0("beta.", i, "[i,k] <- mu.", i, "[i] + delta.", i, "[i,k]"))

      # Add prior for mu
      model <- model.insert(model, pos=which(names(model)=="study"),
                            x=priors[[paste0("mu.", i)]])

      # Add reference arm = 0 for delta
      model <- model.insert(model, pos=which(names(model)=="study"),
                            x=paste0("delta.", i, "[i,1] <- 0"))


      # For treatment effects across network
      if (!(fun$params[i] %in% UME)) {

        # Set Placebo effect equal to zero
        model <- model.insert(model, pos=which(names(model)=="start"),
                              x=paste0("d.", i, "[1] <- 0"))
      }

      if (!((fun$params[i] %in% names(class.effect)) | (fun$params[i] %in% UME))) {

        # Insert treatment effect prior
        model <- model.insert(model, pos=which(names(model)=="trt.prior"),
                              x=priors[[paste0("d.", i)]])

        # CLASS EFFECTS
      } else if (fun$params[i] %in% names(class.effect)) {

        # Set reference class effect equal to zero
        model <- model.insert(model, pos=which(names(model)=="start"),
                              x=paste0("D.", i, "[1] <- 0"))

        # Insert class effect priors
        model <- model.insert(model, pos=which(names(model)=="class.prior"),
                              x=priors[[paste0("D.", i)]])

        if (class.effect[[fun$params[i]]]=="common") {
          # Set trt effect equal to class effect
          model <- model.insert(model, pos=which(names(model)=="trt.prior"),
                                x=paste0("d.", i, "[k]  <- D.", i, "[class[k]]"))

        } else if (class.effect[[fun$params[i]]]=="random") {
          # Draw trt effect from class effect distribution
          model <- model.insert(model, pos=which(names(model)=="trt.prior"),
                                paste0("d.", i, "[k]  ~ dnorm(D.", i, "[class[k]], tau.D.", i, ")"))

          # Insert sd.D prior
          model <- model.insert(model, pos=which(names(model)=="end"),
                                x=c(priors[[paste0("sd.D.", i)]],
                                  paste0("tau.D.", i, " <- pow(sd.D.", i, ", -2)")
                                  )
                                )
        }

      }

      if (fun$params[i] %in% UME) {
        # Inset UME prior
        model <- model.insert(model, pos=which(names(model)=="ume.prior"),
                              x=priors[[paste0("dume.", i)]])

        # Set reference UME to zero
        model <- model.insert(model, pos=which(names(model)=="start"),
                              x=paste0("d.", i, "[1,1] <- 0"))

        if (fun$amethod[i]=="common") {
          # Insert UME common effect
          model <- model.insert(model, pos=which(names(model)=="te"),
                                x=paste0("delta.", i, "[i,k] <- d.", i, "[treat[i,1],treat[i,k]]"))

        } else if (fun$amethod[i]=="random") {
          # Insert UME random effect
          model <- model.insert(model, pos=which(names(model)=="te"),
                                x=paste0("delta.", i, "[i,k] ~ dnorm(d.", i, "[treat[i,1],treat[i,k]], tau.", i, ")"))

          # Insert prior for heterogeneity
          model <- model.insert(model, pos=which(names(model)=="end"),
                                x=c(priors[[paste0("sd.beta.", i)]],
                                  paste0("tau.", i, " <- pow(sd.beta.", i, ", -2)"))
                                )
        }
      } else {
        if (fun$amethod[i]=="common") {
          # Insert common treatment effect
          model <- model.insert(model, pos=which(names(model)=="te"),
                                x=paste0("delta.", i, "[i,k] <- d.", i, "[treat[i,k]] - d.", i, "[treat[i,1]]"))

        } else if (fun$amethod[i]=="random") {
          # Insert random treatment effect
          model <- model.insert(model, pos=which(names(model)=="te"),
                                x=c(paste0("delta.", i, "[i,k] ~ dnorm(md.", i, "[i,k], taud.", i, "[i,k])"),
                                  paste0("md.", i, "[i,k] <- d.", i, "[treat[i,k]] - d.", i, "[treat[i,1]] + sw.", i, "[i,k]"),
                                  paste0("taud.", i, "[i,k] <- tau.", i, " * 2*(k-1)/k"),
                                  paste0("w.", i, "[i,k] <- (delta.", i, "[i,k] - d.", i, "[treat[i,k]] + d.", i, "[treat[i,1]])"),
                                  paste0("sw.", i, "[i,k] <- sum(w.", i, "[i,1:(k-1)])/(k-1)")
                                )
                                )

          # Insert prior for heterogeneity
          model <- model.insert(model, pos=which(names(model)=="end"),
                                x=c(priors[[paste0("sd.beta.", i)]],
                                  paste0("tau.", i, " <- pow(sd.beta.", i, ", -2)"))
          )

          # Insert multi-arm correction zero
          model <- model.insert(model, pos=which(names(model)=="study"),
                                x=paste0("w.", i, "[i,1] <- 0"))
        }
      }
    }
  }

  # Absolute time-course parameters
  for (i in seq_along(fun$amethod)) {
    if ("abs" %in% fun$apool[i]) {

      if (grepl("[A-z]", fun$amethod[i])) {
      # if (is.character(fun$amethod[i])) {
        # Insert prior for absolute effect
        model <- model.insert(model, pos=which(names(model)=="end"),
                              x=priors[[paste0("beta.", i)]])


        if (fun$amethod[i]=="random") {
          # Insert distribution for random absolute effect
          model <- model.insert(model, pos=which(names(model)=="arm"),
                                x=paste0("i.beta.", i, "[i,k] ~ dnorm(beta.", i, ", prec.beta.", i, ")"))

          # Insert sd prior for random absolute effect
          model <- model.insert(model, pos=which(names(model)=="end"),
                                x=c(paste0("prec.beta.", i, " <- pow(sd.beta.", i, ", -2)"),
                                    priors[[paste0("sd.beta.", i)]])
          )
        }
      } else if (grepl("[0-9]", fun$amethod[i])) {
        # Insert fixed value for absolute effect
        model <- model.insert(model, pos=which(names(model)=="start"),
                              x=paste0("beta.", i, " <- ", as.numeric(fun$amethod[i])))

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
write.cor <- function(model, fun, omega=NULL, class.effect=list()) {

  if (length(class.effect)>0) {
    message("Class effects cannot be modelled with correlation between time-course relative effect parameters - correlation will be ignored")
  } else {

    sufparams <- which(fun$apool=="rel")

    # Check the number of parameters modelled using relative effects
    mat.size <- length(sufparams)
    if (mat.size>=2) {
      model <- write.cov.mat(model, sufparams=sufparams,
                             cor="estimate", cor.prior="rho",
                             omega=omega, fun=fun)
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
                          omega=NULL, fun) {

  priors <- default.priors(fun)

  if (any(c("itp", "emax") %in% fun$name) & FALSE %in% fun$p.expon) {
    d.zero <- 0.00001
  } else {
    d.zero <- 0
  }

  jagswish <- c(
    "for (r in 1:mat.size) {",
    paste0("d.prior[r] <- ", d.zero),
    "}",
    "",
    priors[["inv.R"]]
  )

  # jagsrho <- c(
  #   "for (r in 1:mat.size) {",
  #   paste0("d.prior[r] <- ", d.zero),
  #   "R[r,r] <- 1000    # Covariance matrix diagonals",
  #   "}",
  #   "",
  #   "for (r in 1:(mat.size-1)) {  # Covariance matrix upper/lower triangles",
  #   "for (c in (r+1):mat.size) {",
  #   "R[r,c] <- 1000*rho[1]   # Lower triangle",
  #   "R[c,r] <- 1000*rho[1]   # Upper triangle",
  #   "}",
  #   "}"
  # )
  jagsrho <- c(
    "for (r in 1:mat.size) {",
    paste0("d.prior[r] <- ", d.zero),
    "R[r,r] <- 1    # Covariance matrix diagonals",
    "L[r,r] <- pow(R[r,r],0.5) # Cholesky decomposition of covariance matrix",
    "for (c in 1:(r-1)) {  # Covariance matrix upper/lower triangles",
    "R[r,c] <- 1*rhoparam   # Lower triangle",
    "R[c,r] <- 1*rhoparam   # Upper triangle",
    "L[r,c] <- R[r,c] / L[c,c]",
    "L[c,r] <- 0",
    "}",
    "}"
  )

  mat.size <- length(sufparams)
  for (i in seq_along(sufparams)) {
    # Change d.1[k] ~ dnorm(0,0.001)  to   d.1[k] <- d.mult[1,k]
    model <- gsub(paste0("^d\\.", sufparams[i], "\\[k\\].+$"),
                  #paste0("^d\\.", sufparams[i], "\\[k\\] ~ [a-z]+\\([0-9]+(\\.[0-9]+)?,[0-9]+(\\.?[0-9]+)?\\)"),
                  paste0("d.", sufparams[i], "[k] <- mult[", i, ",k]"),
                  model
    )

    model <- gsub(paste0("^mu\\.", sufparams[i], "\\[i\\].+$"),
                  #paste0("^mu\\.", sufparams[i], "\\[i\\] ~ [a-z]+\\([0-9]+(\\.[0-9]+)?,[0-9]+(\\.?[0-9]+)?\\)"),
                  paste0("mu.", sufparams[i], "[i] <- mumult[i,", i, "]"),
                  model
    )
  }

  if (cor.prior=="wishart") {
    addcode <- jagswish

    # Insert multivariate normal dist (Wishart)
    model <- model.insert(model, pos=which(names(model)=="trt.prior"),
                          x=paste0("mult[1:", mat.size, ",k] ~ dmnorm(d.prior[], inv.R[1:", mat.size, ", 1:", mat.size, "])"))

    # Insert multivariate normal dist for mu (Wishart)
    model <- model.insert(model, pos=which(names(model)=="study"),
                          x=paste0("mumult[i,1:", mat.size, "] ~ dmnorm(d.prior[], muinv.R[1:", mat.size, ", 1:", mat.size, "])"))

    model <- model.insert(model, pos=which(names(model)=="end"),
                          x=paste0("muinv.R ~ dwish(omega[,], ", mat.size, ")"))


  } else if (cor.prior=="rho") {
    addcode <- jagsrho

    jagsrho.mu <- c(
      paste0("for (r in 1:", mat.size, ") {"),
      "mumult[i,r] <- d.prior[r] + L[r,1:r] %*% mu.z[i,1:r]",
      "}"
    )

    jagsrho.d <- c(
      paste0("for (r in 1:", mat.size, ") {"),
      "mult[r,k] <- d.prior[r] + L[r,1:r] %*% z[1:r,k]",
      "}"
    )

    # Insert correlated univariate chunks
    model <- model.insert(model, pos=which(names(model)=="study"),
                          x=jagsrho.mu)
    model <- model.insert(model, pos=which(names(model)=="trt.prior"),
                          x=jagsrho.d)

    # Insert prior for corparam
    model <- model.insert(model, pos=which(names(model)=="end"),
                          x=priors[["rhoparam"]])


    for (i in seq_along(sufparams)) {

      # Insert correlated univariate priors for mu
      model <- model.insert(model, pos=which(names(model)=="study"),
                            x=priors[[paste0("mu.z.",sufparams[i])]])

      # Insert correlated univariate priors for d
      model <- model.insert(model, pos=which(names(model)=="trt.prior"),
                            x=priors[[paste0("z.",sufparams[i])]])
    }
  }


  addcode <- gsub("mat\\.size", mat.size, addcode)

  # Insert covariance matrix
  model <- model.insert(model, pos=which(names(model)=="end"),
                        x=addcode)

  return(model)
}






#' Removes any loops from MBNMA model JAGS code that do not contain any
#' expressions
#'
#' @inheritParams write.beta
#'
#' @return A character vector of JAGS MBNMA model code that has had empty loops
#'   removed from it
#'
remove.loops <- function(model) {

  loops <- c(te=1, trt.prior=1, class.prior=1, umeloop=1)
  loops <- c("te", "trt.prior", "class.prior")

  for (i in seq_along(loops)) {
    index <- which(names(model)==loops[i])
    if (model[index+1]=="}") {
      model <- model[-(index:(index+1))]
    }
  }

  index <- which(names(model)=="ume.prior")
  if (model[index+1]=="}") {
    model <- model[-((index-1):(index+2))]
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
#' result <- mb.run(network, fun=tpoly(degree=1,
#'     pool.1="rel", method.1="random"))
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
  checkmate::assertCharacter(model)
  if (!(any(grepl("Begin Model Code", model)) & any(grepl("Model ends", model)))) {
    stop("'model' is not a character vector of MBNMA JAGS model code")
  }

  priorcode <- model[c(grep("^.+~ [A-z]+\\([-?0-9]", model),
                       grep("^.+~ [A-z]+\\(omega", model))]

  priorlist <- strsplit(priorcode, split=" +?~ +?")
  priors <- list()
  for (i in seq_along(priorlist)) {
    priorname <- unlist(strsplit(priorlist[[i]][1], split="\\["))[1]

    if (priorname %in% names(priors)) {
      priors[[priorname]] <- append(priors[[priorname]], priorlist[[i]][2])
    } else {
      priors[[priorname]] <- priorlist[[i]][2]
    }
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
#'   **using distributions as specified in JAGS syntax** (see \insertCite{jagsmanual;textual}{MBNMAtime}).
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
  checkmate::assertCharacter(model, null.ok=TRUE, add=argcheck)
  checkmate::assertList(priors, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (!is.null(mbnma) & !is.null(model)) {
    stop("Must provide EITHER an existing MBNMA model (using `mbnma`) OR MBNMA JAGS code (using `model`)")
  }

  if (!is.null(mbnma)) {
    model <- mbnma$model.arg$jagscode
  } else if (is.null(model)) {
    stop("Must provide EITHER an existing MBNMA model (using `mbnma`) OR MBNMA JAGS code (using `model`)")
  }

  if (!(any(grepl("Begin Model Code", model)) & any(grepl("Model ends", model)))) {
    stop("'model' is not a character vector of MBNMA JAGS model code")
  }

  for (i in seq_along(priors)) {
    # Checks
    if (length(grep(paste0("^( +)?", names(priors)[i]), model))==0) {
      stop("Prior named ", names(priors)[i], " not found in the model code. Check priors currently present in model code using get.prior()")
      # } else if (length(grep(paste0("^( +)?", names(priors)[i]), model))>1) {
      #   stop("Prior named ", names(priors)[i], " has matched on multiple instances in the model code. Check priors currently present in model code using get.prior()")
    }

    if (!all(grepl("^d[a-z.]*\\(", priors[[i]]))) {
      stop("Prior named ", names(priors)[i], " does not follow JAGS distribution syntax\nSee JAGS manual: https://people.stat.sc.edu/hansont/stat740/jags_user_manual.pdf")
    }

    line <- grep(paste0("^( +)?", names(priors)[i], ".+~"), model)
    state <- model[line]

    if (length(priors[[i]])==1) {
      model[line] <- gsub("(^.+~ )(.+$)", paste0("\\1", priors[[i]]), state)

    } else {
      # What if length of priors[[i]]>1 ?
      # Find previous { in code and add priors as new lines there

      # Identifies loop above which to insert
      insert <- max(grep("\\{", model)[grep("\\{", model) < line])

      # Indentifies starting index in the loop (e.g. from 1: or 2:)
      loopind <- as.numeric(gsub("\\D", "", model[insert]))

      # Creates vector of priors
      priors.insert <- paste0(names(priors)[i], "[",
                              loopind:(length(priors[[i]])+loopind-1),
                              "] ~ ", priors[[i]])

      # Drop previous prior line
      model <- model[-line]
      model <- c(model[1:(insert-1)],
                 priors.insert,
                 model[insert:length(model)])

    }
  }

  # Cut irrelevant section from JAGS code
  start <- grep("^model\\{", model)
  end <- grep("# Model ends", model) + 1

  #model <- paste(model[start:end], collapse="\n")
  model <- model[start:end]

  return(model)
}





#' Add named function parameters to the model
#' @noRd
add.funparams <- function(model, fun) {
  for (i in seq_along(fun$params)) {
    if (!grepl("beta", fun$params[i])) {
      if ("rel" %in% fun$apool[i]) {
        model <- gsub(paste0("(?<![A-z])d\\.", i, "\\["), paste0(names(fun$apool)[i], "["), model, perl=TRUE)
      } else if ("abs" %in% fun$apool[i] & is.character(fun$amethod[i])) {
        model <- gsub(paste0("beta\\.", i), names(fun$apool)[i], model)
      }
      model <- gsub(paste0("sd\\.beta.", i), paste0("sd.", names(fun$apool)[i]), model)

      # Check for class effects
      model <- gsub(paste0("D\\.", i, "\\["), paste0(toupper(names(fun$apool)[i]), "["), model)
      model <- gsub(paste0("sd\\.D\\.", i), paste0("sd.", toupper(names(fun$apool)[i])), model)

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
#' # Write a log-linear time-course MBNMA synthesis model with:
#' # Common effects for synthesis of mu
#' # Modelled as ratio of means
#' model <- write.ref.synth(fun=tloglin(pool.rate="rel", method.rate="common"),
#'   mu.synth="common", link="log")
#'
#' cat(model) # Concatenates model representations making code more easily readable
#'
#' @export
write.ref.synth <- function(fun=tpoly(degree = 1), link="identity",
                            positive.scale=TRUE, intercept=TRUE, rho=0, covar="varadj",
                            mu.synth="random",
                            priors=NULL) {

  ####### VECTORS #######

  write.check(fun=fun, positive.scale=positive.scale, intercept=intercept,
              rho=rho, covar=covar)

  model <- write.model()

  alphacode <- write.timecourse(model=model, fun=fun, intercept=intercept, positive.scale=positive.scale)

  timecourse <- alphacode[["timecourse"]]
  model <- alphacode[["model"]]

  model <- write.likelihood(model=model, timecourse=timecourse, rho=rho, covar=covar, link=link, fun=fun)

  model <- write.beta.ref(model=model, timecourse=timecourse, fun=fun, mu.synth=mu.synth)

  model <- remove.loops(model)

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






#' Adds sections of JAGS code for an MBNMA reference synthesis model that
#' correspond to beta parameters
#' @noRd
write.beta.ref <- function(model, timecourse, fun,
                           mu.synth="random"
) {

  priors <- default.priors(fun)

  for (i in seq_along(fun$apool)) {
    if ("rel" %in% fun$apool[i]) {

      model <- model.insert(model, pos=which(names(model)=="end"),
                            x=priors[[paste0("m.mu.",i)]])

      if (mu.synth=="common") {
        model <- gsub(paste0("beta\\.",i, "\\[i\\,k\\]"), paste0("mu.",i), model)
      } else if (mu.synth=="random") {
        model <- gsub(paste0("beta\\.",i), paste0("i.mu.",i), model)

        model <- model.insert(model, pos=which(names(model)=="arm"),
                              x=paste0("i.mu.", i, "[i,k] ~ dnorm(mu.", i, ", tau.mu.", i, ")"))

        model <- model.insert(model, pos=which(names(model)=="end"),
                              x=priors[[paste0("sd.mu.", i)]])
        model <- model.insert(model, pos=which(names(model)=="end"),
                              x=paste0("tau.mu.", i, " <- pow(sd.mu.", i, ", -2)"))
      }
    } else if ("abs" %in% fun$apool[i]) {
      # Insert prior for absolute effect
      model <- model.insert(model, pos=which(names(model)=="end"),
                            x=priors[[paste0("beta.", i)]])


      if (fun$amethod[i]=="random") {
        # Insert distribution for random absolute effect
        model <- model.insert(model, pos=which(names(model)=="arm"),
                              x=paste0("i.beta.", i, "[i,k] ~ dnorm(beta.", i, ", prec.beta.", i, ")"))

        # Insert sd prior for random absolute effect
        model <- model.insert(model, pos=which(names(model)=="end"),
                              x=c(paste0("prec.beta.", i, " <- pow(sd.beta.", i, ", -2)"),
                                  priors[[paste0("sd.beta.", i)]])
        )
      }
    }
  }

  # Create dummy variables to avoid JAGS warnings
  model <- model.insert(model, pos=which(names(model)=="start"),
                        x="dummy1 <- NT")
  model <- model.insert(model, pos=which(names(model)=="start"),
                        x="dummy2 <- treat")

  return(model)
}










#' Write standard NMA model in JAGS
#'
#' @inheritParams mb.run
#' @inheritParams plot.mb.predict
#'
#' @noRd
write.nma <- function(method="common", link="identity", sdscale=FALSE) {
  model <- c(
    start="model{ 			# Begin Model Code",
    "d[1] <- 0",
    study="for(i in 1:NS){ # Run through all NS trials",
    "mu[i] ~ dnorm(0,0.0001)",
    "delta[i,1] <- 0",
    arm="for (k in 1:narm[i]){ # Run through all arms within a study",
    "resdev[i,k] <- pow((y[i,k] - theta[i,k]),2) * prec[i,k] # residual deviance for normal likelihood",
    "}",
    "resstudydev[i] <- sum(resdev[i, 1:narm[i]])",
    te="for(k in 2:narm[i]){ # Treatment effects",
    "}",
    "}",
    "",
    trt.prior="for (k in 2:NT){ # Priors on relative treatment effects",
    "d[k] ~ dnorm(0,0.0001)",
    "}",
    "totresdev <- sum(resstudydev[])",
    "",
    end="",
    "# Model ends",
    "}"
  )

  arm.insert <- c("y[i,k] ~ dnorm(theta[i,k], prec[i,k])",
                  "prec[i,k] <- pow(se[i,k], -2)")

  if (link=="identity") {
    arm.insert <- append(arm.insert,
                         c("theta[i,k] <- mu[i] + delta[i,k]"))

  } else if (link=="log") {
    arm.insert <- append(arm.insert,
                         c("log(theta[i,k]) <- mu[i] + delta[i,k]"))

  } else if (link=="smd") {
    arm.insert <- c("y[i,k] ~ dnorm(phi[i,k], prec[i,k])",
                    "prec[i,k] <- pow(se[i,k], -2)",
                    "phi[i,k] <- theta[i,k] * pool.sd[i]",
                    "theta[i,k] <- mu[i] + delta[i,k]")

    if (sdscale==FALSE) {
      arm.insert <- append(arm.insert, c("sd[i,k] <- se[i,k] * pow(n[i,k],0.5)",
                                         "nvar[i,k] <- (n[i,k]-1) * pow(sd[i,k],2)"
                                         ))

      insert <- c("df[i] <- sum(n[i,1:narm[i]]) - narm[i]",
                  "pool.var[i] <- sum(nvar[i,1:narm[i]])/df[i]",
                  "pool.sd[i] <- pow(pool.var[i], 0.5)")
      model <- model.insert(model, pos=which(names(model)=="start"),
                            x=insert)
    }
  }

  model <- model.insert(model, pos=which(names(model)=="arm"),
                        x=arm.insert)

  if (method=="common") {
    te.insert <- c("delta[i,k] <- md[i,k]",
                   "md[i,k] <- d[treat[i,k]] - d[treat[i,1]]")

  } else if (method=="random") {
    te.insert <- c("delta[i,k] ~ dnorm(md[i,k], taud[i,k])",
                   "md[i,k] <- d[treat[i,k]] - d[treat[i,1]] + sw[i,k]",
                   "taud[i,k] <- tau *2*(k-1)/k",
                   "w[i,k] <- (delta[i,k] - d[treat[i,k]] + d[treat[i,1]])",
                   "sw[i,k] <- sum(w[i,1:(k-1)])/(k-1)")

    # Insert at start
    model <- model.insert(model, pos=which(names(model)=="study"),
                          x="w[i,1] <- 0")

    #Insert at end
    model <- model.insert(model, pos=which(names(model)=="end"),
                          x=c("tau <- pow(sd,-2)",
                              "sd ~ dnorm(0,0.05) T(0,)"))
  }
  model <- model.insert(model, pos=which(names(model)=="te"),
                        x=te.insert)

  return(model)
}






#' Write non-parametric random walk model JAGS code
#'
#' Writes JAGS code for a Bayesian non-parametric model that splits the
#' data into different time-bins and assumes a random walk process for treatment
#' effects between them.
#'
#' @param method A character object that can take the value `"common"` or `"random"` that
#'   specifies the the type of pooling to use for synthesis.
#' @inheritParams mb.run
#'
#' @return A single long character string containing the JAGS model generated
#'   based on the arguments passed to the function.
#'
#' @examples
#' # Write a common effects non-paramtric random walk model
#' write.rw(method="common")
#'
#' @noRd
write.rw <- function(method="common", link="identity") {

  model <- c(
    "model{ \t\t\t# Begin Model Code",
    "d[1] <- 0",
    "",
    "for(i in 1:NS){ # Run through all NS trials",
    "mu[i,1] ~ dnorm(0,0.0001)",
    "delta[i,1,1] <- 0",
    "for (k in 1:narm[i]){ # Run through all arms within a study",
    "for (m in 1:nbin) { # Run through each bin",
    "y[i,k,m] ~ dnorm(theta[i,k,m], prec[i,k,m])",
    "prec[i,k,m] <- pow(se[i,k,m], -2)",
    "theta[i,k,m] <- mu[i,m] + delta[i,k,m]",
    "dev[i,k,m] <- pow((y[i,k,m] - theta[i,k,m]),2) * prec[i,k,m] # residual deviance for normal likelihood",
    "}",
    "",
    "resdev[i,k] <- sum(dev[i,k, 1:nbin])",
    "}",
    "",
    "for (m in 2:nbin) { # Run through each bin",
    "delta[i,1,m] <- 0",
    "mu[i,m] ~ dnorm(mu[i,m-1],tau.rw)",
    "",
    "for(k in 2:narm[i]){ # Treatment effects",
    "delta[i,k,m] ~ dnorm(delta[i,k,m-1],tau.rw)",
    "}",
    "}",
    "",
    "resstudydev[i] <- sum(resdev[i, 1:narm[i]])",
    "",

    # Common effects
    "for(k in 2:narm[i]){ # Treatment effects",
    "delta[i,k,1] <- md[i,k]",
    "md[i,k] <- d[treat[i,k]] - d[treat[i,1]]",
    "}",

    # # Random effects
    # "for(k in 2:narm[i]){ # Treatment effects",
    # "delta[i,k,1] ~ dnorm(md[i,k], taud[i,k])",
    # "md[i,k] <- d[treat[i,k]] - d[treat[i,1]] + sw[i,k]",
    # "taud[i,k] <- tau *2*(k-1)/k",
    # "w[i,k] <- (delta[i,k,1] - d[treat[i,k]] + d[treat[i,1]])",
    # "sw[i,k] <- sum(w[i,1:(k-1)])/(k-1)",
    # "}",

    "",
    "}",
    "",
    "for (k in 2:NT){ # Priors on relative treatment effects",
    "d[k] ~ dnorm(0,0.0001)",
    "}",
    "totresdev <- sum(resstudydev[])",
    "",
    "tau.rw <- pow(sd.rw,-2)",
    "sd.rw ~ dnorm(0,0.05) T(0,)",
    "",

    "tau <- pow(sd,-2)",
    "sd ~ dnorm(0,0.05) T(0,)",
    "# Model ends",
    "}"
    )


return(model)
}




#' Sets default priors for JAGS model code
#'
#' This function creates JAGS code snippets for default MBNMA model priors.
#'
#' @inheritParams mb.run
#'
#' @return A list, each element of which is a named JAGS snippet
#'   corresponding to a prior in the MBNMA JAGS code.
#'
#' @examples
#' \donttest{
#' default.priors(fun=temax())
#'
#' default.priors(fun=titp(p.expon=TRUE))
#' }
#'
#' @export
default.priors <- function(fun=tloglin()) {

  sufparams <- which(fun$apool=="rel")

  priors <- list(
    rho = "rho ~ dunif(0,1)",
    alpha = "alpha[i] ~ dnorm(0,0.0001)",
    inv.R = "inv.R ~ dwish(omega[,], mat.size)",
    rhoparam = "rhoparam ~ dunif(-1,1)"
  )

  for (i in 1:4) {
    priors[[paste0("mu.",i)]] <- paste0("mu.", i, "[i] ~ dnorm(0,0.0001)")
    priors[[paste0("m.mu.",i)]] <- paste0("mu.", i, " ~ dnorm(0,0.0001)")
    priors[[paste0("d.",i)]] <- paste0("d.", i, "[k] ~ dnorm(0,0.001)")
    priors[[paste0("dume.",i)]] <- paste0("d.", i, "[c,k] ~ dnorm(0,0.001)")
    priors[[paste0("D.",i)]] <- paste0("D.", i, "[k] ~ dnorm(0,0.001)")
    priors[[paste0("beta.",i)]] <- paste0("beta.", i, " ~ dnorm(0,0.0001)")

    priors[[paste0("sd.mu.",i)]] <- paste0("sd.mu.", i, " ~ dnorm(0,0.05) T(0,)")
    priors[[paste0("sd.d.",i)]] <- paste0("sd.d.", i, " ~ dnorm(0,0.05) T(0,)")
    priors[[paste0("sd.D.",i)]] <- paste0("sd.D.", i, " ~ dnorm(0,0.05) T(0,)")
    priors[[paste0("sd.beta.",i)]] <- paste0("sd.beta.", i, " ~ dnorm(0,0.05) T(0,)")
  }

  if ((fun$name %in% c("itp") | (fun$name %in% "emax"))) {

    for (i in 2:3) {
      priors[[paste0("mu.",i)]] <- paste0("mu.", i, "[i] ~ dnorm(0.00001,0.0001) T(0,)")
      priors[[paste0("m.mu.",i)]] <- paste0("mu.", i, " ~ dnorm(0.00001,0.0001) T(0,)")
      priors[[paste0("d.",i)]] <- paste0("d.", i, "[k] ~ dnorm(0.00001,0.0001) T(0,)")
      priors[[paste0("dume.",i)]] <- paste0("d.", i, "[c,k] ~ dnorm(0.00001,0.0001) T(0,)")
      priors[[paste0("D.",i)]] <- paste0("D.", i, "[k] ~ dnorm(0.00001,0.0001) T(0,)")
      priors[[paste0("beta.",i)]] <- paste0("beta.", i, " ~ dnorm(0.00001,0.0001) T(0,)")
    }
  }

  # For rho correlation between parameters
  for (i in seq_along(sufparams)) {

    if (((fun$name %in% c("itp") | (fun$name %in% "emax"))) & sufparams[i] %in% c(2:3)) {
      priors[[paste0("mu.z.",sufparams[i])]] <- paste0("mu.z[i,", i, "] ~ dnorm(0,0.0001) T(-d.prior[", i, "],)")
      priors[[paste0("z.",sufparams[i])]] <- paste0("z[", i, ",k] ~ dnorm(0,0.0001) T(-d.prior[", i, "],)")

    } else {
      priors[[paste0("mu.z.",sufparams[i])]] <- paste0("mu.z[i,", i, "] ~ dnorm(0,0.0001)")
      priors[[paste0("z.",sufparams[i])]] <- paste0("z[", i, ",k] ~ dnorm(0,0.0001)")
    }
  }

  return(priors)
}
