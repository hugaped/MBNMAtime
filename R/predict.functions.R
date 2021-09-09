# Functions for predicting study mean responses
# Author: Hugo Pedder
# Date created: 2018-09-10






#' Get MBNMA model values
#'
#' Extracts specific information required for prediction from a time-course
#' MBNMA model
#'
#' @inheritParams predict.mbnma
#' @inheritParams ref.synth
#'
#' @return A list containing named elements that correspond to different
#'   time-course parameters in `mbnma`. These elements contain MCMC results
#'   either taken directly from `mbnma` or (in the case of random time-course
#'   parameters specified as `method="random"`) randomly
#'   generated using parameter values estimated in `mbnma`.
#'
#'   Additional elements contain the following values:
#'   * `timecourse` A character object that specifies the time-course used in `mbnma` in terms of
#'   alpha, beta, mu, d and time. Consistency relative time-course parameters
#'   are specified in terms of mu and d.
#'   * `time.params` A character vector
#'   that indicates the different time-course parameters that are required for
#'   the prediction
#'
#'   @noRd
get.model.vals <- function(mbnma, E0=0, level="treatments") {

  # Check that correct parameters are monitored
  genparams <- gen.parameters.to.save(fun=mbnma$model.arg$fun, model=mbnma$model.arg$jagscode)
  if (!all(genparams %in% mbnma$parameters.to.save)) {
    stop(crayon::red(crayon::bold("Parameters required for estimation of time-course relationship not monitored in model.\nMust include time-course parameters in 'parameters.to.save'")))
  }

  model.vals <- list()
  time.params <- "alpha"
  #mu.prior <- vector()

  n <- mbnma$BUGSoutput$n.sims

  # Assign E0 to alpha in model.vals
  alpha <- eval(parse(text=E0))
  if (length(alpha)==1) {
    model.vals[["alpha"]] <- rep(alpha, n)
  } else if (length(alpha)>1) {
    model.vals[["alpha"]] <- alpha
  }

  # Remove indices from timecourse for betas and time
  timecourse <- gsub("\\[i\\,[a-z]\\]", "", mbnma$model.arg$fun$jags) # Remove [i,k] from betas and [i,m] from time
  timecourse <- gsub("i\\.", "", timecourse) # Remove i from time and spline
  timecourse <- paste0("alpha + ", timecourse)


  sims.matrix <- mbnma$BUGSoutput$sims.matrix

  fun <- mbnma$model.arg$fun
  params <- mbnma$model.arg$fun$params

  for (i in seq_along(params)) {
    if ("abs" %in% fun$apool[i]) {
      if ("common" %in% fun$amethod[i]) {
        # Store matrix of MCMC iterations to list
        model.vals[[fun$bname[i]]] <-
          sims.matrix[,grepl(paste0("^", params[i]), colnames(sims.matrix))]

        # Add beta parameters to the vector of time-course parameters
        time.params <- append(time.params, fun$bname[i])
      } else if ("random" %in% fun$amethod[i]) {
        # Store matrix of beta values generated from random distribution determined by model parameters
        len <- sum(grepl(paste0("^", params[i]), colnames(sims.matrix)))
        mat <- array(dim=c(n, len, 2))
        mat[,,1] <- sims.matrix[,grepl(paste0("^", params[i]), colnames(sims.matrix))]
        mat[,,2] <- sims.matrix[,grepl(paste0("^sd\\.", params[i]), colnames(sims.matrix))]
        mat <- apply(mat, MARGIN=c(1,2), FUN=function(x) stats::rnorm(1, x[1], x[2]))

        model.vals[[fun$bname[i]]] <- mat

        # Add beta parameters to the vector of time-course parameters
        time.params <- append(time.params, fun$bname[i])

      } else if (is.numeric(fun$amethod[i])) {
        stop(paste0("Common time-course parameter for", params[i], " currently not supported for prediction"))
      }
    } else if ("rel" %in% fun$apool[i]) {

      # Ammend time-course equation
      timecourse <- gsub(fun$bname[i], paste0("(mu.", i, " + d.", i, ")"), timecourse)

      # Add mu to list of parameters
      time.params <- append(time.params, paste0("mu.", i))
      time.params <- append(time.params, paste0("d.", i))
      #mu.prior <- append(mu.prior, paste0("mu.", i))

      # If class effects are present
      if (params[i] %in% names(mbnma$model.arg$class.effect)) {
        findd <- ifelse(grepl("beta", params[i]), paste0("^D\\.", i), paste0("^", toupper(params[i])))

        if (level=="treatments") {

          tperc <- table(mbnma$network$classkey$class)
          mat <- sims.matrix[,grepl(findd, colnames(sims.matrix))]
          if (ncol(mat)!=length(tperc)) {
            stop("Classes in 'network$classkey' do not match with those monitored in 'mbnma'")
          }

          # Duplicate class columns for treatments within a class
          mcmcmat <- as.matrix(mat[,1], ncol=1)
          for (k in 2:length(tperc)) {
            mcmcmat <- cbind(mcmcmat, matrix(rep(mat[,k], tperc[k]), ncol=tperc[k]))
          }

          if ("random" %in% mbnma$model.arg$class.effect[[params[i]]]) {
            # Store matrix of beta values generated from random distribution determined by model parameters
            mcmcarray <- array(dim=c(n, ncol(mcmcmat), 2))
            mcmcarray[,,1] <- mcmcmat
            mcmcarray[,,2] <- sims.matrix[,grepl(paste0("^sd\\.", substr(findd, 2, nchar(findd))), colnames(sims.matrix))]
            mcmcarray[,2:ncol(mcmcmat),] <-
              apply(mcmcarray[,2:ncol(mcmcmat),], MARGIN=c(1,2), FUN=function(x) stats::rnorm(1, x[1], x[2]))

            mcmcmat <- mcmcarray[,,1]
          }

          # Store MCMC results for relevant parameters
          model.vals[[paste0("d.", i)]] <- mcmcmat
        }

        if (level=="classes") {

          if ("common" %in% mbnma$model.arg$class.effect[[params[i]]]) {
            # Store MCMC results for relevant parameters
            model.vals[[paste0("d.", i)]] <- sims.matrix[,grepl(findd, colnames(sims.matrix))]
          } else if ("random" %in% mbnma$model.arg$class.effect[[params[i]]]) {
            # Store matrix of beta values generated from random distribution determined by model parameters
            len <- sum(grepl(findd, colnames(sims.matrix)))
            mat <- array(dim=c(n, len, 2))
            mat[,,1] <- sims.matrix[,grepl(findd, colnames(sims.matrix))]
            mat[,,2] <- sims.matrix[,grepl(paste0("^sd\\.", substr(findd, 2, nchar(findd))), colnames(sims.matrix))]
            mat[,2:len,] <- apply(mat[,2:len,], MARGIN=c(1,2), FUN=function(x) stats::rnorm(1, x[1], x[2]))

            model.vals[[paste0("d.", i)]] <- mat[,,1]
          }
        }

      # If class effects are not present
      } else {
        time.params <- append(time.params, paste0("d.", i))
        findd <- ifelse(grepl("beta", params[i]), paste0("^d\\.", i), paste0("^", params[i]))

        if ("common" %in% fun$amethod[i]) {

          # Store MCMC results for relevant parameters
          model.vals[[paste0("d.", i)]] <- sims.matrix[,grepl(findd, colnames(sims.matrix))]

        } else if ("random" %in% fun$amethod[i]) {

          # Store matrix of beta values generated from random distribution determined by model parameters
          len <- sum(grepl(findd, colnames(sims.matrix)))
          mat <- array(dim=c(n, len, 2))
          mat[,,1] <- sims.matrix[,grepl(findd, colnames(sims.matrix))]
          mat[,,2] <- sims.matrix[,grepl(paste0("^sd\\.", params[i]), colnames(sims.matrix))]
          mat[,2:len,] <- apply(mat[,2:len,], MARGIN=c(1,2), FUN=function(x) stats::rnorm(1, x[1], x[2]))

          model.vals[[paste0("d.", i)]] <- mat[,,1]

        }
      }
    } else {
      stop(paste0(params[i], " has not been monitored in the model but is necessary for prediction."))
    }
  }

  #model.vals[["mu.prior"]] <- mu.prior
  model.vals[["timecourse"]] <- timecourse
  model.vals[["time.params"]] <- unique(time.params)

  return(model.vals)
}











#' Synthesise single arm studies with repeated observations of the same
#' treatment over time
#'
#' Synthesises single arm studies with repeated measures by applying a
#' particular time-course function. Used in predicting mean responses from a
#' time-course MBNMA. The same parameterisation of the time course must be used
#' as in the MBNMA.
#'
#' @inheritParams predict.mbnma
#' @inheritParams R2jags::jags
#' @inheritParams mb.run
#' @param mbnma An S3 object of class `"mbnma"` generated by running
#' a time-course MBNMA model
#' @param data.ab A data frame of arm-level data in "long" format containing the
#'   columns:
#'   * `studyID` Study identifiers
#'   * `time` Numeric data indicating follow-up times
#'   * `y` Numeric data indicating the mean response for a given observation
#'   * `se` Numeric data indicating the standard error for a given observation
#'
#' @details `data.ab` can be a collection of studies that closely resemble the
#'   population of interest intended for the prediction, which could be
#'   different to those used to estimate the MBNMA model, and could be include
#'   single arms of RCTs or observational studies. If other data is not
#'   available, the data used to estimate the MBNMA model can be used by
#'   selecting only the studies and arms that specify the network reference
#'   treatment responses.
#'
#' @return A list of named elements corresponding to each time-course parameter
#'   within an MBNMA model that contain the median posterior value for the
#'   network reference treatment response.
#'
#' @examples
#' \donttest{
#' # Create an mb.network object from a dataset
#' network <- mb.network(osteopain)
#'
#' # Run an MBNMA model with an Emax time-course
#' emax <- mb.run(network,
#'   fun=temax(pool.emax="rel", method.emax="common",
#'     pool.et50="abs", method.et50="random"))
#'
#' # Generate a set of studies with which to estimate the network reference treatment response
#' paindata.ref <- osteopain[osteopain$treatname=="Placebo_0",]
#'
#' # Estimate the network reference treatment effect using common effects meta-analysis
#' ref.synth(data.ab=paindata.ref, mbnma=emax, synth="common")
#'
#' # Estimate the network reference treatment effect using random effects meta-analysis
#' ref.synth(data.ab=paindata.ref, mbnma=emax, synth="random")
#' }
#'
#' @export
ref.synth <- function(data.ab, mbnma, synth="common",
                      link=mbnma$model.arg$link,
                      n.iter=mbnma$BUGSoutput$n.iter,
                      n.burnin=mbnma$BUGSoutput$n.burnin,
                      n.thin=mbnma$BUGSoutput$n.thin,
                      n.chains=mbnma$BUGSoutput$n.chains,
                      ...) {

  # First need to validate data.frame to check dataset is in correct format...maybe another function for this
  # Change it to correct format if it is not already
  #data.ab <- ref.validate(data.ab)[["data"]]
  data.ab <- ref.validate(data.ab)[["data.ab"]]

  # Run checks
  argcheck <- checkmate::makeAssertCollection()

  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertChoice(synth, choices=c("random", "common"), add=argcheck)
  checkmate::assertInt(n.iter, lower=1, add=argcheck)
  checkmate::assertInt(n.burnin, lower=1, add=argcheck)
  checkmate::assertInt(n.thin, lower=1, add=argcheck)
  checkmate::assertInt(n.chains, lower=1, add=argcheck)

  checkmate::reportAssertions(argcheck)

  # To get model for meta-analysis of placebo must create v similar model
  #to study model
  # Do all the mb.write bits but without the consistency bits

  jagsmodel <- write.ref.synth(fun=mbnma$model.arg$fun, positive.scale=mbnma$model.arg$positive.scale,
                               intercept=mbnma$model.arg$intercept,
                               rho=mbnma$model.arg$rho, covar=mbnma$model.arg$covar,
                               mu.synth=synth,
                               priors=mbnma$model.arg$priors
  )

  parameters.to.save <- vector()
  for (i in seq_along(mbnma$model.arg$fun$apool)) {
    if ("rel" %in% mbnma$model.arg$fun$apool[i]) {
      parameters.to.save <- append(parameters.to.save, paste0("mu.",i))
      if (synth=="random") {
        parameters.to.save <- append(parameters.to.save, paste0("sd.mu.",i))
      }
    }
  }

  jags.result <- mb.jags(data.ab, link=mbnma$model.arg$link,
                            model=jagsmodel, fun=mbnma$model.arg$fun,
                            rho=mbnma$model.arg$rho, covar=mbnma$model.arg$covar,
                            parameters.to.save=parameters.to.save,
                            n.iter=n.iter, n.burnin=n.burnin,
                            n.thin=n.thin, n.chains=n.chains,
                            ...)[["jagsoutput"]]

  if (any(jags.result$BUGSoutput$summary[,
                                         colnames(jags.result$BUGSoutput$summary)=="Rhat"
                                         ]>1.02)) {
    warning("Rhat values for parameter(s) in reference treatment synthesis model are >1.02. Suggest running for more iterations.")
  }

  return(jags.result)

}







#' Checks the validity of ref.resp if given as data frame
#'
#' Ensures `ref.resp` takes the correct form to allow for synthesis of network
#' reference treatment response if data is provided for meta-analysis
#'
#' @inheritParams ref.synth
ref.validate <- function(data.ab) {

  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, any.missing=FALSE, add=argcheck)
  checkmate::assertNames(names(data.ab), must.include = c("studyID", "y", "se", "time"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Sort data.ab
  data.ab <- dplyr::arrange(data.ab, studyID, time)

  # if (anyMissing(data.ab)) {
  #   stop("Data frame for synthesis of reference treatment contains NA values")
  # }

  message("Data frame must contain only data from reference treatment")

  #### Prepare data frame ####
  # Add arm index (=1 since only one arm in each study)
  data.ab[["arm"]] <- 1
  data.ab[["narm"]] <- 1
  data.ab[["treatment"]] <- 1

  # Ensuring studies are numbered sequentially
  if (!is.numeric(data.ab[["studyID"]])) {
    #message("Studies being recoded to allow sequential numbering")
    data.ab <- transform(data.ab,studyID=as.numeric(factor(studyID, levels=as.character(unique(data.ab$studyID)))))
    data.ab <- dplyr::arrange(data.ab, studyID, time)
  } else if (all(abs(diff(data.ab[["studyID"]])) != TRUE)) {
    #message("Studies being recoded to allow sequential numbering")
    data.ab <- transform(data.ab,studyID=as.numeric(factor(studyID, levels=as.character(unique(data.ab$studyID)))))
    data.ab <- dplyr::arrange(data.ab, studyID, time)
  }

  data.ab <- add_index(data.ab, reference=1)

  return(data.ab)

}



