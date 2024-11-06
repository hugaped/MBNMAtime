# Functions for time-course
# Author: Hugo Pedder
# Date created: 2021-01-04




#' Integrated Two-Component Prediction (ITP) function
#'
#' Similar parameterisation to the Emax model but with non-asymptotic maximal effect (Emax). Proposed
#' by proposed by \insertCite{fumanner;textual}{MBNMAtime}
#'
#' \deqn{{E_{max}}\times\frac{(1-exp(-{rate}\times{x}))}{(1-exp(-{rate}\times{max(x)}))}}
#'
#' @param pool.emax Pooling for exponential Emax parameter. Can take `"rel"` or `"abs"` (see details).
#' @param method.emax Method for synthesis of exponential Emax parameter. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param pool.rate Pooling for parameter controlling rate of onset. Default is `NULL` which avoids including
#' this parameter (i.e. fixes it to 1 for all treatments). Can take `"rel"` or `"abs"` (see details).
#' @param method.rate Method for synthesis of parameter controlling rate of onset. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @inheritParams temax
#'
#' @return An object of `class("timefun")`
#'
#' @section Time-course parameters:
#' Time-course parameters in the model must be specified using a `pool` and a `method` prefix.
#'
#' `pool` is used to define the approach used for pooling of a given time-course parameter and
#' can take any of:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Indicates that \emph{relative} effects should be pooled for this time-course parameter. Relative effects preserve randomisation within included studies, are likely to vary less between studies (only due to effect modification), and allow for testing of consistency between direct and indirect evidence. Pooling follows the general approach for Network Meta-Analysis proposed by \insertCite{lu2004;textual}{MBNMAtime}. |
#' | `"abs"` | Indicates that study arms should be pooled across the whole network for this time-course parameter  *independently of assigned treatment* to estimate an \emph{absolute} effect. This implies estimating a single value across the network for this time-course parameter, and may therefore be making strong assumptions of similarity. |
#'
#'
#' `method` is used to define the model used for meta-analysis for a given time-course parameter
#' and can take any of the following values:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"common"` | Implies that all studies estimate the same true effect (often called a "fixed effect" meta-analysis) |
#' | `"random"` | Implies that all studies estimate a separate true effect, but that each of these true effects vary randomly around a true mean effect. This approach allows for modelling of between-study heterogeneity. |
#' | `numeric()` | Assigned a numeric value, indicating that this time-course parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific time-course parameters (e.g. Hill parameters in Emax functions, power parameters in fractional polynomials) to a single value. |
#'
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' titp(pool.emax="rel", method.emax="random")
#' titp(pool.emax="abs")
#'
#' @export
titp <- function(pool.emax="rel", method.emax="common",
                 pool.rate="rel", method.rate="common", p.expon=FALSE) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(pool.emax, choices=c("rel", "abs"), add=argcheck)
  #checkmate::assertChoice(method.emax, choices=c("common", "random"), add=argcheck)
  checkmate::assertChoice(pool.rate, choices=c("rel", "abs"), null.ok = TRUE, add=argcheck)
  #checkmate::assertChoice(method.onset, choices=c("common", "random"), null.ok = TRUE, add=argcheck)
  checkmate::assertLogical(p.expon, add=argcheck)
  checkmate::reportAssertions(argcheck)

  params <- list(method.emax=method.emax, method.rate=method.rate)
  for (i in seq_along(params)) {
    if (!is.null(params[[i]])) {
      err <- TRUE
      if (length(params[[i]])==1) {
        if (any(c("common", "random") %in% params[[i]])) {
          err <- FALSE
        } else if (is.numeric(params[[i]])) {
          err <- FALSE
        }
      }
      if (err) {
        stop(paste0(names(params)[i], " must take either 'common', 'random' or be assigned a numeric value"))
      }
    }
  }

  # Define time-course function
  if (p.expon==TRUE) {
    fun <- ~ emax * (1 - exp(-exp(rate)*time)) / (1 - exp(-exp(rate)*max(time)))
    jags <- "beta.1 * ((1-exp(-exp(beta.2)*time[i,m])) / (1-exp(-exp(beta.2)*maxtime)))"
    latex <- "\beta_1 * (1-exp(-exp(\beta_2)*x_m)) / (1-exp(-exp(\beta_2)*max(x_m)))"

  } else if (p.expon==FALSE) {
    fun <- ~ emax * (1 - exp(-rate*time)) / (1 - exp(-rate*max(time)))
    jags <- "beta.1 * ((1-exp(-beta.2*time[i,m])) / (1-exp(-beta.2*maxtime)))"
    latex <- "\beta_1 * (1-exp(-\beta_2*x_m)) / (1-exp(-\beta_2*max(x_m)))"
  }


  f <- function(time, beta.1, beta.2) {
    y <- beta.1 * (1-exp(-beta.2*time)) / (1-exp(-beta.2*max(time)))
    return(y)
  }

  if (pool.emax=="rel") {
    jags <- gsub("beta\\.1", "beta.1[i,k]", jags)
  } else if (pool.emax=="abs" & method.emax=="random") {
    jags <- gsub("beta\\.1", "i.beta.1[i,k]", jags)
  }
  if (pool.rate=="rel") {
    jags <- gsub("beta\\.2", "beta.2[i,k]", jags)
  } else if (pool.rate=="abs" & method.rate=="random") {
    jags <- gsub("beta\\.2", "i.beta.2[i,k]", jags)
  }


  # Generate output values
  paramnames <- c("emax", "rate")
  nparam <- length(paramnames)

  apool <- pool.emax
  amethod <- method.emax

  apool <- append(apool, pool.rate)
  amethod <- append(amethod, method.rate)

  names(apool) <- paramnames
  names(amethod) <- paramnames

  bname <- paste0("beta.", 1:nparam)
  names(bname) <- paramnames

  bpool <- paste0("pool.", 1:nparam)
  names(bpool) <- paramnames
  bmethod <- paste0("method.", 1:nparam)
  names(bmethod) <- paramnames

  out <- list(name="itp", fun=fun, f=f, latex=latex,
              params=paramnames, nparam=nparam, jags=jags,
              apool=apool, amethod=amethod, bname=bname,
              bpool=bpool, bmethod=bmethod, p.expon=p.expon)

  class(out) <- "timefun"

  if (p.expon==TRUE) {
    message("'rate' parameters are on exponential scale to ensure they take positive values on the natural scale")

  } else if (p.expon==FALSE) {
    message("'rate' parameters must take positive values.\n Default half-normal prior restricts posterior to positive values.")
  }

  return(out)
}






#' Log-linear (exponential) time-course function
#'
#' \eqn{rate\times{log(x + 1)}}
#'
#' @param pool.rate Pooling for rate parameter. Can take `"rel"` or `"abs"` (see details).
#' @param method.rate Method for synthesis of rate parameter. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#'
#' @return An object of `class("timefun")`
#'
#'
#' @section Time-course parameters:
#' Time-course parameters in the model must be specified using a `pool` and a `method` prefix.
#'
#' `pool` is used to define the approach used for pooling of a given time-course parameter and
#' can take any of:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Indicates that \emph{relative} effects should be pooled for this time-course parameter. Relative effects preserve randomisation within included studies, are likely to vary less between studies (only due to effect modification), and allow for testing of consistency between direct and indirect evidence. Pooling follows the general approach for Network Meta-Analysis proposed by \insertCite{lu2004;textual}{MBNMAtime}. |
#' | `"abs"` | Indicates that study arms should be pooled across the whole network for this time-course parameter  *independently of assigned treatment* to estimate an \emph{absolute} effect. This implies estimating a single value across the network for this time-course parameter, and may therefore be making strong assumptions of similarity. |
#'
#'
#' `method` is used to define the model used for meta-analysis for a given time-course parameter
#' and can take any of the following values:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"common"` | Implies that all studies estimate the same true effect (often called a "fixed effect" meta-analysis) |
#' | `"random"` | Implies that all studies estimate a separate true effect, but that each of these true effects vary randomly around a true mean effect. This approach allows for modelling of between-study heterogeneity. |
#' | `numeric()` | Assigned a numeric value, indicating that this time-course parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific time-course parameters (e.g. Hill parameters in Emax functions, power parameters in fractional polynomials) to a single value. |
#'
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' tloglin(pool.rate="rel", method.rate="random")
#' tloglin(pool.rate="abs")
#'
#' @export
tloglin <- function(pool.rate="rel", method.rate="common") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(pool.rate, choices=c("rel", "abs"), add=argcheck)
  # checkmate::assertChoice(method.rate, choices=c("common", "random"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  params <- list(method.rate=method.rate)
  for (i in seq_along(params)) {
    if (!is.null(params[[i]])) {
      err <- TRUE
      if (length(params[[i]])==1) {
        if (any(c("common", "random") %in% params[[i]])) {
          err <- FALSE
        } else if (is.numeric(params[[i]])) {
          err <- FALSE
        }
      }
      if (err) {
        stop(paste0(names(params)[i], " must take either 'common', 'random' or be assigned a numeric value"))
      }
    }
  }

  # Define time-course function
  fun <- ~ rate * log(time + 1)
  latex <- "\beta_1 * log(x_m + 1)"
  jags <- paste0("beta.1 * log(time[i,m] + 1)")

  if (pool.rate=="rel") {
    jags <- gsub("beta\\.1", "beta.1[i,k]", jags)
  } else if (pool.rate=="abs" & method.rate=="random") {
    jags <- gsub("beta\\.1", "i.beta.1[i,k]", jags)
  }

  f <- function(time, beta.1) {
    y <- beta.1 * log(time + 1)
    return(y)
  }

  # Generate output values
  paramnames <- "rate"
  nparam <- 1

  apool <- pool.rate
  names(apool) <- paramnames
  amethod <- method.rate
  names(amethod) <- paramnames
  bname <- paste0("beta.", 1:nparam)
  names(bname) <- paramnames

  bpool <- paste0("pool.", 1:nparam)
  names(bpool) <- paramnames
  bmethod <- paste0("method.", 1:nparam)
  names(bmethod) <- paramnames

  out <- list(name="loglin", fun=fun, f=f, latex=latex,
              params=paramnames, nparam=nparam, jags=jags,
              apool=apool, amethod=amethod, bname=bname,
              bpool=bpool, bmethod=bmethod)

  class(out) <- "timefun"

  return(out)

}





#' Emax time-course function
#'
#' ** For version 0.2.3: to ensure positive posterior values, et50 and hill parameters are now
#' modeled on the natural scale using a half-normal prior rather than a symmetrical prior
#' on the exponential scale to improve model stability **
#'
#' @param pool.emax Pooling for Emax parameter. Can take `"rel"` or `"abs"` (see details).
#' @param method.emax Method for synthesis of Emax parameter. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param pool.et50 Pooling for ET50 parameter. Can take `"rel"` or `"abs"` (see details).
#' @param method.et50 Method for synthesis of ET50 parameter. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param pool.hill Pooling for Hill parameter. Can take `"rel"` or `"abs"` (see details).
#' @param method.hill Method for synthesis of Hill parameter. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param p.expon Should parameters that can only take positive values be modeled on the exponential scale (`TRUE`)
#' or should they be assigned a prior that restricts the posterior to positive values (`FALSE`)
#'
#' @return An object of `class("timefun")`
#'
#' @details
#'
#' * Emax represents the maximum response.
#' * ET50 represents the time at which 50% of the maximum response is achieved. This can only take
#' positive values and so is modeled on the exponential scale and assigned a symmetrical normal prior
#' Alternatively it can be assigned a normal prior truncated at zero (half-normal) (this
#' will be the default in MBNMAtime version >=0.2.3).
#' * Hill is the Hill parameter, which allows for a sigmoidal function. This can only take
#' positive values and so is modeled on the exponential scale and assigned a symmetrical normal prior
#' Alternatively it can be assigned a normal prior truncated at zero (half-normal) (this
#' will be the default in MBNMAtime version >=0.2.3).
#'
#' Without Hill parameter:
#' \deqn{\frac{E_{max}\times{x}}{ET_{50}+x}}
#'
#' With Hill parameter:
#' \deqn{\frac{E_{max}\times{x^{hill}}}{ET_{50}\times{hill}+x^{hill}}}
#'
#'
#' @section Time-course parameters:
#' Time-course parameters in the model must be specified using a `pool` and a `method` prefix.
#'
#' `pool` is used to define the approach used for pooling of a given time-course parameter and
#' can take any of:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Indicates that \emph{relative} effects should be pooled for this time-course parameter. Relative effects preserve randomisation within included studies, are likely to vary less between studies (only due to effect modification), and allow for testing of consistency between direct and indirect evidence. Pooling follows the general approach for Network Meta-Analysis proposed by \insertCite{lu2004;textual}{MBNMAtime}. |
#' | `"abs"` | Indicates that study arms should be pooled across the whole network for this time-course parameter  *independently of assigned treatment* to estimate an \emph{absolute} effect. This implies estimating a single value across the network for this time-course parameter, and may therefore be making strong assumptions of similarity. |
#'
#'
#' `method` is used to define the model used for meta-analysis for a given time-course parameter
#' and can take any of the following values:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"common"` | Implies that all studies estimate the same true effect (often called a "fixed effect" meta-analysis) |
#' | `"random"` | Implies that all studies estimate a separate true effect, but that each of these true effects vary randomly around a true mean effect. This approach allows for modelling of between-study heterogeneity. |
#' | `numeric()` | Assigned a numeric value, indicating that this time-course parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific time-course parameters (e.g. Hill parameters in Emax functions, power parameters in fractional polynomials) to a single value. |
#'
#'
#' When relative effects are modelled on more than one time-course parameter,
#' correlation between them is automatically estimated using a vague inverse-Wishart prior.
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' # Model without a Hill parameter
#' temax(pool.emax="rel", method.emax="random", pool.et50="abs", method.et50="common")
#'
#' # Model including a Hill parameter and defaults for Emax and ET50 parameters
#' temax(pool.hill="abs", method.hill="common")
#'
#' @export
temax <- function(pool.emax="rel", method.emax="common", pool.et50="rel", method.et50="common",
                 pool.hill=NULL, method.hill=NULL, p.expon=FALSE) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(pool.emax, choices=c("rel", "abs"), add=argcheck)
  # checkmate::assertChoice(method.emax, choices=c("common", "random"), add=argcheck)
  checkmate::assertChoice(pool.et50, choices=c("rel", "abs"), add=argcheck)
  # checkmate::assertChoice(method.et50, choices=c("common", "random"), add=argcheck)
  checkmate::assertChoice(pool.hill, choices=c("rel", "abs"), null.ok = TRUE, add=argcheck)
  # checkmate::assertChoice(method.hill, choices=c("common", "random"), null.ok = TRUE, add=argcheck)
  checkmate::assertLogical(p.expon, add=argcheck)
  checkmate::reportAssertions(argcheck)

  params <- list(method.emax=method.emax, method.et50=method.et50, method.hill=method.hill)
  for (i in seq_along(params)) {
    if (!is.null(params[[i]])) {
      err <- TRUE
      if (length(params[[i]])==1) {
        if (any(c("common", "random") %in% params[[i]])) {
          err <- FALSE
        } else if (is.numeric(params[[i]])) {
          err <- FALSE
        }
      }
      if (err) {
        stop(paste0(names(params)[i], " must take either 'common', 'random' or be assigned a numeric value"))
      }
    }
  }

  if (is.null(method.hill)) {
    ehill <- FALSE
  } else {
    ehill <- TRUE
    #pool.hill <- "abs"
  }

  if (p.expon==TRUE) {
    if (ehill) {
      fun <- ~ (emax * (time ^ exp(hill))) / ((exp(et50) ^ exp(hill)) + (time ^ exp(hill)))
      jags <- "(beta.1 * (time[i,m] ^ exp(beta.3))) / ((exp(beta.2) ^ exp(beta.3)) + (time[i,m] ^ exp(beta.3)))"
      latex <- "$\\frac{\\beta_1 \\times x_m^{e^\\beta_3}}{{e^\\beta_2}^{e^\\beta_3} + x_m^{e^\\beta_3}}$"
    } else {
      fun <- ~ (emax * time) / (exp(et50) + time)
      jags <- "(beta.1 * time[i,m]) / (exp(beta.2) + time[i,m])"
      latex <- "$\\frac{\\beta_1 \\times x_m}{{e^\\beta_2} + x_m}$"
    }

  } else {
    if (ehill) {
      fun <- ~ (emax * (time ^ hill)) / ((et50 ^ hill) + (time ^ hill))
      jags <- "(beta.1 * (abs(time[i,m]) ^ beta.3)) / ((abs(beta.2) ^ beta.3) + (abs(time[i,m]) ^ beta.3))"
      latex <- "$\\frac{\\beta_1 \\times x_m^{\\beta_3}}{{\\beta_2}^{\\beta_3} + x_m^{\\beta_3}}$"
    } else {
      fun <- ~ (emax * time) / (et50 + time)
      jags <- "(beta.1 * time[i,m]) / (beta.2 + time[i,m])"
      latex <- "$\\frac{\\beta_1 \\times x_m}{{\\beta_2} + x_m}$"
    }
  }



  if (pool.emax=="rel") {
    jags <- gsub("beta\\.1", "beta.1[i,k]", jags)
  } else if (pool.emax=="abs" & method.emax=="random") {
    jags <- gsub("beta\\.1", "i.beta.1[i,k]", jags)
  }

  if (pool.et50=="rel") {
    jags <- gsub("beta\\.2", "beta.2[i,k]", jags)
  } else if (pool.et50=="abs" & method.et50=="random") {
    jags <- gsub("beta\\.2", "i.beta.2[i,k]", jags)
  }

  if (ehill) {
    if (pool.hill=="rel") {
      jags <- gsub("beta\\.3", "beta.3[i,k]", jags)
    } else if (pool.hill=="abs" & method.hill=="random") {
      jags <- gsub("beta\\.3", "i.beta.3[i,k]", jags)
    }
  }


  f <- function(time, beta.1, beta.2, beta.3, p.expon) {

    if (p.expon==TRUE) {
      y <- (beta.1 * (time ^ exp(beta.3)) ) / ((exp(beta.2) ^ exp(beta.3)) + (time ^ exp(beta.3)))
    } else {
      y <- (beta.1 * (time ^ beta.3) ) / ((exp(beta.2) ^ beta.3) + (time ^ beta.3))
    }
    return(y)
  }


  # Generate output values
  paramnames <- c("emax", "et50")
  nparam <- 2

  apool <- c(pool.emax, pool.et50)
  amethod <- c(method.emax, method.et50)

  if (ehill) {
    paramnames <- append(paramnames, "hill")
    nparam <- 3
    apool <- append(apool, pool.hill)
    amethod <- append(amethod, method.hill)
  }
  bname <- paste0("beta.", 1:nparam)

  names(apool) <- paramnames
  names(amethod) <- paramnames
  names(bname) <- paramnames

  bpool <- paste0("pool.", 1:nparam)
  names(bpool) <- paramnames
  bmethod <- paste0("method.", 1:nparam)
  names(bmethod) <- paramnames

  out <- list(name="emax", fun=fun, f=f,
              latex=latex, params=paramnames, nparam=nparam, jags=jags,
              apool=apool, amethod=amethod, bname=bname,
              bpool=bpool, bmethod=bmethod, p.expon=p.expon)
  class(out) <- "timefun"

  if (p.expon==TRUE) {
    message("'et50' parameters are on exponential scale to ensure they take positive values on the natural scale")

    if (ehill) {
      message("'hill' parameters are on exponential scale to ensure they take positive values on the natural scale")
    }
  } else if (p.expon==FALSE) {
    message("'et50' parameters must take positive values.\n Default half-normal prior restricts posterior to positive values.")

    if (ehill) {
      message("'hill' parameters must take positive values.\n Default half-normal prior restricts posterior to positive values.")
    }
  }


  return(out)
}



#' Polynomial time-course function
#'
#' @param degree The degree of the polynomial - e.g. `degree=1` for linear, `degree=2` for quadratic, `degree=3` for cubic.
#' @param pool.1 Pooling for the 1st polynomial coefficient. Can take `"rel"` or `"abs"` (see details).
#' @param method.1 Method for synthesis of the 1st polynomial coefficient.Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param pool.2 Pooling for the 2nd polynomial coefficient. Can take `"rel"` or `"abs"` (see details).
#' @param method.2 Method for synthesis of the 2nd polynomial coefficient. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param pool.3 Pooling for the 3rd polynomial coefficient. Can take `"rel"` or `"abs"` (see details).
#' @param method.3 Method for synthesis of the 3rd polynomial coefficient. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param pool.4 Pooling for the 4th polynomial coefficient. Can take `"rel"` or `"abs"` (see details).
#' @param method.4 Method for synthesis of the 4th polynomial coefficient. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#'
#' @return An object of `class("timefun")`
#'
#' @details
#' * \eqn{\beta_1} represents the 1st coefficient.
#' * \eqn{\beta_2} represents the 2nd coefficient.
#' * \eqn{\beta_3} represents the 3rd coefficient.
#' * \eqn{\beta_4} represents the 4th coefficient.
#'
#' Linear model:
#' \deqn{\beta_1{x}}
#'
#' Quadratic model:
#' \deqn{\beta_1{x} + \beta_2{x^2}}
#'
#' Cubic model:
#' \deqn{\beta_1{x} + \beta_2{x^2} + \beta_3{x^3}}
#'
#' Quartic model:
#' \deqn{\beta_1{x} + \beta_2{x^2} + \beta_3{x^3} + \beta_4{x^4}}
#'
#'
#'
#' @section Time-course parameters:
#' Time-course parameters in the model must be specified using a `pool` and a `method` prefix.
#'
#' `pool` is used to define the approach used for pooling of a given time-course parameter and
#' can take any of:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Indicates that \emph{relative} effects should be pooled for this time-course parameter. Relative effects preserve randomisation within included studies, are likely to vary less between studies (only due to effect modification), and allow for testing of consistency between direct and indirect evidence. Pooling follows the general approach for Network Meta-Analysis proposed by \insertCite{lu2004;textual}{MBNMAtime}. |
#' | `"abs"` | Indicates that study arms should be pooled across the whole network for this time-course parameter  *independently of assigned treatment* to estimate an \emph{absolute} effect. This implies estimating a single value across the network for this time-course parameter, and may therefore be making strong assumptions of similarity. |
#'
#'
#' `method` is used to define the model used for meta-analysis for a given time-course parameter
#' and can take any of the following values:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"common"` | Implies that all studies estimate the same true effect (often called a "fixed effect" meta-analysis) |
#' | `"random"` | Implies that all studies estimate a separate true effect, but that each of these true effects vary randomly around a true mean effect. This approach allows for modelling of between-study heterogeneity. |
#' | `numeric()` | Assigned a numeric value, indicating that this time-course parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific time-course parameters (e.g. Hill parameters in Emax functions, power parameters in fractional polynomials) to a single value. |
#'
#'
#' When relative effects are modelled on more than one time-course parameter,
#' correlation between them is automatically estimated using a vague inverse-Wishart prior.
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' # Linear model with random effects
#' tpoly(pool.1="rel", method.1="random")
#'
#' # Quadratic model with a single absolute parameter estimated for the 2nd coefficient
#' tpoly(pool.1="rel", method.1="common", pool.2="abs", method.2="random")
#'
#' @export
tpoly <- function(degree=1, pool.1="rel", method.1="common", pool.2="rel", method.2="common",
                  pool.3="rel", method.3="common", pool.4="rel", method.4="common") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertIntegerish(degree, lower=1, upper = 4, add=argcheck)
  for (i in 1:4) {
    checkmate::assertChoice(get(paste0("pool.", i)), choices=c("rel", "abs"), add=argcheck)
    # checkmate::assertChoice(get(paste0("method.", i)), choices=c("common", "random"), add=argcheck)
  }
  checkmate::reportAssertions(argcheck)

  params <- list(method.1=method.1, method.2=method.2, method.3=method.3, method.4=method.4)
  for (i in seq_along(params)) {
    err <- TRUE
    if (length(params[[i]])==1) {
      if (any(c("common", "random") %in% params[[i]])) {
        err <- FALSE
      } else if (is.numeric(params[[i]])) {
        err <- FALSE
      }
    }
    if (err) {
      stop(paste0(names(params)[i], " must take either 'common', 'random' or be assigned a numeric value"))
    }
  }

  # Define time-course function
  fun <- "beta.1 * time"
  latex <- "beta_1 * x_m"
  for (i in 1:3) {
    if (degree>i) {
      fun <- paste0(fun, " + beta.", i+1, " * (time^", i+1, ")")
      latex <- paste0(latex, " + beta_", i+1, " * time^", i+1)
    }
  }
  jags <- gsub("time", "time[i,m]", fun)
  fun <- stats::as.formula(paste0("~", fun))


  for (i in 1:degree) {
    if (get(paste0("pool.",i))=="rel") {
      jags <- gsub(paste0("beta\\.", i), paste0("beta.", i, "[i,k]"), jags)
    } else if (get(paste0("pool.",i))=="abs" & get(paste0("method.",i))=="random") {
      jags <- gsub(paste0("beta\\.", i), paste0("i.beta.", i, "[i,k]"), jags)
    }
  }


  # Generate output values
  paramnames <- paste0("beta.", 1:degree)
  nparam <- degree

  apool <- vector()
  amethod <- vector()
  for (i in 1:nparam) {
    apool <- append(apool, get(paste0("pool.",i)))
    amethod <- append(amethod, get(paste0("method.",i)))
  }
  bname <- paste0("beta.", 1:nparam)

  names(apool) <- paramnames
  names(amethod) <- paramnames
  names(bname) <- paramnames

  bpool <- paste0("pool.", 1:nparam)
  bmethod <- paste0("method.", 1:nparam)
  names(bpool) <- paramnames
  names(bmethod) <- paramnames

  out <- list(name="poly", fun=fun, latex=latex, params=paramnames, nparam=nparam, jags=jags,
              apool=apool, amethod=amethod, bname=bname,
              bpool=bpool, bmethod=bmethod)
  class(out) <- "timefun"

  return(out)
}



#' Fractional polynomial time-course function
#'
#' As first described for use in Network Meta-Analysis by \insertCite{jansen2015;textual}{MBNMAtime}.
#'
#' @param degree The degree of the fractional polynomial as defined in  \insertCite{royston1994;textual}{MBNMAtime}
#' @param pool.1 Pooling for the 1st fractional polynomial coefficient. Can take `"rel"` or `"abs"` (see details).
#' @param method.1 Method for synthesis of the 1st fractional polynomial coefficient. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param pool.2 Pooling for the 2nd fractional polynomial coefficient. Can take `"rel"` or `"abs"` (see details).
#' @param method.2 Method for synthesis of the 2nd fractional polynomial coefficient. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param method.power1 Value for the 1st fractional polynomial power. Must take any numeric value in the set `-2, -1, -0.5, 0, 0.5, 1, 2, 3`.
#'   `pool` for this parameter is set to `"abs"`.
#' @param method.power2 Value for the 2nd fractional polynomial power. Must take any numeric value in the set `-2, -1, -0.5, 0, 0.5, 1, 2, 3`.
#'   `pool` for this parameter is set to `"abs"`.
#'
#' @return An object of `class("timefun")`
#'
#' @details
#' * \eqn{\beta_1} represents the 1st coefficient.
#' * \eqn{\beta_2} represents the 2nd coefficient.
#' * \eqn{p_1} represents the 1st power
#' * \eqn{p_2} represents the 2nd power
#'
#' For a polynomial of `degree=1`:
#' \deqn{{\beta_1}x^{p_1}}
#'
#' For a polynomial of `degree=2`:
#' \deqn{{\beta_1}x^{p_1}+{\beta_2}x^{p_2}}
#'
#' \eqn{x^{(p)}} is a regular power except where \eqn{p=0}, where \eqn{x^{(0)}=ln(x)}.
#' If a fractional polynomial power \eqn{p_m} repeats within the function it is multiplied by another \eqn{ln(x)}.
#'
#'
#' @section Time-course parameters:
#' Time-course parameters in the model must be specified using a `pool` and a `method` prefix.
#'
#' `pool` is used to define the approach used for pooling of a given time-course parameter and
#' can take any of:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Indicates that \emph{relative} effects should be pooled for this time-course parameter. Relative effects preserve randomisation within included studies, are likely to vary less between studies (only due to effect modification), and allow for testing of consistency between direct and indirect evidence. Pooling follows the general approach for Network Meta-Analysis proposed by \insertCite{lu2004;textual}{MBNMAtime}. |
#' | `"abs"` | Indicates that study arms should be pooled across the whole network for this time-course parameter  *independently of assigned treatment* to estimate an \emph{absolute} effect. This implies estimating a single value across the network for this time-course parameter, and may therefore be making strong assumptions of similarity. |
#'
#'
#' `method` is used to define the model used for meta-analysis for a given time-course parameter
#' and can take any of the following values:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"common"` | Implies that all studies estimate the same true effect (often called a "fixed effect" meta-analysis) |
#' | `"random"` | Implies that all studies estimate a separate true effect, but that each of these true effects vary randomly around a true mean effect. This approach allows for modelling of between-study heterogeneity. |
#' | `numeric()` | Assigned a numeric value, indicating that this time-course parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific time-course parameters (e.g. Hill parameters in Emax functions, power parameters in fractional polynomials) to a single value. |
#'
#'
#' When relative effects are modelled on more than one time-course parameter,
#' correlation between them is automatically estimated using a vague inverse-Wishart prior.
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' # 1st order fractional polynomial with random effects
#' tfpoly(pool.1="rel", method.1="random")
#'
#' # 2nd order fractional polynomial
#' # with a single absolute parameter estimated for the 2nd coefficient
#' # 1st power equal to zero
#' tfpoly(degree=2, pool.1="rel", method.1="common",
#'   pool.2="abs", method.2="random",
#'   method.power1=0)
#'
#' @export
tfpoly <- function(degree=1, pool.1="rel", method.1="common", pool.2="rel", method.2="common",
                   method.power1=0, method.power2=0) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertIntegerish(degree, len=1, lower=1, upper = 2, add=argcheck)
  for (i in 1:2) {
    checkmate::assertChoice(get(paste0("pool.", i)), choices=c("rel", "abs"), add=argcheck)
    # checkmate::assertChoice(get(paste0("method.", i)), choices=c("common", "random"), add=argcheck)
    checkmate::assertChoice(get(paste0("method.power", i)), choices=c(-2,-1,-0.5,0,0.5,1,2,3), add=argcheck)
  }
  checkmate::reportAssertions(argcheck)

  params <- list(method.1=method.1, method.2=method.2, method.power1=method.power1, method.power2=method.power2)
  for (i in seq_along(params)) {
    err <- TRUE
    if (length(params[[i]])==1) {
      if (any(c("common", "random") %in% params[[i]])) {
        err <- FALSE
      } else if (is.numeric(params[[i]])) {
        err <- FALSE
      }
    }
    if (err) {
      stop(paste0(names(params)[i], " must take either 'common', 'random' or be assigned a numeric value"))
    }
  }


  pool.power1 <- "abs"

  # Define time-course function
  if (degree==1) {
    fun <- ~ beta.1 * ifelse(time>0, ifelse(beta.2==0, log(time), time^beta.2), 0)
    jags <- "beta.1 * ifelse(time[i,m]>0, ifelse(beta.2==0, log(time[i,m]), time[i,m]^beta.2), 0)"
    latex <- "TO BE WRITTEN"
  } else if (degree==2) {
    pool.power2 <- "abs"
    fun <- ~ beta.1 * ifelse(time>0, ifelse(beta.3==0, log(time), time^beta.3), 0) + (beta.2 * ifelse(beta.4==beta.3, ifelse(time>0, ifelse(beta.4==0, log(time)^2, (time^beta.4) * log(time)), 0), ifelse(time>0, ifelse(beta.4==0, log(time), time^beta.4), 0)))
    jags <- "beta.1 * ifelse(time[i,m]>0, ifelse(beta.3==0, log(time[i,m]), time[i,m]^beta.3), 0) + (beta.2 * ifelse(beta.4==beta.3, ifelse(time[i,m]>0, ifelse(beta.4==0, log(time[i,m])^2, (time[i,m]^beta.4) * log(time[i,m])), 0), ifelse(time[i,m]>0, ifelse(beta.4==0, log(time[i,m]), time[i,m]^beta.4), 0)))"
    latex <- "TO BE WRITTEN"
  }

  # Set parameters
  if (pool.1=="rel") {
    jags <- gsub("beta\\.1", "beta.1[i,k]", jags)
  } else if (pool.1=="abs" & method.1=="random") {
    jags <- gsub("beta\\.1", "i.beta.1[i,k]", jags)
  }

  if (degree==1) {
    if (method.power1=="random") {
      jags <- gsub("beta\\.2", "i.beta.2[i,k]", jags)
    }
  } else if (degree==2) {
    if (pool.2=="rel") {
      jags <- gsub("beta\\.2", "beta.2[i,k]", jags)
    } else if (pool.2=="abs" & method.2=="random") {
      jags <- gsub("beta\\.2", "i.beta.2[i,k]", jags)
    }
    if (method.power1=="random") {
      jags <- gsub("beta\\.3", "i.beta.3[i,k]", jags)
    }
    if (method.power2=="random") {
      jags <- gsub("beta\\.4", "i.beta.4[i,k]", jags)
    }
  }


  # Write function
  f1 <- function(time, beta.1, beta.2) {
    if (time>0) {
      if (beta.2==0) {
        y <- log(time)
      } else {
        y <- time^beta.2
      }
    } else {
      y <- 0
    }
    return(beta.1 * y)
  }

  f2 <- function(time, beta.1, beta.2, beta.3, beta.4) {
    if (time>0) {
      if (beta.3==0) {
        y1 <- log(time)
      } else {
        y1 <- time^beta.3
      }
    } else {
      y1 <- 0
    }
    y1 <- beta.1 * y1

    if (beta.4==beta.3) {
      if (time>0) {
        if (beta.4==0) {
          y2 <- log(time)^2
        } else {
          y2 <- time^beta.4 * log(time)
        }
      } else {
        if (time >0) {
          if (beta.4==0) {
            y2 <- log(time)
          } else {
            y2 <- time^beta.4
          }
        }
      }
    } else {
      y2 <- 0
    }
    return(y1 + y2)
  }

  if (degree==1) {
    f <- f1
  } else if (degree==2) {
    f <- f2
  }


  # Generate output values
  paramnames <- c(paste0("beta.", 1:degree), paste0("power", 1:degree))
  nparam <- degree*2

  apool <- vector()
  amethod <- vector()
  for (i in 1:degree) {
    apool <- append(apool, get(paste0("pool.",i)))
    amethod <- append(amethod, get(paste0("method.",i)))
  }
  for (i in 1:degree) {
    apool <- append(apool, get(paste0("pool.power",i)))
    amethod <- append(amethod, get(paste0("method.power",i)))
  }
  bname <- paste0("beta.", 1:nparam)

  #names(apool) <- paramnames[1:degree]
  names(apool) <- paramnames
  names(amethod) <- paramnames
  names(bname) <- paramnames

  #bpool <- paste0("pool.", 1:degree)
  bpool <- paste0("pool.", 1:nparam)
  bmethod <- paste0("method.", 1:nparam)
  #names(bpool) <- paramnames[1:degree]
  names(bpool) <- paramnames[1:nparam]
  names(bmethod) <- paramnames

  out <- list(name="fpoly", fun=fun, f=f,
              latex=latex, params=paramnames, nparam=nparam, jags=jags,
              apool=apool, amethod=amethod, bname=bname,
              bpool=bpool, bmethod=bmethod)
  class(out) <- "timefun"

  return(out)

}













#' Spline time-course functions
#'
#' Used to fit B-splines, natural cubic splines, and
#' piecewise linear splines\insertCite{perperoglu2019}{MBNMAtime}. Note that
#' B-splines with `degree=1` and linear splines are equivalent in fit, but are
#' parameterised differently which can allow different informative prior specification.
#'
#' @param type The type of spline. Can take `"bs"` (\href{https://mathworld.wolfram.com/B-Spline.html}{B-spline}),
#'   `"ns"` (\href{https://mathworld.wolfram.com/CubicSpline.html}{natural cubic spline}) or `"ls"` (piecewise linear spline)
#' @param degree The degree of the piecewise B-spline polynomial - e.g. `degree=1` for linear, `degree=2` for quadratic, `degree=3` for cubic.
#' @param pool.1 Pooling for the 1st coefficient. Can take `"rel"` or `"abs"` (see details).
#' @param method.1 Method for synthesis of the 1st coefficient. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param pool.2 Pooling for the 2nd coefficient. Can take `"rel"` or `"abs"` (see details).
#' @param method.2 Method for synthesis of the 2nd coefficient. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param pool.3 Pooling for the 3rd coefficient. Can take `"rel"` or `"abs"` (see details).
#' @param method.3 Method for synthesis of the 3rd coefficient. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param pool.4 Pooling for the 4th coefficient. Can take `"rel"` or `"abs"` (see details).
#' @param method.4 Method for synthesis of the 4th coefficient. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param pool.5 Pooling for the 5th coefficient. Can take `"rel"` or `"abs"` (see details).
#' @param method.5 Method for synthesis of the 5th coefficient. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param pool.6 Pooling for the 6th coefficient. Can take `"rel"` or `"abs"` (see details).
#' @param method.6 Method for synthesis of the 6th coefficient. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @inheritParams genspline
#'
#' @return An object of `class("timefun")`
#'
#' @section Time-course parameters:
#' Time-course parameters in the model must be specified using a `pool` and a `method` prefix.
#'
#' `pool` is used to define the approach used for pooling of a given time-course parameter and
#' can take any of:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Indicates that \emph{relative} effects should be pooled for this time-course parameter. Relative effects preserve randomisation within included studies, are likely to vary less between studies (only due to effect modification), and allow for testing of consistency between direct and indirect evidence. Pooling follows the general approach for Network Meta-Analysis proposed by \insertCite{lu2004;textual}{MBNMAtime}. |
#' | `"abs"` | Indicates that study arms should be pooled across the whole network for this time-course parameter  *independently of assigned treatment* to estimate an \emph{absolute} effect. This implies estimating a single value across the network for this time-course parameter, and may therefore be making strong assumptions of similarity. |
#'
#'
#' `method` is used to define the model used for meta-analysis for a given time-course parameter
#' and can take any of the following values:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"common"` | Implies that all studies estimate the same true effect (often called a "fixed effect" meta-analysis) |
#' | `"random"` | Implies that all studies estimate a separate true effect, but that each of these true effects vary randomly around a true mean effect. This approach allows for modelling of between-study heterogeneity. |
#' | `numeric()` | Assigned a numeric value, indicating that this time-course parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific time-course parameters (e.g. Hill parameters in Emax functions, power parameters in fractional polynomials) to a single value. |
#'
#'
#' When relative effects are modelled on more than one time-course parameter,
#' correlation between them is automatically estimated using a vague inverse-Wishart prior.
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' # Second order B spline with 2 equally spaced knots and random effects on
#' #the 2nd coefficient
#' tspline(type="bs", nknots=2, degree=2,
#'   pool.1="rel", method.1="common",
#'   pool.2="rel", method.2="random")
#'
#' # Piecewise linear spline with knots at times of 5 and 10
#' # Single parameter independent of treatment estimated for 1st coefficient
#' #with random effects
#' tspline(type="ls", knots=c(5,10),
#'   pool.1="abs", method.1="random",
#'   pool.2="rel", method.2="common")
#'
#' @export
tspline <- function(type="bs", knots=NULL, nknots=1, degree=1, pool.1="rel", method.1="common",
                      pool.2="rel", method.2="common", pool.3="rel", method.3="common",
                      pool.4="rel", method.4="common", pool.5="rel", method.5="common",
                      pool.6="rel", method.6="common") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(type, choices=c("bs", "ns", "ls"), add=argcheck)
  checkmate::assertNumeric(nknots, null.ok=FALSE, len=1, lower=1, add=argcheck)
  checkmate::assertNumeric(knots, null.ok=TRUE, lower=0, add=argcheck)
  checkmate::assertIntegerish(degree, add=argcheck)

  for (i in 1:6) {
    checkmate::assertChoice(get(paste0("pool.", i)), choices=c("rel", "abs"), add=argcheck)
    # checkmate::assertChoice(get(paste0("method.", i)), choices=c("common", "random"), add=argcheck)
  }
  checkmate::reportAssertions(argcheck)

  params <- list(method.1=method.1, method.2=method.2, method.3=method.3, method.4=method.4)
  for (i in seq_along(params)) {
    err <- TRUE
    if (length(params[[i]])==1) {
      if (any(c("common", "random") %in% params[[i]])) {
        err <- FALSE
      } else if (is.numeric(params[[i]])) {
        err <- FALSE
      }
    }
    if (err) {
      stop(paste0(names(params)[i], " must take either 'common', 'random' or be assigned a numeric value"))
    }
  }

  # Check knots and degrees
  x <- c(0:100)
  if (!is.null(knots)) {nknots <- length(knots)}
  x <- genspline(x, spline=type, knots = knots, nknots=nknots, degree=degree)

  nparam <- ncol(x)

  # Define time-course function
  base <- "beta.1 * spline.1"
  basetex <- "\beta_1 * X[m,1]"
  jags <- base
  latex <- basetex
  if (nparam>1) {
    for (i in 2:(nparam)) {
      temp <- gsub("1", i, base)
      jags <- paste(jags, "+", temp)

      temptex <- gsub("1", i, basetex)
      latex <- paste(latex, "+", temptex)
    }
  }
  fun <- stats::as.formula(paste("~", jags))
  jags <- gsub("(spline)\\.([0-9])", "\\1[i,m,\\2]", jags)


  # Define parameters
  for (i in 1:(nparam)) {
    if (get(paste0("pool.",i))=="rel") {
      jags <- gsub(paste0("beta\\.", i), paste0("beta.", i, "[i,k]"), jags)
    } else if (get(paste0("pool.",i))=="abs" & get(paste0("method.",i))=="random") {
      jags <- gsub(paste0("beta\\.", i), paste0("i.beta.", i, "[i,k]"), jags)
    }
  }


  # Generate output values
  paramnames <- paste0("beta.", 1:nparam)

  apool <- vector()
  amethod <- vector()
  for (i in 1:nparam) {
    apool <- append(apool, get(paste0("pool.",i)))
    amethod <- append(amethod, get(paste0("method.",i)))
  }
  bname <- paste0("beta.", 1:nparam)

  names(apool) <- paramnames
  names(amethod) <- paramnames
  names(bname) <- paramnames

  bpool <- paste0("pool.", 1:nparam)
  bmethod <- paste0("method.", 1:nparam)
  names(bpool) <- paramnames
  names(bmethod) <- paramnames

  out <- list(name=type, fun=fun, latex=latex, params=paramnames,
              nparam=nparam, knots=knots, nknots=nknots, degree=degree, jags=jags,
              apool=apool, amethod=amethod, bname=bname,
              bpool=bpool, bmethod=bmethod)
  class(out) <- "timefun"

  return(out)

}






#' User-defined time-course function
#'
#' @param fun A formula specifying any relationship including `time` and
#'   one/several of: `beta.1`, `beta.2`, `beta.3`, `beta.4`.
#' @param pool.1 Pooling for `beta.1`. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param method.1 Method for synthesis of `beta.1`. Can take `"common` or `"random"` (see details).
#' @param pool.2 Pooling for `beta.2`. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param method.2 Method for synthesis of `beta.2. Can take `"common` or `"random"` (see details).
#' @param pool.3 Pooling for `beta.3`. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param method.3 Method for synthesis of `beta.3. Can take `"common` or `"random"` (see details).
#' @param pool.4 Pooling for `beta.4`. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#' @param method.4 Method for synthesis of `beta.4`. Can take `"common`, `"random"`, or be assigned a numeric value (see details).
#'
#' @return An object of `class("timefun")`
#'
#' @section Time-course parameters:
#' Time-course parameters in the model must be specified using a `pool` and a `method` prefix.
#'
#' `pool` is used to define the approach used for pooling of a given time-course parameter and
#' can take any of:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Indicates that \emph{relative} effects should be pooled for this time-course parameter. Relative effects preserve randomisation within included studies, are likely to vary less between studies (only due to effect modification), and allow for testing of consistency between direct and indirect evidence. Pooling follows the general approach for Network Meta-Analysis proposed by \insertCite{lu2004;textual}{MBNMAtime}. |
#' | `"abs"` | Indicates that study arms should be pooled across the whole network for this time-course parameter  *independently of assigned treatment* to estimate an \emph{absolute} effect. This implies estimating a single value across the network for this time-course parameter, and may therefore be making strong assumptions of similarity. |
#'
#'
#' `method` is used to define the model used for meta-analysis for a given time-course parameter
#' and can take any of the following values:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"common"` | Implies that all studies estimate the same true effect (often called a "fixed effect" meta-analysis) |
#' | `"random"` | Implies that all studies estimate a separate true effect, but that each of these true effects vary randomly around a true mean effect. This approach allows for modelling of between-study heterogeneity. |
#' | `numeric()` | Assigned a numeric value, indicating that this time-course parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific time-course parameters (e.g. Hill parameters in Emax functions, power parameters in fractional polynomials) to a single value. |
#'
#'
#' When relative effects are modelled on more than one time-course parameter,
#' correlation between them is automatically estimated using a vague inverse-Wishart prior.
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#'
#' timecourse <- ~ beta.1 * (1/(time+1)) + beta.2 * time^2
#' tuser(fun=timecourse,
#'   pool.1="abs", method.1="common",
#'   pool.2="rel", method.2="common")
#'
#' @export
tuser <- function(fun, pool.1="rel", method.1="common",
                      pool.2="rel", method.2="common", pool.3="rel", method.3="common",
                      pool.4="rel", method.4="common") {


  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertFormula(fun, add=argcheck)
  for (i in 1:4) {
    checkmate::assertChoice(get(paste0("pool.", i)), choices=c("rel", "abs"), add=argcheck)
    # checkmate::assertChoice(get(paste0("method.", i)), choices=c("common", "random"), add=argcheck)
  }
  checkmate::reportAssertions(argcheck)

  params <- list(method.1=method.1, method.2=method.2, method.3=method.3, method.4=method.4)
  for (i in seq_along(params)) {
    err <- TRUE
    if (length(params[[i]])==1) {
      if (any(c("common", "random") %in% params[[i]])) {
        err <- FALSE
      } else if (is.numeric(params[[i]])) {
        err <- FALSE
      }
    }
    if (err) {
      stop(paste0(names(params)[i], " must take either 'common', 'random' or be assigned a numeric value"))
    }
  }

  # Check user function
  user.str <- as.character(fun[2])
  if (grepl("beta\\.2", user.str)==TRUE & grepl("beta\\.1", user.str)==FALSE) {
    stop("'fun' cannot contain beta.2 if beta.1 is not present")
  } else if (grepl("beta\\.3", user.str)==TRUE & grepl("beta\\.2", user.str)==FALSE | grepl("beta\\.1", user.str)==FALSE) {
    stop("'fun' cannot contain beta.3 if beta.2 and beta.1 are not present")
  }
  if (!(grepl("time", user.str))) {
    stop("'fun' must be a function of beta parameters and time")
  }
  if (grepl("alpha", user.str)) {
    stop("The intercept (alpha) should not be included in 'fun'")
  }
  jags <- gsub("time", "time[i,m]", user.str)


  # Get number of parameters
  nparam <- 1
  if (grepl("beta\\.4", user.str)) {
    nparam <- 4
  } else if (grepl("beta\\.3", user.str)) {
    nparam <- 3
  } else if (grepl("beta\\.2", user.str)) {
    nparam <- 2
  }

  # Define parameters
  for (i in 1:nparam) {
    if (get(paste0("pool.",i))=="rel") {
      jags <- gsub(paste0("beta\\.", i), paste0("beta.", i, "[i,k]"), jags)
    } else if (get(paste0("pool.",i))=="abs" & get(paste0("method.",i))=="random") {
      jags <- gsub(paste0("beta\\.", i), paste0("i.beta.", i, "[i,k]"), jags)
    }
  }

  # Generate output values
  paramnames <- paste0("beta.", 1:nparam)

  apool <- vector()
  amethod <- vector()
  for (i in 1:nparam) {
    apool <- append(apool, get(paste0("pool.",i)))
    amethod <- append(amethod, get(paste0("method.",i)))
  }
  bname <- paste0("beta.", 1:nparam)

  names(apool) <- paramnames
  names(amethod) <- paramnames
  names(bname) <- names(bname)

  bpool <- paste0("pool.", 1:nparam)
  bmethod <- paste0("method.", 1:nparam)
  names(bpool) <- paramnames
  names(bmethod) <- paramnames

  out <- list(name="user", fun=fun, latex=NULL, params=paramnames, nparam=nparam, jags=jags,
              apool=apool, amethod=amethod, bname=bname,
              bpool=bpool, bmethod=bmethod)
  class(out) <- "timefun"

  return(out)
}
