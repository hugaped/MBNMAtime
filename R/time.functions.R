# Functions for time-course
# Author: Hugo Pedder
# Date created: 2021-01-04

texp <- function(pool.rate="rel", method.rate="common") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(pool.rate, choices=c("rel", "abs"), add=argcheck)
  checkmate::assertChoice(method.rate, choices=c("common", "random"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Define time-course function
  fun <- ~ rate * (1 - exp(-time))
  latex <- "\beta_1 * (1 - exp(-x_m))"
  jags <- "beta.1 * (1 - exp(- time[i,m]))"

  if (pool.rate=="rel") {
    jags <- gsub("beta\\.1", "beta.1[i,k]", jags)
  } else if (pool.rate=="abs" & method.rate=="random") {
    jags <- gsub("beta\\.1", "i.beta.1[i,k]", jags)
  }

  f <- function(time, beta.1) {
    y <- beta.1 * (1-exp(-time))
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

  out <- list(name="exp", fun=fun, f=f, latex=latex,
              params=paramnames, nparam=nparam, jags=jags,
              apool=apool, amethod=amethod, bname=bname,
              bpool=bpool, bmethod=bmethod)

  class(out) <- "timefun"

  return(out)
}




tloglin <- function(pool.rate="rel", method.rate="common", off=.1) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(pool.rate, choices=c("rel", "abs"), add=argcheck)
  checkmate::assertChoice(method.rate, choices=c("common", "random"), add=argcheck)
  checkmate::assertNumeric(off, len=1, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Define time-course function
  fun <- ~ rate * log(time + off)
  latex <- "\beta_1 * log(x_m + off)"
  jags <- paste0("beta.1 * log(time[i,m] + ", off, ")")

  if (pool.rate=="rel") {
    jags <- gsub("beta\\.1", "beta.1[i,k]", jags)
  } else if (pool.rate=="abs" & method.rate=="random") {
    jags <- gsub("beta\\.1", "i.beta.1[i,k]", jags)
  }

  f <- function(time, beta.1, beta.2, off) {
    y <- beta.1 * log(time + off)
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

  out <- list(name="loglin", fun=fun, f=f, latex=latex, off=off,
              params=paramnames, nparam=nparam, jags=jags,
              apool=apool, amethod=amethod, bname=bname,
              bpool=bpool, bmethod=bmethod)

  class(out) <- "timefun"

  return(out)

}




temax <- function(pool.emax="rel", method.emax="common", pool.et50="rel", method.et50="common",
                 pool.hill=NULL, method.hill=NULL) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(pool.emax, choices=c("rel", "abs"), add=argcheck)
  checkmate::assertChoice(method.emax, choices=c("common", "random"), add=argcheck)
  checkmate::assertChoice(pool.et50, choices=c("rel", "abs"), add=argcheck)
  checkmate::assertChoice(method.et50, choices=c("common", "random"), add=argcheck)
  checkmate::assertChoice(pool.hill, choices=c("rel", "abs"), null.ok = TRUE, add=argcheck)
  checkmate::assertChoice(method.hill, choices=c("common", "random"), null.ok = TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (is.null(method.hill)) {
    ehill <- FALSE
  } else {
    ehill <- TRUE
    pool.hill <- "abs"
  }

  if (ehill) {
    fun <- ~ (emax * (time ^ hill)) / ((exp(et50) ^ hill) + (time ^ hill))
    jags <- "(beta.1 * (time[i,m] ^ exp(beta.3))) / ((exp(beta.2) ^ exp(beta.3)) + (time[i,m] ^ exp(beta.3)))"
    latex <- "$\\frac{\\beta_1 \\times x_m^{e^\\beta_3}}{{e^\\beta_2}^{e^\\beta_3} + x_m^{e^\\beta_3}}$"
    plotmath <- "frac(beta[1] %*% x[m]^e^beta[3], e^beta[2]^e^beta[3] + x[m]^e^beta[3])"
  } else {
    fun <- ~ (emax * time) / (exp(et50) + time)
    jags <- "(beta.1 * time[i,m]) / (exp(beta.2) + time[i,m])"
    latex <- "$\\frac{\\beta_1 \\times x_m}{{e^\\beta_2} + x_m}$"
    plotmath <- "frac(beta[1] %*% x[m], e^beta[2] + x[m])"

    # plot.new()
    # text(0.5,0.5, eval(parse(text=paste0("expression(", plotmath, ")"))))
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
    if (pool.hill=="abs" & method.et50=="random") {
      jags <- gsub("beta\\.3", "i.beta.3[i,k]", jags)
    }
  }


  f <- function(time, beta.1, beta.2, beta.3) {
    y <- (beta.1 * (time ^ beta.3) ) / ((exp(beta.2) ^ beta.3) + (time ^ beta.3))
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
              bpool=bpool, bmethod=bmethod)
  class(out) <- "timefun"

  message("'et50' parameters are on exponential scale to ensure they take positive values on the natural scale")

  if (ehill) {
    message("'hill' parameters are on exponential scale to ensure they take positive values on the natural scale")
  }
  return(out)
}



#'
#' @param degree the degree of the polynomial - e.g. `degree=1` for linear, `degree=2` for quadratic,
tpoly <- function(degree=1, pool.1="rel", method.1="common", pool.2="rel", method.2="common",
                  pool.3="rel", method.3="common", pool.4="rel", method.4="common") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertIntegerish(degree, lower=1, upper = 4, add=argcheck)
  for (i in 1:4) {
    checkmate::assertChoice(get(paste0("pool.", i)), choices=c("rel", "abs"), add=argcheck)
    checkmate::assertChoice(get(paste0("method.", i)), choices=c("common", "random"), add=argcheck)
  }
  checkmate::reportAssertions(argcheck)

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
  fun <- as.formula(paste0("~", fun))


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



#' @param degree the degree of the fractional polynomial
tfpoly <- function(degree=1, pool.1="rel", method.1="common", pool.2="rel", method.2="common",
                   method.power1="common", method.power2="common") {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertIntegerish(degree, len=1, lower=1, upper = 2, add=argcheck)
  for (i in 1:2) {
    checkmate::assertChoice(get(paste0("pool.", i)), choices=c("rel", "abs"), add=argcheck)
    checkmate::assertChoice(get(paste0("method.", i)), choices=c("common", "random"), add=argcheck)
    checkmate::assertChoice(get(paste0("method.power", i)), choices=c("common", "random"), add=argcheck)
  }
  checkmate::reportAssertions(argcheck)

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
      jags <- gsub("beta\\.2", "beta.2[i,k]", jags)
    }
  } else if (degree==2) {
    if (pool.2=="rel") {
      jags <- gsub("beta\\.2", "beta.2[i,k]", jags)
    } else if (pool.2=="abs" & method.2=="random") {
      jags <- gsub("beta\\.2", "i.beta.2[i,k]", jags)
    }
    if (method.power1=="random") {
      jags <- gsub("beta\\.3", "beta.3[i,k]", jags)
    }
    if (method.power2=="random") {
      jags <- gsub("beta\\.4", "beta.4[i,k]", jags)
    }
  }


  # Write function
  if (degree==1) {
    f <- function(time, beta.1, beta.2) {
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
  } else if (degree==2) {
    f <- function(time, beta.1, beta.2, beta.3, beta.4) {
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
    apool <- append(apool, get(paste0("pool.power",degree)))
    amethod <- append(amethod, get(paste0("method.power",degree)))
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


#'
#'
#' @param knots The number/location of spline internal knots. If a single number is given it indicates the number of knots (they will
#'   be equally spaced across the range of time points). If a numeric vector is given it indicates the location of the knots.
tspline <- function(type="bs", knots=1, degree=1, pool.1="rel", method.1="common",
                      pool.2="rel", method.2="common", pool.3="rel", method.3="common",
                      pool.4="rel", method.4="common") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(type, choices=c("rcs", "bs", "ns", "ls"), add=argcheck)
  checkmate::assertNumeric(knots, null.ok=FALSE, add=argcheck)
  checkmate::assertIntegerish(degree, add=argcheck)
  for (i in 1:4) {
    checkmate::assertChoice(get(paste0("pool.", i)), choices=c("rel", "abs"), add=argcheck)
    checkmate::assertChoice(get(paste0("method.", i)), choices=c("common", "random"), add=argcheck)
  }
  checkmate::reportAssertions(argcheck)

  # Check knots and degrees
  x <- c(0:100)
  x <- genspline(x, spline=type, knots = knots, degree=degree)

  nparam <- ncol(x)

  # knoterr <- "Minimum number of `knots` for fun=`rcs` is 3"
  # if (length(knots)==1) {
  #   if (knots>=1) {
  #     if (knots<3 & type=="rcs") {
  #       stop(knoterr)
  #     }
  #     if (type=="rcs") {
  #       nparam <- knots-1
  #     } else if (type=="bs") {
  #       nparam <- knots
  #     }
  #
  #   } else {
  #     if (type=="rcs") {
  #       stop(knoterr)
  #     }
  #   }
  # } else if (length(knots)>1) {
  #   if (length(knots)<3 & type=="rcs") {
  #     stop(knoterr)
  #   }
  #   if (type=="rcs") {
  #     nparam <- length(knots)-1
  #   } else if (type=="bs") {
  #
  #   }
  #
  # }
  #
  # if (length(knots)>1) {
  #   if (!(all(knots<=1 & all(knots>=0)))) {
  #     stop("`knots` specified as quantiles must be between 0 and 1")
  #   }
  # }


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
  fun <- as.formula(paste("~", jags))
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
              nparam=nparam, knots=knots, degree=degree, jags=jags,
              apool=apool, amethod=amethod, bname=bname,
              bpool=bpool, bmethod=bmethod)
  class(out) <- "timefun"

  return(out)

}





# tpiece <- function(knots, pool.1="rel", method.1="common",
#                    pool.2="rel", method.2="common", pool.3="rel", method.3="common",
#                    pool.4="rel", method.4="common") {
#
#   # Run checks
#   argcheck <- checkmate::makeAssertCollection()
#   checkmate::assertNumeric(knots, min.len=1, max.len=3, null.ok = FALSE, add=argcheck)
#   for (i in 1:4) {
#     checkmate::assertChoice(get(paste0("pool.", i)), choices=c("rel", "abs"), add=argcheck)
#     checkmate::assertChoice(get(paste0("method.", i)), choices=c("common", "random"), add=argcheck)
#   }
#   checkmate::reportAssertions(argcheck)
#
#
#   # Define time-course function
#   # Diffiulty of including alpha...
#   latex <- "TO BE WRITTEN"
#
#   "((time < beta.3) * alpha) + ((time < beta.3) * (beta.1 * time)) + ((time >= beta.3) * (alpha + (beta.1 * beta.3))) + ((time >= beta.3) * (beta.2 * time))"
#
#
#
#   # Define parameters
#   for (i in 1:degree) {
#     if (get(paste0("pool."),i)=="rel") {
#       jags <- gsub(paste0("beta\\.", i), paste0("beta.", i, "[i,k]"), jags)
#     } else if (get(paste0("pool."),i)=="abs" & get(paste0("method."),i)=="random") {
#       jags <- gsub(paste0("beta\\.", i), paste0("i.beta.", i, "[i,k]"), jags)
#     }
#   }
#
#
#   # Generate output values
#   paramnames <- paste0("beta.", 1:degree)
#   nparam <- degree
#
#   apool <- vector()
#   amethod <- vector()
#   for (i in 1:nparam) {
#     apool <- append(get(paste0("pool.",i)))
#     amethod <- append(get(paste0("method.",i)))
#   }
#
#   names(apool) <- paramnames
#   names(amethod) <- paramnames
#
#   bpool <- paste0("pool.", 1:nparam)
#   bmethod <- paste0("method.", 1:nparam)
#   names(bpool) <- paramnames
#   names(bmethod) <- paramnames
#
#   out <- list(name=type, fun=fun, latex=latex, params=paramnames, nparam=nparam, jags=jags,
#               apool=apool, amethod=amethod,
#               bpool=bpool, bmethod=bmethod)
#   class(out) <- "timefun"
#
#   return(out)
# }






tuser <- function(fun, pool.1="rel", method.1="common",
                      pool.2="rel", method.2="common", pool.3="rel", method.3="common",
                      pool.4="rel", method.4="common") {


  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertFormula(fun, add=argcheck)
  for (i in 1:4) {
    checkmate::assertChoice(get(paste0("pool.", i)), choices=c("rel", "abs"), add=argcheck)
    checkmate::assertChoice(get(paste0("method.", i)), choices=c("common", "random"), add=argcheck)
  }
  checkmate::reportAssertions(argcheck)

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
    if (get(paste0("pool."),i)=="rel") {
      jags <- gsub(paste0("beta\\.", i), paste0("beta.", i, "[i,k]"), jags)
    } else if (get(paste0("pool."),i)=="abs" & get(paste0("method."),i)=="random") {
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
