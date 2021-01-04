# Functions for time-course
# Author: Hugo Pedder
# Date created: 2021-01-04

linear <- function(pool="rel", method="common") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(pool, choices=c("rel", "abs"), add=argcheck)
  checkmate::assertChoice(method, choices=c("common", "random"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  fun <- ~ beta.1 * time
  jags <- "beta.1 * time[i,m]"

  if (pool=="rel") {
    jags <- gsub("beta\\.1", "beta.1[i,k]", jags)
  } else if (pool=="abs" & method=="random") {
    jags <- gsub("beta\\.1", "s.beta.1[i,k]", jags)
  }

  return(list(name="linear", fun=fun, pool=pool, method=method, jags=jags))
}





emax <- function(pool.emax="rel", method.emax="common", pool.et50="rel", method.et50="common") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(pool.emax, choices=c("rel", "abs"), add=argcheck)
  checkmate::assertChoice(method.emax, choices=c("common", "random"), add=argcheck)
  checkmate::assertChoice(pool.et50, choices=c("rel", "abs"), add=argcheck)
  checkmate::assertChoice(method.et50, choices=c("common", "random"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  fun <- ~ (beta.1 * time) / (exp(beta.2) + time)
  jags <- "(beta.1 * time[i,k]) / (exp(beta.2) + time[i,m])"

  if (pool.emax=="rel") {
    jags <- gsub("beta\\.1", "beta.1[i,k]", jags)
  } else if (pool.emax=="abs" & method.emax=="random") {
    jags <- gsub("beta\\.1", "s.beta.1[i,k]", jags)
  }

  if (pool.et50=="rel") {
    jags <- gsub("beta\\.2", "beta.2[i,k]", jags)
  } else if (pool.et50=="abs" & method.et50=="random") {
    jags <- gsub("beta\\.2", "s.beta.2[i,k]", jags)
  }

  return(list(name="emax", fun=fun, pool.emax=pool.emax, method.emax=method.emax,
              pool.et50=pool.et50, method.et50=method.et50, jags=jags))
}






funspline <- function(type="rcs", knots=3, pool.1="rel", method.1="common",
                      pool.2="rel", method.2="common", pool.3="rel", method.3="common",
                      pool.4="rel", method.4="common") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(type, choices=c("rcs", "bs", "ns"), add=argcheck)
  checkmate::assertNumeric(knots, null.ok=FALSE, add=argcheck)
  for (i in 1:4) {
    checkmate::assertChoice(get(paste0("pool.", i)), choices=c("rel", "abs"), add=argcheck)
    checkmate::assertChoice(get(paste0("method.", i)), choices=c("common", "random"), add=argcheck)
  }
  checkmate::reportAssertions(argcheck)

  # Check knots
  knoterr <- "Minimum number of `knots` for fun=`rcs` is 3"
  if (length(knots)==1 & knots>=1) {
    if (knots<3 & type=="rcs") {
      stop(knoterr)
    }
    nknot <- knots
  } else if (length(knots)>1) {
    if (length(knots)<3 & type=="rcs") {
      stop(knoterr)
    }
    nknot <- length(knots)
  }
  if (length(knots)>1) {
    if (!(all(knots<=1 & all(knots>=0)))) {
      stop("`knots` specified as quantiles must be between 0 and 1")
    }
  }

  # Define function
  base <- "beta.1 * spline[i,m,1]"
  jags <- base
  if (nknot>1) {
    for (i in 2:nknot) {
      temp <- gsub("1", i, base)
      jags <- paste(jags, "+", temp)
    }
  }

  fun <- as.formula(paste("~", jags))

  for (i in 1:nknot) {
    if (get(paste0("pool.", i))=="rel") {
      jags <- gsub(paste0("beta\\.", i), paste0("beta.", i, "[i,k]"), jags)
    } else if (get(paste0("pool.", i))=="abs" & get(paste0("method.", i))=="random") {
      jags <- gsub(paste0("beta\\.", i), paste0("s.beta.", i, "[i,k]"), jags)
    }
  }

  return(list(name=type, knots=knots, nparam=nknot, fun=fun, jags=jags,
              pool.1=pool.1, method.1=method.1, pool.2=pool.2, method.2=method.2,
              pool.3=pool.3, method.3=method.3, pool.4=pool.4, method.4=method.4))

}







funuser <- function(fun, pool.1="rel", method.1="common",
                      pool.2="rel", method.2="common", pool.3="rel", method.3="common",
                      pool.4="rel", method.4="common") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertFormula(fun, add=argcheck)
  checkmate::assertNumeric(knots, null.ok=FALSE, add=argcheck)
  for (i in 1:4) {
    checkmate::assertChoice(get(paste0("pool.", i)), choices=c("rel", "abs"), add=argcheck)
    checkmate::assertChoice(get(paste0("method.", i)), choices=c("common", "random"), add=argcheck)
  }
  checkmate::reportAssertions(argcheck)

  # Check user function
  user.str <- as.character(user.fun[2])
  if (grepl("beta\\.2", user.str)==TRUE & grepl("beta\\.1", user.str)==FALSE) {
    stop("'fun' cannot contain beta.2 if beta.1 is not present")
  } else if (grepl("beta\\.3", user.str)==TRUE & grepl("beta\\.2", user.str)==FALSE | grepl("beta\\.1", user.str)==FALSE) {
    stop("'fun' cannot contain beta.3 if beta.2 and beta.1 are not present")
  }
  if (!(grepl("time", user.str))) {
    stop("'fun' must be a function of beta parameters and time")
  }

  nparam <- 1
  if (grepl("beta\\.4", user.str)) {
    nparam <- 4
  } else if (grepl("beta\\.3", user.str)) {
    nparam <- 3
  } else if (grepl("beta\\.2", user.str)) {
    nparam <- 2
  }

  # Define function
  base <- "beta.1 * time[i,m]"
  jags <- base
  if (nparam>1) {
    for (i in 2:nparam) {
      temp <- gsub("1", i, base)
      jags <- paste(jags, "+", temp)
    }
  }

  for (i in 1:nparam) {
    if (get(paste0("pool.", i))=="rel") {
      jags <- gsub(paste0("beta\\.", i), paste0("beta.", i, "[i,k]"), jags)
    } else if (get(paste0("pool.", i))=="abs" & get(paste0("method.", i))=="random") {
      jags <- gsub(paste0("beta\\.", i), paste0("s.beta.", i, "[i,k]"), jags)
    }
  }

  return(list(fun=fun, nparam=nparam, jags=jags,
              pool.1=pool.1, method.1=method.1, pool.2=pool.2, method.2=method.2,
              pool.3=pool.3, method.3=method.3, pool.4=pool.4, method.4=method.4))

}
