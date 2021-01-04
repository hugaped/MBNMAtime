# Functions for time-course
# Author: Hugo Pedder
# Date created: 2021-01-04

linear <- function(pool.slope="rel", method.slope="common") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(pool.slope, choices=c("rel", "abs"), add=argcheck)
  checkmate::assertChoice(method.slope, choices=c("common", "random"), add=argcheck)
  checkmate::reportAssertions(argcheck)


  fun <- ~ slope * time
  jags <- "beta.1 * time[i,m]"

  if (pool.slope=="rel") {
    jags <- gsub("beta\\.1", "beta.1[i,k]", jags)
  } else if (pool.slope=="abs" & method.slope=="random") {
    jags <- gsub("beta\\.1", "i.beta.1[i,k]", jags)
  }

  # Generate output values
  paramnames <- "slope"
  nparam <- 1

  apool <- poolslope
  names(apool) <- paramnames
  amethod <- method.slope
  names(amethod) <- paramnames

  bpool <- paste0("pool.", 1:nparam)
  names(bpool) <- paramnames
  bmethod <- paste0("method.", 1:nparam)
  names(bmethod) <- paramnames

  out <- list(name="linear", fun=fun, params=paramnames, nparam=nparam, jags=jags,
              apool=apool, amethod=amethod,
              bpool=bpool, bmethod=bmethod)

  class(out) <- "timefun"

  return(out)
}





emax <- function(pool.emax="rel", method.emax="common", pool.et50="rel", method.et50="common") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(pool.emax, choices=c("rel", "abs"), add=argcheck)
  checkmate::assertChoice(method.emax, choices=c("common", "random"), add=argcheck)
  checkmate::assertChoice(pool.et50, choices=c("rel", "abs"), add=argcheck)
  checkmate::assertChoice(method.et50, choices=c("common", "random"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  fun <- ~ (emax * time) / (exp(et50) + time)
  jags <- "(beta.1 * time[i,m]) / (exp(beta.2) + time[i,m])"

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

  # Generate output values
  paramnames <- c("emax", "et50")
  nparam <- 2

  apool <- c(pool.emax, pool.et50)
  names(apool) <- paramnames
  amethod <- c(method.emax, method.et50)
  names(amethod) <- paramnames

  bpool <- paste0("pool.", 1:nparam)
  names(bpool) <- paramnames
  bmethod <- paste0("method.", 1:nparam)
  names(bmethod) <- paramnames

  out <- list(name="emax", fun=fun, params=paramnames, nparam=nparam, jags=jags,
              apool=apool, amethod=amethod,
              bpool=bpool, bmethod=bmethod)
  class(out) <- "timefun"

  message("'et50' parameters are on exponential scale to ensure they take positive values on the natural scale")
  return(out)
}





fractpoly <- function(degree=1, pool.1="rel", method.1="common", pool.2="rel", method.2="common",
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

  if (degree==1) {
    fun <- ~ beta.1 * ifelse(time>0, ifelse(beta.2==0, log(time), time^beta.2), 0)
    jags <- "beta.1 * ifelse(time[i,m]>0, ifelse(beta.2==0, log(time[i,m]), time[i,m]^beta.2), 0)"
  } else if (degree==2) {
    fun <- ~ beta.1 * ifelse(time>0, ifelse(beta.3==0, log(time), time^beta.3), 0) + (beta.2 * ifelse(beta.4==beta.3, ifelse(time>0, ifelse(beta.4==0, log(time)^2, (time^beta.4) * log(time)), 0), ifelse(time>0, ifelse(beta.4==0, log(time), time^beta.4), 0)))
    jags <- "beta.1 * ifelse(time[i,m]>0, ifelse(beta.3==0, log(time[i,m]), time[i,m]^beta.3), 0) + (beta.2 * ifelse(beta.4==beta.3, ifelse(time[i,m]>0, ifelse(beta.4==0, log(time[i,m])^2, (time[i,m]^beta.4) * log(time[i,m])), 0), ifelse(time[i,m]>0, ifelse(beta.4==0, log(time[i,m]), time[i,m]^beta.4), 0)))"
  }

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
  } else {
    stop("'degree' can take either 1 or 2")
  }

  params <- ifelse(degree==1, c("beta.1", "power1"), c("beta.1", "beta.2", "power1", "power2"))

  out <- list(name="fractpoly", fun=fun, degree=degree,
              params=params, nparam=length(params), jags=jags,
              pool.1=pool.1, method.1=method.1,
              pool.2=pool.2, method.2=method.2,
              method.power1=method.power1, method.power2=method.power2,
              bpool=paste0("pool.", 1:(nparam-degree)), bmethod=paste0("method.", 1:nparam))
  class(out) <- "timefun"

  return(out)
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
      jags <- gsub(paste0("beta\\.", i), paste0("i.beta.", i, "[i,k]"), jags)
    }
  }

  out <- list(name=type, fun=fun, knots=knots, params=paste0("beta.", 1:nparam), nparam=nknot, jags=jags,
              pool.1=pool.1, method.1=method.1, pool.2=pool.2, method.2=method.2,
              pool.3=pool.3, method.3=method.3, pool.4=pool.4, method.4=method.4,
              bpool=paste0("pool.", 1:nparam), bmethod=paste0("method.", 1:nparam))
  class(out) <- "timefun"

  return(out)

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
      jags <- gsub(paste0("beta\\.", i), paste0("i.beta.", i, "[i,k]"), jags)
    }
  }

  out <- list(name="funuser", fun=fun, params=paste0("beta.", 1:nparam), nparam=nparam, jags=jags,
              pool.1=pool.1, method.1=method.1, pool.2=pool.2, method.2=method.2,
              pool.3=pool.3, method.3=method.3, pool.4=pool.4, method.4=method.4,
              bpool=paste0("pool.", 1:nparam), bmethod=paste0("method.", 1:nparam))
  class(out) <- "timefun"

  return(out)

}
