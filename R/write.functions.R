# Functions for writing MBNMA models
# Author: Hugo Pedder
# Date created: 2018-09-10



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












