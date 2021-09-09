# Functions for printing and summarising results from MBNMA models
# Author: Hugo Pedder
# Date created: 2019-07-24


get.beta.names <- function(mbnma) {
  betanames <- list()
  for (i in 1:4) {
    if (!is.null(mbnma$model.arg[[paste0("beta.",i)]])) {
      if (is.null(mbnma$model.arg$arg.params)) {
        betanames[[paste0("beta.", i)]] <- paste0("beta.", i)
      } else if (!is.null(mbnma$model.arg$arg.params)) {
        temp <- mbnma$model.arg$arg.params$wrap.params[
          mbnma$model.arg$arg.params$run.params==paste0("beta.", i)
          ]
        betanames[[paste0("beta.", i)]] <- temp
      }
    }
  }

  return(betanames)
}




neatCrI <- function(vals, digits=3) {
  vals <- signif(vals, digits = digits)
  neat <- paste0(vals[2], " (", vals[1], ", ", vals[3], ")")
  return(neat)
}




get.timeparam.str <- function(mbnma, beta=NULL, param="d") {
  betanames <- get.beta.names(mbnma)

  if (grepl("beta", betanames[[beta]])) {
    temp <- strsplit(betanames[[beta]], split="\\.")[[1]][2]
  } else {
    temp <- betanames[[beta]]
  }

  match <- paste0("^", param, "\\.", temp, "(\\[[0-9]+\\])?")

  sum.mat <- mbnma$BUGSoutput$summary[grepl(match, rownames(mbnma$BUGSoutput$summary)),
                                      c(3,5,7)]

  # Check for UME
  # if (any(grepl("\\[[0-9]+,[0-9]+\\]", rownames(sum.mat)))) {
  #   sum.mat <- NULL
  # }

  # Check if size of sum.mat is too great for printing
  if (is.matrix(sum.mat)) {
    if (nrow(sum.mat) > length(mbnma$network$treatments)) {
      sum.mat <- NULL
    }
  }

  if (length(sum.mat)>0) {
    tab.str <- c()

    if (is.matrix(sum.mat)) {
      if (any(grepl("^d\\..+\\[1\\]", rownames(sum.mat)[1]) |
              grepl("^D\\..+\\[1\\]", rownames(sum.mat)[1]))) {
        tab.str <- paste(tab.str,
                         paste(rownames(sum.mat)[1], "Network reference",
                               sep="\t"),
                         sep="\n"
        )
        count <- 2
      } else {
        tab.str <- paste(tab.str,
                         paste(rownames(sum.mat)[1],
                               neatCrI(sum.mat),
                               sep="\t"),
                         sep="\n"
        )
        count <- 1
      }

      for (i in count:nrow(sum.mat)) {
        tab.str <- paste(tab.str,
                         paste(rownames(sum.mat)[i], neatCrI(sum.mat[i,]),
                               sep="\t"),
                         sep="\n"
        )
      }
    } else if (is.vector(sum.mat)) {
      tab.str <- paste(tab.str,
                       paste(rownames(mbnma$BUGSoutput$summary)[grepl(match, rownames(mbnma$BUGSoutput$summary))],
                             neatCrI(sum.mat),
                             sep="\t\t"),
                       sep="\n"
      )
    }

    return(tab.str)
  } else {
    return(NULL)
  }

}



print.treat.str <- function(mbnma, digits=4, ...) {
  fun <- mbnma$model.arg$fun
  treat.df <- mbnma$BUGSoutput$summary

  treat.sect <- c()

  for (i in seq_along(fun$apool)) {
    cat(crayon::underline(crayon::bold(paste0(names(fun$apool)[i], " parameter\n"))))

    paramdet <- vector()
    if (fun$apool[i]=="rel") {
      paramdet <- append(paramdet, paste0("Pooling: relative effects"))
    } else if (fun$apool[i]=="abs") {
      paramdet <- append(paramdet, paste0("Pooling: absolute effects"))
    }

    if (is.numeric(fun$amethod[i])) {
      paramdet <- append(paramdet, paste0("Assigned a numeric value: ", fun$amethod[i]))
    } else {
      if (fun$amethod[i]=="common") {
        paramdet <- append(paramdet, paste0("Method: common treatment effects"))
      } else if (fun$amethod[i]=="random") {
        paramdet <- append(paramdet, paste0("Method: random treatment effects"))
      }
    }

    if (names(fun$apool)[i] %in% c("et50", "hill")) {
      paramdet <- append(paramdet, "Parameter modelled on exponential scale to ensure it takes positive values on the natural scale")
    }

    if (names(fun$amethod)[i] %in% names(mbnma$model.arg$class.effect)) {
      paramdet <- append(paramdet,  "Class effects modelled for this parameter")
    }

    if (any(c(names(fun$apool)[i], TRUE) %in% mbnma$model.arg$UME)) {
      paramdet <- append(paramdet, "Unrelated Mean Effect results modelled for this parameter.\nToo many parameters to display in this summary")

      cat(paste(paramdet, collapse="\n"))
    } else {
      cat(paste(paramdet, collapse="\n"))

      # Select summary rows that correspond to parameter of interest
      index <- grepl(paste0("^", names(fun$apool)[i]), rownames(treat.df)) |
        grepl(paste0("^d\\.", i), rownames(treat.df)) |
        grepl(paste0("^beta\\.", i), rownames(treat.df))

      param.df <- data.frame(treat=mbnma$network$treatments,
                             param=rownames(treat.df)[index],
                             median=treat.df[index,'50%'],
                             l95=treat.df[index,'2.5%'],
                             u95=treat.df[index,'97.5%']
                             )
      rownames(param.df) <- NULL
      print(knitr::kable(param.df, col.names = c("Treatment", "Parameter", "Median", "2.5%", "97.5%"), digits=digits, ...))
    }

    # Add between-study SD for parameter
    if (fun$amethod[i]=="random") {
      cat("\nBetween-study SD modelled for this parameter:")

      index <- grepl(paste0("^sd\\.", names(fun$apool)[i]), rownames(treat.df)) |
        grepl(paste0("^sd\\.d\\.", i), rownames(treat.df)) |
        grepl(paste0("^sd\\.beta\\.", i), rownames(treat.df))

      param.df <- data.frame(param=rownames(treat.df)[index],
                             median=treat.df[index,'50%'],
                             l95=treat.df[index,'2.5%'],
                             u95=treat.df[index,'97.5%']
      )
      rownames(param.df) <- NULL
      print(knitr::kable(param.df, col.names = c("Parameter", "Median", "2.5%", "97.5%"), digits=digits, ...))
    }
    cat("\n\n")
  }
}





print.class.str <- function(mbnma, digits=4, ...) {
  if (length(mbnma$model.arg$class.effect)>0) {
    treat.df <- mbnma$BUGSoutput$summary

    classes <- mbnma$model.arg$class.effect

    cat(crayon::bold(crayon::underline("Class Effects\n")))

    for (i in seq_along(classes)) {

      cat(paste0("\nClass effects for ", names(classes)[i], "\n"))

      if ("common" %in% classes[[i]]) {
        cat("Common (fixed) class effects")
      } else if ("random" %in% classes[[i]]) {
        cat("Random (exchangeable) class effects")
      }

      # Select summary rows that correspond to parameter of interest
      index <- grepl(paste0("^", toupper(names(classes)[i])), rownames(treat.df)) |
        grepl(paste0("^D\\.", which(mbnma$model.arg$fun$params==names(classes)[i])), rownames(treat.df))

      param.df <- data.frame(treat=mbnma$network$classes,
                             param=rownames(treat.df)[index],
                             median=treat.df[index,'50%'],
                             l95=treat.df[index,'2.5%'],
                             u95=treat.df[index,'97.5%']
      )
      rownames(param.df) <- NULL
      print(knitr::kable(param.df, col.names = c("Class", "Parameter", "Median", "2.5%", "97.5%"), digits=digits, ...))

      # Add between-study SD for parameter
      if (classes[[i]]=="random") {
        cat("\nWithin-class SD modelled for this parameter:")

        index <- grepl(paste0("^sd\\.", toupper(names(classes)[i])), rownames(treat.df)) |
          grepl(paste0("^sd\\.D\\.", which(mbnma$model.arg$fun$params==names(classes)[i])), rownames(treat.df))

        param.df <- data.frame(param=rownames(treat.df)[index],
                               median=treat.df[index,'50%'],
                               l95=treat.df[index,'2.5%'],
                               u95=treat.df[index,'97.5%']
        )
        rownames(param.df) <- NULL
        print(knitr::kable(param.df, col.names = c("Parameter", "Median", "2.5%", "97.5%"), digits=digits, ...))
      }
      cat("\n")

    }
  }
}






print.overall.str <- function(mbnma) {

  # Print title
  title <- crayon::bold("========================================\nTime-course MBNMA\n========================================\n")

  # Prep time-course function
  timefun <- mbnma$model.arg$fun$name

  if (timefun %in% c("poly", "fpoly")) {
    timefun <- paste0(timefun, " (degree = ", mbnma$model.arg$fun$nparam, ")")
  }
  if (timefun=="rcs") {
    timefun <- paste0("Restricted cubic spline (knots = ", paste(mbnma$model.arg$fun$knots, collapse=", "), ")")
  } else if (timefun=="ns") {
    timefun <- paste0("Natural cubic spline (knots = ", paste(mbnma$model.arg$fun$knots, collapse=", "), ")")
  } else if (timefun=="ls") {
    timefun <- paste0("Piecewise linear spline (knots = ", paste(mbnma$model.arg$fun$knots, collapse=", "), ")")
  } else if (timefun=="bs") {
    timefun <- paste0("B-spline (knots = ", paste(mbnma$model.arg$fun$knots, collapse=", "), "; degree = ",
                      mbnma$model.arg$fun$degree, ")")
  }

  # Print time-course function
  overall.sect <- paste("Time-course function:", timefun, sep=" ")

  # Add info on intercept
  if (mbnma$model.arg$intercept==TRUE) {
    overall.sect <- paste(overall.sect, "Data modelled with intercept", sep="\n")
  } else if (mbnma$model.arg$intercept==FALSE) {
    overall.sect <- paste(overall.sect, "Data modelled without intercept (change from baseline data assumed)", sep="\n")
  }

  # positive.scale
  if (mbnma$model.arg$positive.scale==TRUE) {
    overall.sect <- paste(overall.sect, "Responses restricted to taking positive values", sep="\n")
  }

  overall.sect <- paste(title, overall.sect, sep="\n")

  cat(overall.sect)
}




print.cor.str <- function(mbnma, digits=4) {
  if (!is.null(mbnma$model.arg$rho)) {
    cat(paste0("\n", crayon::underline(crayon::bold("Correlation between time points")), "\n"))
    cat(paste("Covariance structure:", mbnma$model.arg$covar, "\n", sep=" "))

    if (is.numeric(mbnma$model.arg$rho)) {
      cat(paste("Rho assigned a numeric value:", mbnma$model.arg$rho, sep=" "))
    } else if (is.character(mbnma$model.arg$rho)) {
      if ("rho" %in% mbnma$parameters.to.save) {
        cat("Rho estimated from the data:")
        rho <- mbnma$BUGSoutput$summary["rho",]

        rho <- data.frame(param="rho",
                               median=rho['50%'],
                               l95=rho['2.5%'],
                               u95=rho['97.5%']
        )
        rownames(rho) <- NULL
        print(knitr::kable(rho, col.names = c("Parameter", "Median", "2.5%", "97.5%"), digits=digits))

      } else {
        cat("<rho not monitored in parameters.to.save>")
      }
    }
    cat("\n\n")
  }
}




print.modfit.str <- function(x) {
  totresdev.str <- c()

  head <- crayon::bold("#### Model Fit Statistics ####\n")

  # pD
  pd.str <- "Effective number of parameters:"
  if (x$model.arg$pd=="pv") {
    pd <- "pD (pV) calculated using the rule, pD = var(deviance)/2 ="
  } else if (x$model.arg$pd=="plugin") {
    pd <- "pD calculated using the plug-in method ="
  } else if (x$model.arg$pd=="pd.kl") {
    pd <- "pD calculated using the Kullback-Leibler divergence ="
  } else if (x$model.arg$pd=="popt") {
    pd <- "pD calculated using an optimism adjustment ="
  }
  pd.str <- paste(pd.str, paste(pd, round(x$BUGSoutput$pD,0), sep=" "), sep="\n")

  # Deviance
  dev <- x$BUGSoutput$summary[
    rownames(x$BUGSoutput$summary)=="deviance", 5]
  dev.str <- paste("Deviance =", round(dev, 0), sep=" ")

  # Totresdev
  if ("totresdev" %in% x$parameters.to.save) {
    totresdev <- round(
      x$BUGSoutput$summary[
        rownames(x$BUGSoutput$summary)=="totresdev", 5],
      0)
  } else {
    totresdev <- "NOT MONITORED IN MODEL"
  }
  totresdev.str <- paste("Residual deviance =", totresdev, sep=" ")

  dic <- x$BUGSoutput$DIC
  dic.str <- paste("Deviance Information Criterion (DIC) =", round(dic, 0), "\n", sep=" ")

  modfit.sect <- paste(head, pd.str, dev.str, totresdev.str, dic.str, sep="\n")
  return(modfit.sect)
}








# Gets summary columns for parameters
getsum <- function(params, mbnma) {
  if (length(params)==1) {
    rows <- mbnma$BUGSoutput$summary[grepl(params, rownames(mbnma$BUGSoutput$summary)),]

  } else if (length(params)==2) {
    rows <- mbnma$BUGSoutput$summary[
      grepl(params[1], rownames(mbnma$BUGSoutput$summary)) |
        grepl(params[2], rownames(mbnma$BUGSoutput$summary))
      ,]
  }

  if (is.matrix(rows)) {
    rows <- rows[,colnames(rows) %in% c("mean", "2.5%", "97.5%")]
  } else if (is.vector(rows)) {
    rows <- rows[names(rows) %in% c("mean", "2.5%", "97.5%")]
  }

  rows <- round(rows, digits = max(3, getOption("digits")-3))

  return(rows)

}




rhat.warning <- function(mbnma, cutoff=1.02) {
  rhats <- mbnma$BUGSoutput$summary[,colnames(mbnma$BUGSoutput$summary)=="Rhat"]
  rhats <- names(rhats)[rhats>cutoff]
  if (length(rhats)>0) {
    msg <- paste0("The following parameters have Rhat values > ",
                  cutoff,
                  " which could be due to convergence issues:\n")
    cat(crayon::red((paste0(msg, paste(rhats, collapse="\n")))))
  }
}




