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



print.treat.str <- function(mbnma) {
  betanames <- get.beta.names(mbnma)

  treat.sect <- c()

  # Time course parameters (treatment-level) (generate treat.str)
  for (i in seq_along(betanames)) {
    sect.head <- paste("####", betanames[[i]], "time-course parameter pooling ####", sep=" ")

    data.head <- paste("Parameter", "Median (95%CrI)", sep="\t\t")
    data.head <- paste(data.head, "---------------------------------", sep="\n")

    if (mbnma$model.arg[[names(betanames)[i]]]$pool=="rel") {
      param <- "d"
    } else if (mbnma$model.arg[[names(betanames)[i]]]$pool %in% c("arm", "const")) {
      param <- "beta"
    }

    data.tab <- get.timeparam.str(mbnma, beta=names(betanames)[i], param = param)
    if (!is.null(data.tab)) {
      data.str <- paste(data.head,
                        data.tab,
                        sep="")
    } else if (names(betanames)[i] %in% names(mbnma$model.arg$class.effect)) {
      data.str <- "--- CLASS EFFECTS MODELLED (results shown in separate section below) ---"
    } else if (is.null(data.tab)) {
      data.str <- "--- DATA TOO LONG FOR SUMMARY ---"
    }

    sd.str <- NULL

    # UME
    if ((mbnma$model.arg$UME==TRUE | betanames[[i]] %in% mbnma$model.arg$UME) &
        mbnma$model.arg[[names(betanames)[i]]]$pool=="rel") {
      sect.head <- paste(sect.head,
                         "Unrelated Mean Effects modelled for this parameter (not shown here)", sep="\n")
    }

    # Parameters on exponential scale
    if (mbnma$model.arg$fun=="emax" | mbnma$model.arg$fun=="emax.hill") {
      if (names(betanames)[i] %in% c("beta.2", "et50")) {
        sect.head <- paste(sect.head,
                           "Parameter modelled on exponential scale to ensure it takes positive values on the natural scale", sep="\n")
      }
    }

    # String for pooling
    if (mbnma$model.arg[[names(betanames)[i]]]$pool=="rel") {
      pool <- "relative effects"
    } else if (mbnma$model.arg[[names(betanames)[i]]]$pool=="arm") {
      pool <- "arm-based"
    } else if (mbnma$model.arg[[names(betanames)[i]]]$pool=="const") {
      pool <- "constant (single parameter constant across all studies and treatments within the network)"
    }
    pool.str <- paste("Pooling:", pool, sep=" ")

    # String for method
    if (mbnma$model.arg[[names(betanames)[i]]]$method=="common") {
      method <- "common (fixed) effects\n\nEstimated from the data:\n"
    } else if (mbnma$model.arg[[names(betanames)[i]]]$method=="random") {
      method <- "random effects\n\nEstimated from the data:\n"

      if (grepl("beta", betanames[[i]])) {
        suffix <- strsplit(betanames[[i]], split="\\.")[[1]][2]
      } else {
        suffix <- betanames[[i]]
      }

      sd.str <- "\n# Between-study SD\n"
      match.1 <- paste0("^sd\\.(beta\\.)?", suffix)
      #match.2 <- paste0("^sd\\.beta", betanames[[i]])
      temp <- mbnma$BUGSoutput$summary[grepl(match.1, rownames(mbnma$BUGSoutput$summary)),
                                       c(3,5,7)]
      nametemp <- rownames(mbnma$BUGSoutput$summary)[grepl(match.1, rownames(mbnma$BUGSoutput$summary))]
      if (!is.vector(temp)) {stop("temp should only be length 1")}

      sd.str <- paste(sd.str, data.head,
                      paste(nametemp, neatCrI(temp), sep="\t\t"),
                      sep="\n")

    } else if (is.numeric(mbnma$model.arg[[names(betanames)[i]]]$method)) {
      method <- paste("none\n\nAssigned a numeric value:",
                      mbnma$model.arg[[names(betanames)[i]]]$method,
                      sep=" ")
      data.str <- ""
    }
    method.str <- paste("Method:", method, sep=" ")

    treat.str <- paste(sect.head, pool.str, method.str, data.str, sep="\n")

    if (!is.null(sd.str)) {
      treat.str <- paste(treat.str, sd.str, sep="\n")
    }

    treat.sect <- paste(treat.sect, treat.str, "", sep="\n\n")
  }
  return(treat.sect)
}





print.class.str <- function(mbnma) {
  if (length(mbnma$model.arg$class.effect)>0) {
    classes <- mbnma$model.arg$class.effect

    head <- "\n#### Class effects ####\n"
    data.head <- paste("Parameter", "Median (95%CrI)", sep="\t")
    data.head <- paste(data.head, "---------------------------------", sep="\n")

    data.str <- c(data.head)
    class.str <- c()
    sd.str <- NULL
    for (i in seq_along(classes)) {
      betaparam <- names(classes)[i]

      if (!is.null(mbnma$model.arg$arg.params)) {
        wrapparam <- mbnma$model.arg$arg.params$wrap.params[
          mbnma$model.arg$arg.params$run.params==names(classes)[i]]

        class.str <- c(class.str, paste("Class effect on",
                                        paste0(wrapparam,":"), classes[[i]], "\n",
                                        sep=" "))
      } else {
        wrapparam <- ""
        class.str <- c(class.str, paste("Class effect on",
                                        paste0(betaparam,":"), classes[[i]], "\n",
                                        sep=" "))
        }


      # betaparam <- mbnma$model.arg$arg.params$run.params[
      #   mbnma$model.arg$arg.params$wrap.params==names(classes)[i]]

      if (mbnma$model.arg[[betaparam]]$pool=="rel") {
        data.str <- paste0(data.str,
                      get.timeparam.str(mbnma, beta = names(classes)[i], param = "D"))
      } else if (mbnma$model.arg[[betaparam]]$pool=="arm") {
        data.str <- c(data.str,
                      get.timeparam.str(mbnma, beta = names(classes)[i], param = "BETA"))
      }

      if (classes[[i]]=="random") {
        sd.str <- "\n# Within-class SD\n"
        match.1 <- paste0("(", strsplit(names(classes)[i], split="\\.")[[1]][2], ")?")
        match.2 <- paste0("(", wrapparam, ")?")
        match <- paste0("^sd\\.[A-Z]+\\.", match.1, match.2)
        temp <- mbnma$BUGSoutput$summary[grepl(match, rownames(mbnma$BUGSoutput$summary)),
                                         c(3,5,7)]
        if (!is.vector(temp)) {stop("temp should only be length 1")}

        sd.str <- paste(sd.str, data.head,
                        paste(rownames(mbnma$BUGSoutput$summary)[grepl(match, rownames(mbnma$BUGSoutput$summary))],
                          neatCrI(temp), sep="\t\t"),
                        sep="\n")
        # sd.str <- "\n# Within-class SD\n"
        # match <- grepl("^sd\\.[A-Z]+\\.", names(classes)[i])
        # temp <- mbnma$BUGSoutput$summary[grepl(match, rownames(mbnma$BUGSoutput$summary)),
        #                                  c(3,5,7)]
        # sd.str <- paste(sd.str, data.head, neatCrI(temp))
      }
    }
    class.sect <- c(head,
                    paste(class.str, collapse="\n"),
                    paste(data.str, collapse="\n"))

    if (!is.null(sd.str)) {
      class.sect <- c(class.sect, sd.str, "\n")
    }

    class.sect <- paste(class.sect, collapse="\n")

    return(class.sect)
  }
}






print.overall.str <- function(mbnma) {

  # Print title
  title <- "========================================\nTime-course MBNMA\n========================================\n"

  # Print time-course function
  overall.sect <- paste("Time-course function:", mbnma$model.arg$fun, sep=" ")

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

  return(overall.sect)
}




print.cor.str <- function(mbnma) {
  if (!is.null(mbnma$model.arg$rho)) {
    head <- "\n#### Correlation between time points ####\n"
    cov.str <- paste("Covariance structure:", mbnma$model.arg$covar, sep=" ")

    if (is.numeric(mbnma$model.arg$rho)) {
      rho.str <- paste("Rho assigned a numeric value:", mbnma$model.arg$rho, sep=" ")
    } else if (mbnma$model.arg$rho=="estimate") {
      data.head <- paste("Parameter", "Median (95%CrI)", sep="\t")
      data.head <- paste(data.head, "---------------------------------", sep="\n")
      if ("rho" %in% mbnma$parameters.to.save) {
        data.str <- neatCrI(mbnma$BUGSoutput$summary[
          rownames(mbnma$BUGSoutput$summary)=="rho",
          c(3,5,7)])
        data.str <- paste("rho", data.str, sep="\t\t")
      } else {
        data.str <- "<rho not monitored in parameters.to.save>"
      }
      data.str <- paste(data.head, data.str, sep="\n")

      rho.str <- paste("Rho estimated from the data:\n", data.str, sep="\n")
    }
    cor.sect <- paste(head, cov.str, rho.str, "", sep="\n")
    return(cor.sect)
  }
}




print.modfit.str <- function(x) {
  totresdev.str <- c()

  head <- "#### Model Fit Statistics ####\n"

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




#' Print summary MBNMA results to the console
#'
#' @inheritParams predict.mbnma
#' @param ... further arguments passed to or from other methods
#'
#' @export
summary.mbnma <- function(object, ...) {
  checkmate::assertClass(object, "mbnma")

  # State that function does not work if "parameters.to.save" has been specified
  if (!is.null(object$model.arg$parameters.to.save)) {
    stop("Cannot use `summary()` method if `parameters.to.save` have been assigned. Use `print()` instead.")
  }

  # Check for rhat < 1.02
  rhat.warning(object)

  # Overall section
  overall.sect <- print.overall.str(object)

  # Print treatment-level section
  treat.sect <- print.treat.str(object)

  # Correlation section
  cor.sect <- print.cor.str(object)

  # Class-effect section
  class.sect <- print.class.str(object)

  # Model fit statistics section
  modfit.sect <- print.modfit.str(object)

  output <- paste(overall.sect, treat.sect, cor.sect, class.sect, modfit.sect, sep="\n")
  cat(output, ...)
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
    warning(paste0(msg, paste(rhats, collapse="\n")))
  }
}







#' Takes node-split results and produces summary data frame
#'
#' @param object An object of class `"mb.nodesplit"` generated by `mb.nodeplit()`
#' @param ... further arguments passed to or from other methods
#'
#' @return
#' A data frame of summary node-split results with the following variables:
#' * `Comparison` The treatment comparison on which a node-split has been performed
#' * `Time.Param` The time-course parameter on which a node-split has been performed
#' * `Evidence` The evidence contribution for the given comparison (either "Direct" or "Indirect")
#' * `Median` The posterior median
#' * `2.5%` The lower 95% credible interval limit
#' * `97.5%` The upper 95% credible interval limit
#' * `p.value` The Bayesian p-value for the overlap between direct and indirect evidence for
#' the given comparison (it will therefore have an identical value for direct and indirect evidence
#' within a particular comparison and time-course parameter)
#'
#' @export
summary.mb.nodesplit <- function(object, ...) {
  checkmate::assertClass(object, "mb.nodesplit")

  sum.mat <- matrix(ncol=3)
  comp <- vector()
  time.param <- vector()
  evidence <- vector()
  pvals <- vector()

  for (i in seq_along(object)) {
    for (k in seq_along(object[[i]])) {
      post <- object[[i]][[k]]$quantiles

      sum.mat <- rbind(sum.mat, post$direct)
      evidence <- c(evidence, "Direct")

      sum.mat <- rbind(sum.mat, post$indirect)
      evidence <- c(evidence, "Indirect")

      pvals <- c(pvals, rep(object[[i]][[k]]$p.values, 2))
      time.param <- c(time.param, rep(names(object[[i]])[k], 2))
      comp <- c(comp, rep(names(object)[i], 2))
    }
  }
  sum.mat <- round(sum.mat[-1,], digits = max(3L, getOption("digits") - 5L))
  pvals <- round(pvals, max(3L, getOption("digits") - 5L))

  sum.df <- data.frame(comp, time.param,
                       evidence, sum.mat[,2],
                       sum.mat[,1], sum.mat[,3],
                       pvals, ...
  )

  names(sum.df) <- c("Comparison", "Time.Param", "Evidence", "Median",
                     "2.5%", "97.5%", "p.value")

  return(sum.df)

}



#' Prints basic results from a node-split to the console
#'
#' @param x An object of class `"mb.nodesplit"` generated by `mb.nodeplit()`
#' @param groupby A character object that can take the value `"time.param"` to present
#' results grouped by time-course parameter (the default) or `"comparison"` to present
#' results grouped by treatment comparison.
#' @param ... further arguments passed to or from other methods
#'
#' @export
print.mb.nodesplit <- function(x, groupby="time.param", ...) {

  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mb.nodesplit", add=argcheck)
  checkmate::assertChoice(groupby, choices=c("time.param", "comparison"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  width <- "\t\t"
  output <- "========================================\nNode-splitting analysis of inconsistency\n========================================"

  if (groupby=="time.param") {
    params <- names(x[[1]])
    colnam <- "comparison\tp.value\t\t\tMedian (95% CrI)"
    for (i in seq_along(params)) {
      paramname <- paste("\n\n\n####", params[i], "####\n", sep=" ")
      paramsect <- colnam

      for (k in seq_along(x)) {
        pval <- signif(x[[k]][[params[i]]]$p.values,
                      max(3L, getOption("digits") - 3L))
        tab <- x[[k]][[params[i]]]$quantiles

        heading <- paste(names(x)[k], pval, sep=width)
        direct <- paste("-> direct", "", neatCrI(tab$direct), sep=width)
        indirect <- paste("-> indirect", "", neatCrI(tab$indirect), sep=width)

        out <- paste(heading, direct, indirect, sep="\n")

        paramsect <- append(paramsect, out)
      }
      groupsect <- paste(c(paramname, paramsect), collapse="\n")
      output <- append(output, groupsect)
    }
  } else if (groupby=="comparison") {
    params <- names(x)
    colnam <- "time parameter\tp.value\t\t\tMedian (95% CrI)"
    for (i in seq_along(params)) {
      paramname <- paste("\n\n\n####", params[i], "####\n", sep=" ")
      paramsect <- colnam

      for (k in seq_along(x[[1]])) {
        pval <- signif(x[[params[i]]][[k]]$p.values,
                      max(3L, getOption("digits") - 3L))
        tab <- x[[params[i]]][[k]]$quantiles

        heading <- paste(names(x[[i]])[k], pval, sep=width)
        direct <- paste("-> direct", "", neatCrI(tab$direct), sep=width)
        indirect <- paste("-> indirect", "", neatCrI(tab$indirect), sep=width)

        out <- paste(heading, direct, indirect, sep="\n")

        paramsect <- append(paramsect, out)
      }
      groupsect <- paste(c(paramname, paramsect), collapse="\n")
      output <- append(output, groupsect)
    }
  }
  cat(output, ...)
}






#' Prints summary of mb.predict object
#'
#' Prints a summary table of the mean of MCMC iterations at each time point
#' for each treatment
#'
#' @param object An object of class `"mb.predict"`
#' @param ... further arguments passed to or from other methods
#'
#' @return A matrix containing times at which responses have been predicted (`time`)
#' and an additional column for each treatment for which responses have been predicted.
#' Each row represents mean MCMC predicted responses for each treatment at a particular
#' time.
#'
#' @examples
#' \donttest{
#' # Define network
#' network <- mb.network(obesityBW_CFB, reference="plac")
#'
#' # Run an MBNMA with a quadratic time-course function
#' quad <- mb.quadratic(network,
#'   beta.1=list(pool="rel", method="common"),
#'   beta.2=list(pool="rel", method="common"),
#'   intercept=TRUE)
#'
#' # Predict responses
#' pred <- predict(quad, times=c(0:50), treats=c(1:5),
#'   ref.estimate = network$data.ab[network$data.ab$treatment==1,],
#'   baseline=10)
#'
#' # Generate summary of predictions
#' summary(pred)
#' }
#' @export
summary.mb.predict <- function(object, ...) {
  checkmate::assertClass(object, "mb.predict")

  sumlist <- object[["summary"]]

  time <- sumlist[[1]]$time
  treats <- names(sumlist)
  #treats <- unlist(lapply(treats, FUN=function(x) paste0("treat_", x)))
  sum.df <- time
  for (i in seq_along(sumlist)) {
    sum.df <- cbind(sum.df, sumlist[[i]]$mean)
  }
  #sum.df <- data.frame(sum.df)
  colnames(sum.df) <- c("time", treats)

  #print(sum.df, digits = max(3, getOption("digits")-3), max=ncol(sum.df)*10)

  #print(sum.df)
  #return(invisible(sum.df))
  return(sum.df)
}





#' Prints a summary of rankings for each parameter
#'
#' @inheritParams plot.mb.rank
#' @param ... further arguments passed to or from other methods
#'
#' @export
print.mb.rank <- function(x, ...) {

  output <- "========================================\nTreatment rankings\n========================================"

  for (param in seq_along(names(x))) {
    head <- paste("####", names(x)[param], "####", sep=" ")

    sumtab <- x[[names(x)[param]]]$summary
    data.str <- "Treatment\tMedian rank (95% CrI)"

    for (i in 1:nrow(sumtab)) {
      row <- paste(sumtab[i,1], neatCrI(sumtab[i, c(4,6,8)]), sep="\t\t")
      data.str <- paste(data.str, row, sep="\n")
    }

    sect <- paste(head, data.str, sep="\n")

    output <- paste(output, sect, sep="\n\n")

  }
  return(cat(output, ...))
}





#' Print mb.network information to the console
#'
#' @param x An object of class `mb.network`.
#' @param ... further arguments passed to or from other methods
#'
#' @export
print.mb.network <- function(x,...) {
  nn <- names(x)
  ll <- length(x)
  if (length(nn) != ll)
    nn <- paste("Component", seq.int(ll))
  for (i in seq_len(ll)) {
    cat(nn[i], ":\n")
    if (is.data.frame((x[[i]]))) {
      print(x[[i]], max=ncol(x[[i]])*6, ...)
    } else {
      print(x[[i]], ...)
    }
    cat("\n")
  }
  invisible(x)
}



#' Print summary information from an mb.predict object
#'
#' @param x An object of `class("mb.predict")` generated by `predict.mbnma()`
#' @param ... further arguments passed to or from other methods
#' @inheritParams summary.mb.predict
#'
#' @export
print.mb.predict <- function(x, ...) {

  sum.df <- summary(x)

  sumlist <- x[["summary"]]

  if (!("1" %in% names(sumlist))) {
    cat("Responses have not been predicted for the network reference treatment\n")
  }

  msg <- paste0("Predicted responses at ", nrow(sum.df), " different follow-up times ",
                "for treatments: ", paste(names(sumlist), collapse=", "), "\n\n")
  cat(msg)

  print(sum.df, digits = max(3, getOption("digits")-3), max=ncol(sum.df)*10, ...)
}







#' MBNMA ggplot2 theme style
#' @noRd
theme_mbnma <- function(...) {
  ggplot2::theme_bw(...) +
    ggplot2::theme(
      # change stuff here
      panel.background  = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill="gray96", colour=NA),
      legend.background = ggplot2::element_rect(fill="transparent", colour=NA),
      legend.key = ggplot2::element_rect(fill="transparent", colour=NA),

      # From multinma
      #panel.border = ggplot2::element_rect(colour = "grey70", fill = NA),
      panel.grid.major = ggplot2::element_line(colour = "grey95"),
      panel.grid.minor = ggplot2::element_line(colour = "grey95"),
      strip.background = ggplot2::element_rect(colour = "black",
                                               fill = "lightsteelblue1"),
      strip.text = ggplot2::element_text(colour = "black")
    )
}
