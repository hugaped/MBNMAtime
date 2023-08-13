testthat::context("Testing full set of functions")

# Includes tests for mbnma.run()

# Across a range of datasets and time-course functions:
# Tests running
# Tests default plots of running (including fitplot and devplot)
# Tests default ranking (including cumrank)
# Tests default prediction
# Occasionally tests get.relative


# Tested datasets must have at least 5 treatments
# options are osteopain, goutSUA_CFBcomb, obesityBW_CFB, alog_pcfb, diabetes, hyalarthritis


alldfs <- list(osteopain, goutSUA_CFBcomb, obesityBW_CFB, alog_pcfb, diabetes, hyalarthritis)
datanams <- c("osteopain", "goutSUA_CFBcomb", "obesityBW_CFB", "alog_pcfb", "diabetes", "hyalarthritis")



for (dat in seq_along(alldfs)) {

  datanam <- datanams[dat]
  dataset <- alldfs[[dat]]

  test_that(paste("Testing full set of functions for:", datanam), {

    # Add sdscale
    dataset$standsd <- 2
    dataset$standsd[dataset$studyID %in% unique(dataset$studyID)[c(1,4,6,8)]] <- 0.5

    ### Datasets ####
    network <- mb.network(dataset)


    # Make class data
    df <- dataset

    if ("class" %in% names(dataset)) {
      netclass <- mb.network(df)
    }


    test_that(paste("mb.run functions correctly for:", datanam), {
      skip_on_appveyor()
      skip_on_ci()
      skip_on_cran()

      n.iter=500
      pd <- "pv"

      #set.seed(042189)
      samp <- sample(c(1,2), size=1)
      sdscale <- c(TRUE,FALSE)[samp]

      # NMA
      nma.df <- get.latest.time(network)
      nma <- suppressWarnings(nma.run(nma.df$data.ab, treatments=nma.df$treatments,
                                      method="random", n.iter=500, sdscale=sdscale))


      # Single parameter DR functions
      result <- mb.run(network, fun=tpoly(degree=1, method.1 = "common"),
                          pd="pd.kl", n.iter=n.iter, sdscale=sdscale)
      expect_equal(class(result), c("mbnma", "rjags"))
      expect_equal("d.1" %in% result$parameters.to.save, TRUE)
      expect_equal(result$model.arg$pd, "pd.kl")
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)
      expect_identical(sort(network$studyID), sort(result$model.arg$jagsdata$studyID))





      result <- mb.run(network, fun=tloglin(pool.rate="rel", method.rate="random"),
                          pd="pd.kl", n.iter=n.iter, sdscale=sdscale)
      expect_equal(class(result), c("mbnma", "rjags"))
      expect_equal("sd.rate" %in% result$parameters.to.save, TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result, param=c("rate", "auc")[samp]), NA)
      expect_error(suppressWarnings(devplot(result, dev.type=c("dev", "resdev")[samp])), NA)
      expect_error(fitplot(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)


      if ("class" %in% names(dataset)) {
        result <- mb.run(netclass, fun=tloglin(pool.rate="rel", method.rate="common"),
                            pd="popt", class.effect = list(rate="random"), n.iter=n.iter,
                            sdscale=sdscale)
        expect_equal(class(result), c("mbnma", "rjags"))
        expect_equal("RATE" %in% result$parameters.to.save, TRUE)
        expect_equal("sd.RATE" %in% result$parameters.to.save, TRUE)
        expect_error(plot(result), NA)
        expect_error(rank(result, param=c("RATE", "rate")[samp]), NA)
        expect_error(suppressWarnings(devplot(result)), NA)
        expect_error(fitplot(result), NA)
        expect_error(predict(result), NA)
        expect_error(suppressWarnings(summary(result)), NA)
        expect_error(get.relative(result), NA)
      }



      # Two parameter DR functions
      result <- mb.run(network, fun=temax(method.emax="common", method.et50 = "common"),
                          n.iter=n.iter, pd=pd, sdscale=sdscale)
      expect_equal(all(c("emax", "et50") %in% result$parameters.to.save), TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result, param=c("auc")), NA)
      expect_error(suppressWarnings(devplot(result)), NA)
      expect_error(fitplot(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)


      result <- mb.run(network, fun=temax(pool.et50="abs", method.et50="common", method.emax="random"),
                          n.iter=n.iter, pd=pd, sdscale=sdscale)
      expect_equal("sd.emax" %in% result$parameters.to.save, TRUE)
      expect_equal("et50" %in% rownames(result$BUGSoutput$summary), TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)


      if ("class" %in% names(dataset)) {
        expect_error(mb.run(netclass, fun=temax(p.expon=TRUE, method.et50="random"), corparam=TRUE,
                               class.effect=list(fakeparam="common"), sdscale=sdscale), "The following list element names")

        expect_error(mb.run(netclass, fun=temax(p.expon=TRUE, method.et50="random"), corparam=TRUE,
                                 class.effect=list(et50="common"), n.iter=n.iter, pd=pd, sdscale=sdscale), NA)

        result <- suppressWarnings(mb.run(netclass, fun=temax(method.emax="random"),
                                             class.effect=list(et50="common"), n.iter=n.iter, pd=pd, sdscale=sdscale))
        expect_equal(all(c("emax", "ET50", "et50", "sd.emax") %in% result$parameters.to.save), TRUE)
        expect_error(plot(result), NA)
        expect_error(rank(result), NA)
        expect_error(predict(result), NA)
        expect_error(suppressWarnings(summary(result)), NA)

      }


      result <- mb.run(network, fun=temax(pool.et50="abs", method.et50="random"),
                          n.iter=n.iter, pd=pd, sdscale=sdscale)
      expect_equal(all(c("emax", "et50", "sd.et50") %in% result$parameters.to.save), TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)
      expect_error(get.relative(result), NA)


      # Three parameter DR function
      result <- tryCatch(mb.run(network, fun=temax(pool.emax="rel", pool.et50="abs", pool.hill="abs",
                                                   method.emax="common", method.et50="random", method.hill="common"),
                                   n.iter=n.iter, pd=pd, priors = list(hill="dunif(0.1,5)"),
                                   sdscale=sdscale), error=function(e){})

      if (!is.null(result)) {
        expect_equal(all(c("emax", "et50", "sd.et50", "hill") %in% result$parameters.to.save), TRUE)
        expect_error(plot(result), NA)
        expect_error(rank(result), NA)
        expect_error(predict(result), NA)
        expect_error(suppressWarnings(summary(result)), NA)
        expect_error(get.relative(result), NA)
      }

      result <- mb.run(network, fun=tspline(type="ns", knots=2,
                                            pool.1="rel", pool.2="abs", pool.3="abs",
                                            method.1="common", method.2="random", method.3="common"
                                            ),
                          n.iter=n.iter, pd=pd, sdscale=sdscale)
      expect_equal(all(c("d.1", "beta.2", "sd.beta.2", "beta.3") %in% result$parameters.to.save), TRUE)
      expect_equal(any(grepl("spline", result$model.arg$jagscode)), TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result, param=c("auc", "d.1")[samp]), NA)
      expect_error(suppressWarnings(devplot(result)), NA)
      expect_error(fitplot(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)


      result <- tryCatch(mb.run(network, fun=temax(method.et50="random", pool.hill="abs", method.hill=1.2),
                                   parameters.to.save = "totresdev", sdscale=sdscale,
                                   n.iter=n.iter, pd=pd), error=function(e){})

      if (!is.null(result)) {
        expect_equal("totresdev" %in% result$parameters.to.save, TRUE)
        expect_equal("d.1" %in% result$parameters.to.save, FALSE)
        expect_error(plot(result), "No time-course consistency")
        expect_error(rank(result), "Parameters required for estimation")
        # expect_error(devplot(result), NA)
        # expect_error(fitplot(result), NA)
        expect_error(predict(result), "Parameters required for estimation")
        expect_error(summary(result), "Cannot use")
      }


      # Splines and polynomials
      result <- mb.run(network, fun=tspline(type="bs", knots=2,
                                            pool.1="abs", pool.2="rel", pool.3="abs",
                                            method.1="common", method.2 = "common", method.3="random"),
                          n.iter=n.iter, pd=pd, sdscale=sdscale)
      expect_equal(all(c("beta.1", "d.2", "sd.beta.3", "beta.3") %in% result$parameters.to.save), TRUE)
      expect_equal(all(c("sd.d.2") %in% result$parameters.to.save), FALSE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)
      expect_error(get.relative(result), NA)



      result <- mb.run(network, fun=tspline(type="ns", knots=c(0.2,0.5),
                                            pool.1="abs", pool.2="rel", pool.3="abs",
                                            method.1="common", method.2 = "common", method.3="random"),
                          n.iter=n.iter, pd=pd, sdscale=sdscale)
      expect_equal(all(c("beta.1", "d.2", "sd.beta.3", "beta.3") %in% result$parameters.to.save), TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result, param="d.2"), NA)
      expect_error(suppressWarnings(devplot(result)), NA)
      expect_error(fitplot(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)


      result <- mb.run(network, fun=tpoly(degree=3,
                                          pool.1="abs", pool.2="rel", pool.3="abs",
                                          method.1="common", method.2 = "random", method.3="random"),
                          n.iter=n.iter, pd=pd, sdscale=sdscale)
      expect_equal(all(c("beta.1", "d.2", "sd.beta.3", "beta.3", "sd.beta.2") %in% result$parameters.to.save), TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result, param=c("d.2", "auc")[samp]), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)


      # Test different covariances
      result <- mb.run(network, fun=tpoly(degree=3,
                                          pool.1="abs", pool.2="rel", pool.3="abs",
                                          method.1="common", method.2 = "random", method.3="random"),
                       n.iter=n.iter, pd=pd, sdscale=sdscale,
                       covar="varadj", rho="dunif(-1,1)")
      expect_false(result$BUGSoutput$median$rho==0)
      expect_error(plot(result), NA)
      expect_error(predict(result), NA)
      expect_output(summary(result), "Covariance structure: varadj")

      result <- mb.run(network, fun=tpoly(degree=3,
                                          pool.1="abs", pool.2="rel", pool.3="abs",
                                          method.1="common", method.2 = "random", method.3="random"),
                       n.iter=n.iter/2, pd=pd, sdscale=sdscale,
                       covar="AR1", rho="dunif(0,1)")
      expect_false(result$BUGSoutput$median$rho==0)
      expect_error(plot(result), NA)
      expect_error(predict(result), NA)
      expect_output(summary(result), "Covariance structure: AR1")

      result <- mb.run(network, fun=tpoly(degree=3,
                                          pool.1="abs", pool.2="rel", pool.3="abs",
                                          method.1="common", method.2 = "random", method.3="random"),
                       n.iter=n.iter/2, pd=pd, sdscale=sdscale,
                       covar="CS", rho="dunif(-1,1)")
      expect_false(result$BUGSoutput$median$rho==0)
      expect_error(plot(result), NA)
      expect_error(predict(result), NA)
      expect_output(summary(result), "Covariance structure: CS")
      expect_error(get.relative(result), NA)


      # Test corparam
      expect_equal("inv.R" %in% names(result$model.arg$priors), FALSE)
      expect_equal(result$model.arg$omega, NULL)

      result <- mb.run(network, fun=tpoly(degree=3,
                                          pool.1="rel", pool.2="rel", pool.3="abs",
                                          method.1="common", method.2 = "random", method.3="random"),
                       n.iter=n.iter/2, pd=pd, sdscale=sdscale,
                       covar="CS", rho="dunif(-1,1)", corparam = TRUE)
      expect_equal("rhoparam" %in% names(result$model.arg$priors), TRUE)
      expect_error(get.relative(result), NA)


      # Test UME
      result <- mb.run(network, fun=tpoly(degree=3,
                                          pool.1="abs", pool.2="rel", pool.3="rel",
                                          method.1="common", method.2 = "random", method.3="random"),
                          n.iter=n.iter, pd=pd, UME=TRUE, sdscale=sdscale)
      expect_equal(all(c("beta.1", "d.2", "sd.beta.2", "d.3", "sd.beta.3") %in% result$parameters.to.save), TRUE)
      expect_equal(any(grepl("d\\.2\\[1,2\\]", rownames(result$BUGSoutput$summary))), TRUE)
      expect_equal(any(grepl("d\\.3\\[1,2\\]", rownames(result$BUGSoutput$summary))), TRUE)
      expect_equal(any(grepl("d\\.1\\[1,2\\]", rownames(result$BUGSoutput$summary))), FALSE)

      expect_error(get.relative(result), "cannot be used with UME")
      expect_error(plot(result), "cannot be used with UME")
      expect_error(rank(result), "cannot be used with UME")
      expect_error(suppressWarnings(devplot(result)), NA)
      expect_error(fitplot(result), NA)
      expect_error(predict(result), "UME model can only be used for prediction of direct estimates for a single")
      expect_error(suppressWarnings(summary(result)), NA)


      result <- mb.run(network, fun=tpoly(degree=3,
                                          pool.1="abs", pool.2="rel", pool.3="rel",
                                          method.1="common", method.2 = "random", method.3="random"),
                       n.iter=n.iter, pd=pd, UME="beta.3", sdscale=sdscale)
      expect_equal(any(grepl("d\\.2\\[1,2\\]", rownames(result$BUGSoutput$summary))), FALSE)
      expect_equal(any(grepl("d\\.3\\[1,2\\]", rownames(result$BUGSoutput$summary))), TRUE)
      expect_equal(any(grepl("d\\.1\\[1,2\\]", rownames(result$BUGSoutput$summary))), FALSE)


      # Link functions (include sdscale tests)
      if (datanam %in% c("osteopain", "diabetes", "hyalarthritis")) {

        if (datanam %in% c("diabetes", "hyalarthritis")) {
          expect_error(mbnma.run(network, link="log", n.iter=n.iter, pd=pd),
                       "cannot be used with means (y) that take negative values")
        }

        absdat <- dataset
        absdat$y <- abs(absdat$y)

        absnet <- mb.network(absdat)

        result <- mb.run(absnet, fun=temax(), link="log", n.iter=n.iter, pd=pd,
                            sdscale=sdscale)
        expect_equal(result$model.arg$link, "log")

        expect_error(plot(result), NA)
        expect_error(rank(result), NA)
        expect_error(suppressWarnings(devplot(result)), NA)
        expect_error(fitplot(result), NA)
        expect_error(predict(result), NA)
        expect_error(suppressWarnings(summary(result)), NA)
        expect_error(get.relative(result), NA)
      }


      # Changing priors
      result <- mb.run(network, fun=temax(method.emax = "random"),
                          n.iter=n.iter, pd=pd, sdscale=sdscale)
      prior <- list(sd.emax="dunif(0,5)", et50="dlnorm(1,0.001)")

      runprior <- mb.run(network, fun=temax(method.emax = "random"),
                       n.iter=n.iter, pd=pd, sdscale=sdscale, priors=prior)

      expect_equal(runprior$model.arg$priors$sd.emax, prior$sd.emax)
      expect_equal(runprior$model.arg$priors$et50, prior$et50)
      expect_equal(result$model.arg$priors$et50!=runprior$model.arg$priors$et50, TRUE)
    })

  })
}
