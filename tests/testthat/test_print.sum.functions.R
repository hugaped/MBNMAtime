testthat::context("Testing print.sum.functions")

datalist <- list(osteopain=osteopain, copd=copd, goutSUA_CFBcomb=goutSUA_CFBcomb,
                 hyalarthritis=hyalarthritis, diabetes=diabetes, alog_pcfb=alog_pcfb)

testthat::test_that("print.sum tests pass correctly", {

  testthat::expect_equal(1,1) # Avoids empty tests

  skip_on_appveyor()
  skip_on_ci()
  skip_on_cran()

  n.iter <- 200
  n.burnin <- 100
  n.thin <- 1

  seed <- 890421

  for (i in seq_along(datalist)) {

    print(names(datalist)[i])

    network <- mb.network(datalist[[i]])

    # REMOVE SUPPRESSWARNINGS FROM VERSION 0.2.3 ONWARNS
    suppressWarnings(
    emax1 <- mb.run(network,
                    fun=temax(pool.emax="abs", method.emax="random",
                              pool.et50="rel", method.et50="common"),
                    n.chain=3, n.iter=200, n.burnin=100, jags.seed=seed,
                    rho="dunif(0,1)", covar="AR1",
                    positive.scale = TRUE, intercept = TRUE,
                    UME=TRUE)
    )

    network.gout <- mb.network(goutSUA_CFBcomb)

    # REMOVE SUPPRESSWARNINGS FROM VERSION 0.2.3 ONWARNS
    suppressWarnings(
    piece1 <- mb.run(network.gout, fun=tspline(type="ls",
                                               knots=0.6,
                                               pool.1="abs", method.1="random",
                                               pool.2="rel", method.2="common"),
                     intercept=FALSE,
                     n.chain=3, n.iter=200, n.burnin=100, jags.seed=seed,
                     rho=0.2, covar="CS", pd="popt",
                     class.effect = list(beta.2="random"))
    )


    ##############################################
    ###########    TESTS    ######################
    ##############################################

    testthat::test_that("print.overall.str functions correctly", {

      testthat::expect_output(print.overall.str(emax1), "Data modelled with intercept")
      testthat::expect_output(print.overall.str(emax1), "Time-course function: emax")
      testthat::expect_output(print.overall.str(emax1), "Responses restricted to taking positive values")

      testthat::expect_output(print.overall.str(piece1), "Data modelled without intercept")
      testthat::expect_output(print.overall.str(piece1), "Time-course function: Piecewise linear spline \\(knots = 0\\.6\\)")

    })

    testthat::test_that("print.treat.str functions correctly", {

      # Emax
      testthat::expect_output(print.treat.str(emax1), "Unrelated Mean Effect results modelled for this parameter")
      testthat::expect_output(print.treat.str(emax1), "Too many parameters")
      testthat::expect_output(print.treat.str(emax1), "Pooling: absolute effects")
      testthat::expect_output(print.treat.str(emax1), "Method: random treatment effects")
      testthat::expect_output(print.treat.str(emax1), "Between-study SD")
      testthat::expect_output(print.treat.str(emax1), "Pooling: relative effects")
      testthat::expect_output(print.treat.str(emax1), "Method: common treatment effects")


      testthat::expect_output(print.treat.str(piece1), "Class effects modelled for this")
      testthat::expect_output(print.treat.str(piece1), "Pooling: absolute effects")
      testthat::expect_output(print.treat.str(piece1), "Method: random treatment effects")
      testthat::expect_output(print.treat.str(piece1), "Between-study SD")
      testthat::expect_output(print.treat.str(piece1), "Pooling: relative effects")
      testthat::expect_output(print.treat.str(piece1), "Method: common treatment effects")
      testthat::expect_output(print.treat.str(piece1), "beta\\.1")
      testthat::expect_output(print.treat.str(piece1), "sd\\.beta\\.1")

    })





    testthat::test_that("print.cor.str functions correctly", {

      testthat::expect_output(print.cor.str(emax1), "Correlation between time points")
      testthat::expect_output(print.cor.str(emax1), "Rho estimated from the data")
      testthat::expect_output(print.cor.str(emax1), "Median")
      testthat::expect_output(print.cor.str(emax1), "AR1")

      testthat::expect_output(print.cor.str(piece1), "Correlation between time points")
      testthat::expect_output(print.cor.str(piece1), "CS")
      testthat::expect_output(print.cor.str(piece1), "Rho assigned a numeric value\\: 0\\.2")

    })




    testthat::test_that("print.class.str functions correctly", {

      testthat::expect_output(print.class.str(emax1), NA)

      testthat::expect_output(print.class.str(piece1), "Class effects for beta\\.2")
      testthat::expect_output(print.class.str(piece1), "Random \\(exchangeable")
      testthat::expect_output(print.class.str(piece1), "Median")
      testthat::expect_output(print.class.str(piece1), "Within-class SD")
      testthat::expect_output(print.class.str(piece1), "sd\\.D\\.2")

    })





    testthat::test_that("print.modfit.str functions correctly", {

      test <- print.modfit.str(emax1)
      testthat::expect_equal(grepl("Model Fit Statistics", test), TRUE)
      testthat::expect_equal(grepl("Effective number of parameters:\n", test), TRUE)
      testthat::expect_equal(grepl("pD \\(pV\\) calculated using", test), TRUE)
      testthat::expect_equal(grepl("Deviance =", test), TRUE)
      testthat::expect_equal(grepl("Residual deviance = NOT MONITORED", test), TRUE)
      testthat::expect_equal(grepl("DIC", test), TRUE)

      test <- print.modfit.str(piece1)
      testthat::expect_equal(grepl("Model Fit Statistics", test), TRUE)
      testthat::expect_equal(grepl("Effective number of parameters:\n", test), TRUE)
      testthat::expect_equal(grepl("pV", test), FALSE)
      testthat::expect_equal(grepl("optimism", test), TRUE)
      testthat::expect_equal(grepl("Deviance =", test), TRUE)
      testthat::expect_equal(grepl("Residual deviance = NOT MONITORED", test), TRUE)
      testthat::expect_equal(grepl("DIC", test), TRUE)

    })




  }


})

