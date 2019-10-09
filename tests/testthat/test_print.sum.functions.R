testthat::context("Testing predict.functions")
network <- mb.network(osteopain)

emax1 <- mb.emax(network,
                    emax=list(pool="arm", method="random"),
                    et50=list(pool="rel", method="common"),
                    positive.scale=TRUE,
                    n.chain=3, n.iter=200, n.burnin=100,
                    rho="estimate", covar="AR1",
                    UME=TRUE, priors=list("rho"="dunif(0,1)"))

emax2 <- mb.run(network, fun="emax",
                   beta.1=list(pool="arm", method="random"),
                   beta.2=list(pool="rel", method="common"),
                   positive.scale=TRUE,
                   n.chain=3, n.iter=200, n.burnin=100,
                   rho="estimate", covar="AR1",
                   UME=TRUE, priors=list("rho"="dunif(0,1)"))

emax.mods <- list(emax1, emax2)


network.gout <- mb.network(goutSUA_CFBcomb)
piece1 <- mb.piecelinear(network.gout,
                    slope.1=list(pool="arm", method="random"),
                    slope.2=list(pool="rel", method="common"),
                    knot=list(pool="const", method=0.6),
                    positive.scale=FALSE, intercept = FALSE,
                    n.chain=3, n.iter=200, n.burnin=100,
                    rho=0.2, covar="CS",
                    pd="popt",
                    class.effect=list("slope.2"="random"))

piece2 <- mb.run(network.gout, fun="piecelinear",
                    beta.1=list(pool="arm", method="random"),
                    beta.2=list(pool="rel", method="common"),
                    beta.3=list(pool="const", method=0.6),
                    positive.scale=FALSE, intercept = FALSE,
                    n.chain=3, n.iter=200, n.burnin=100,
                    rho=0.2, covar="CS",
                    pd="popt",
                    class.effect=list("beta.2"="random"))

piece.mods <- list(piece1, piece2)


##############################################
###########    TESTS    ######################
##############################################

testthat::test_that("print.overall.str functions correctly", {

  for (i in seq_along(emax.mods)) {
    test <- print.overall.str(emax.mods[[i]])
    testthat::expect_equal(grepl("Data modelled with intercept", test), TRUE)
    testthat::expect_equal(grepl("Time-course function: emax", test), TRUE)
    testthat::expect_equal(grepl("Responses restricted to taking positive values", test), TRUE)
    #testthat::expect_equal(grepl("Parameter modelled on exponential scale", test), TRUE)
  }

  for (i in seq_along(piece.mods)) {
    test <- print.overall.str(piece.mods[[i]])
    testthat::expect_equal(grepl("Data modelled without intercept", test), TRUE)
    testthat::expect_equal(grepl("Time-course function: piecelinear", test), TRUE)
    testthat::expect_equal(grepl("Responses restricted to taking positive values", test), FALSE)
  }

})

testthat::test_that("print.treat.str functions correctly", {

  # Emax
  for (i in seq_along(emax.mods)) {
    test <- print.treat.str(emax.mods[[i]])
    testthat::expect_equal(grepl("Parameter modelled on exponential scale", test), TRUE)
    testthat::expect_equal(grepl("Unrelated Mean Effects modelled", test), TRUE)
    testthat::expect_equal(grepl("DATA TOO LONG FOR SUMMARY", test), TRUE)
    testthat::expect_equal(grepl("Pooling: arm-based", test), TRUE)
    testthat::expect_equal(grepl("Method: random effects", test), TRUE)
    testthat::expect_equal(grepl("# Between-study SD", test), TRUE)
    testthat::expect_equal(grepl("Pooling: relative effects", test), TRUE)
    testthat::expect_equal(grepl("Method: common \\(fixed\\) effects", test), TRUE)
  }

  test <- print.treat.str(emax1)
  testthat::expect_equal(grepl("beta\\.emax\\[1\\]", test), TRUE)
  testthat::expect_equal(grepl("d\\.emax\\[2\\]", test), FALSE)
  testthat::expect_equal(grepl("sd\\.beta\\.emax", test), TRUE)
  testthat::expect_equal(grepl("#### et50 time-course parameter pooling ####", test), TRUE)

  test <- print.treat.str(emax2)
  testthat::expect_equal(grepl("beta\\.1\\[1\\]", test), TRUE)
  testthat::expect_equal(grepl("d\\.1\\[2\\]", test), FALSE)
  testthat::expect_equal(grepl("sd\\.beta\\.1", test), TRUE)
  testthat::expect_equal(grepl("#### beta.2 time-course parameter pooling ####", test), TRUE)


  # Piecelinear
  for (i in seq_along(piece.mods)) {
    test <- print.treat.str(piece.mods[[i]])
    testthat::expect_equal(grepl("Unrelated Mean Effects modelled", test), FALSE)
    testthat::expect_equal(grepl("DATA TOO LONG FOR SUMMARY", test), FALSE)
    testthat::expect_equal(grepl("CLASS EFFECTS MODELLED", test), TRUE)
    testthat::expect_equal(grepl("Pooling: arm-based", test), TRUE)
    testthat::expect_equal(grepl("Method: random effects", test), TRUE)
    testthat::expect_equal(grepl("# Between-study SD", test), TRUE)
    testthat::expect_equal(grepl("Pooling: relative effects", test), TRUE)
    testthat::expect_equal(grepl("Method: common \\(fixed\\) effects", test), TRUE)
    testthat::expect_equal(grepl("Pooling: constant \\(single parameter", test), TRUE)
    testthat::expect_equal(grepl("Assigned a numeric value: 0.6", test), TRUE)
  }

  test <- print.treat.str(piece1)
  testthat::expect_equal(grepl("beta\\.slope\\.1\\[1\\]", test), TRUE)
  testthat::expect_equal(grepl("d\\.slope.1\\[2\\]", test), FALSE)
  testthat::expect_equal(grepl("Network reference", test), FALSE)
  testthat::expect_equal(grepl("d\\.slope.2\\[2\\]", test), FALSE)
  testthat::expect_equal(grepl("sd\\.beta\\.slope\\.1", test), TRUE)
  testthat::expect_equal(grepl("#### knot time-course parameter pooling ####", test), TRUE)

  test <- print.treat.str(piece2)
  testthat::expect_equal(grepl("beta\\.1\\[1\\]", test), TRUE)
  testthat::expect_equal(grepl("beta\\.2\\[2\\]", test), FALSE)
  testthat::expect_equal(grepl("Network reference", test), FALSE)
  testthat::expect_equal(grepl("d\\.2\\[2\\]", test), FALSE)
  testthat::expect_equal(grepl("sd\\.beta\\.1", test), TRUE)
  testthat::expect_equal(grepl("#### beta.3 time-course parameter pooling ####", test), TRUE)

})





testthat::test_that("print.cor.str functions correctly", {
  for (i in seq_along(emax.mods)) {
    test <- print.cor.str(emax.mods[[i]])
    testthat::expect_equal(grepl("Correlation between time points", test), TRUE)
    testthat::expect_equal(grepl("AR1", test), TRUE)
    testthat::expect_equal(grepl("Rho estimated from the data", test), TRUE)
    testthat::expect_equal(grepl("rho", test), TRUE)
    testthat::expect_equal(grepl("Median \\(95%CrI\\)", test), TRUE)
  }

  for (i in seq_along(piece.mods)) {
    test <- print.cor.str(piece.mods[[i]])
    testthat::expect_equal(grepl("Correlation between time points", test), TRUE)
    testthat::expect_equal(grepl("CS", test), TRUE)
    testthat::expect_equal(grepl("Rho estimated from the data", test), FALSE)
    testthat::expect_equal(grepl("Rho assigned a numeric value: 0.2", test), TRUE)
    testthat::expect_equal(grepl("rho", test), FALSE)
    testthat::expect_equal(grepl("Median \\(95%CrI\\)", test), FALSE)
  }
})




testthat::test_that("print.class.str functions correctly", {
  for (i in seq_along(emax.mods)) {
    test <- print.class.str(emax.mods[[i]])
    testthat::expect_equal(is.null(test), TRUE)
  }

  for (i in seq_along(piece.mods)) {
    test <- print.class.str(piece.mods[[i]])
    testthat::expect_equal(grepl("Class effects ###", test), TRUE)
    testthat::expect_equal(grepl("Median \\(95%CrI\\)", test), TRUE)
    testthat::expect_equal(grepl("Within-class SD", test), TRUE)
    testthat::expect_equal(grepl("Network reference", test), TRUE)
  }

  test <- print.class.str(piece1)
  testthat::expect_equal(grepl("D\\.slope\\.2", test), TRUE)
  testthat::expect_equal(grepl("Class effect on slope\\.2", test), TRUE)
  testthat::expect_equal(grepl("sd\\.D\\.slope\\.2", test), TRUE)

  test <- print.class.str(piece2)
  testthat::expect_equal(grepl("D\\.2", test), TRUE)
  testthat::expect_equal(grepl("Class effect on beta\\.2", test), TRUE)
  testthat::expect_equal(grepl("sd\\.D\\.2", test), TRUE)
})





testthat::test_that("print.modfit.str functions correctly", {
  for (i in seq_along(emax.mods)) {
    test <- print.modfit.str(emax.mods[[i]])
    testthat::expect_equal(grepl("Model Fit Statistics", test), TRUE)
    testthat::expect_equal(grepl("Effective number of parameters:\n", test), TRUE)
    testthat::expect_equal(grepl("pD \\(pV\\) calculated using", test), TRUE)
    testthat::expect_equal(grepl("Deviance =", test), TRUE)
    testthat::expect_equal(grepl("Residual deviance = NOT MONITORED", test), TRUE)
    testthat::expect_equal(grepl("DIC", test), TRUE)
  }

  for (i in seq_along(piece.mods)) {
    test <- print.modfit.str(piece.mods[[i]])
    testthat::expect_equal(grepl("Model Fit Statistics", test), TRUE)
    testthat::expect_equal(grepl("Effective number of parameters:\n", test), TRUE)
    testthat::expect_equal(grepl("pV", test), FALSE)
    testthat::expect_equal(grepl("optimism", test), TRUE)
    testthat::expect_equal(grepl("Deviance =", test), TRUE)
    testthat::expect_equal(grepl("Residual deviance = NOT MONITORED", test), TRUE)
    testthat::expect_equal(grepl("DIC", test), TRUE)
  }

})
