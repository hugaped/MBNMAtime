testthat::context("Testing run.functions")
painnet <- mb.network(osteopain)
copdnet <- mb.network(copd)

classnetwork <- mb.network(goutSUA_CFBcomb)


################ Testing mb.run wrapped functions ################

testthat::test_that("exponential time-course function works correctly", {
  skip_on_ci()
  skip_on_cran()

  mb.result <- mb.run(painnet, fun=texp(pool.rate="rel", method.rate="common"),
                              positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  expect_equal(all(c("rate", "totresdev") %in% mb.result$parameters.to.save), TRUE)

  mb.result <- mb.run(copdnet, link="smd", fun=texp(pool.rate="rel", method.rate="random"),
                      positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  expect_equal(all(c("rate", "sd.rate", "totresdev") %in% mb.result$parameters.to.save), TRUE)

  # Class effects
  mb.result <- mb.run(classnetwork, fun=texp(pool.rate="rel", method.rate="common"),
                              positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                      class.effect = list("rate"="random"))
  expect_equal(all(c("RATE", "rate", "sd.RATE") %in% mb.result$parameters.to.save), TRUE)

  mb.result <- mb.run(classnetwork, fun=texp(pool.rate="rel", method.rate="random"),
                      positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                      class.effect = list("rate"="random"))
  testthat::expect_equal(all(c("RATE", "sd.RATE", "sd.rate") %in% mb.result$parameters.to.save), TRUE)

  # UME
  mb.result <- mb.run(copdnet, link="log", fun=texp(pool.rate="rel", method.rate="common"),
                      positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                      UME = "rate")
  testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("rate", colnames(mb.result$BUGSoutput$sims.matrix))]),
                         4)

})


testthat::test_that("emax time-course function works correctly", {
  skip_on_ci()
  skip_on_cran()

  mb.result <- suppressWarnings(mb.run(painnet, fun=temax(pool.emax="rel", method.emax="common",
                                         pool.et50="abs", method.et50="common"),
                      positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200, pd="plugin"))
  testthat::expect_equal(all(c("emax", "et50", "totresdev") %in% mb.result$parameters.to.save), TRUE)

  mb.result <- mb.run(copdnet, fun=temax(pool.emax="rel", method.emax="common",
                                         pool.et50="rel", method.et50="random"),
                      positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200, pd="popt")
  testthat::expect_equal(all(c("emax", "et50", "sd.et50", "totresdev") %in% mb.result$parameters.to.save), TRUE)


  result <- mb.run(copdnet, fun=temax(),
                   rho="dunif(0,1)", covar="CS", n.iter=500, n.burnin=200, pd="pv")
  test <- all(c("emax", "et50", "rho") %in% result$parameters.to.save)
  testthat::expect_equal(test, TRUE)

  # Class effects
  mb.result <- mb.run(classnetwork, fun=temax(pool.emax="rel", method.emax="common",
                                              pool.et50="abs", method.et50="common"),
                      positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                      class.effect = list("emax"="common"))
  testthat::expect_equal(all(c("EMAX") %in% mb.result$parameters.to.save), TRUE)
  testthat::expect_equal(all(c("sd.EMAX") %in% mb.result$parameters.to.save), FALSE)

  testthat::expect_error(mb.run(classnetwork, fun=temax(pool.emax="rel", method.emax="common",
                                                        pool.et50="abs", method.et50="common"),
                                class.effect = list("et50"="common")),
                         "Class effects can only"
  )

  # UME
  mb.result <- mb.run(copdnet, fun=temax(pool.emax="rel", method.emax="common",
                                         pool.et50="rel", method.et50="common"),
                      positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                      UME=TRUE)
  testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("emax", colnames(mb.result$BUGSoutput$sims.matrix))]),
                         4)
  testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("et50", colnames(mb.result$BUGSoutput$sims.matrix))]),
                         4)


  # Include a Hill parameter
  mb.result <- mb.run(painnet, fun=temax(pool.emax="rel", method.emax="common",
                                         pool.et50="rel", method.et50="common",
                                         pool.hill = "abs", method.hill = "common"),
                      n.chain=3, n.iter=500, n.burnin=200, pd="pv", priors = list(hill="dnorm(0,0.1) T(-0.5,0.5)"))
  testthat::expect_equal(all(c("hill") %in% mb.result$parameters.to.save), TRUE)

})








testthat::test_that("polynomial time-course function works correctly", {
  skip_on_ci()
  skip_on_cran()

  mb.result <- mb.run(painnet, fun=tpoly(degree = 2, pool.1 = "rel", method.1="common",
                                         pool.2="abs", method.2="common"),
                      n.chain=3, n.iter=500, n.burnin=200, pd="pv")
  testthat::expect_equal(all(c("beta.2", "d.1", "totresdev") %in% mb.result$parameters.to.save), TRUE)


  mb.result <- mb.run(painnet, fun=tpoly(degree = 4, pool.1 = "rel", method.1="common",
                                         pool.2="rel", method.2="common",
                                         pool.3="abs", method.3="random",
                                         pool.4="rel", method.4="random"),
                      n.chain=3, n.iter=500, n.burnin=200, pd="pv", rho=0.8)
  testthat::expect_equal(all(c("beta.3", "d.1", "d.2", "d.4",
                               "sd.beta.3", "sd.beta.4",
                               "totresdev", "rho") %in% mb.result$parameters.to.save), TRUE)



  # Class effects
  mb.result <- mb.run(classnetwork, fun=tpoly(degree = 3, pool.1 = "rel", method.1="common",
                                              pool.2="abs", method.2="common",
                                              pool.3="rel", method.3="random"),
                      positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200, pd="pv",
                      class.effect = list("beta.3"="random"))
  testthat::expect_equal(all(c("D.3") %in% mb.result$parameters.to.save), TRUE)
  testthat::expect_equal(all(c("sd.D.3") %in% mb.result$parameters.to.save), TRUE)

  # UME
  mb.result <- mb.run(copdnet, fun=tpoly(degree = 1, pool.1 = "rel", method.1="common"),
                      positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                      UME=TRUE)
  testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("d.1", colnames(mb.result$BUGSoutput$sims.matrix))]),
                         4)

})





testthat::test_that("Fractional polynomial time-course function works correctly", {
  skip_on_ci()
  skip_on_cran()

  mb.result <- mb.run(painnet, fun=tfpoly(degree = 2, pool.1 = "rel", method.1="random",
                                          pool.2="abs", method.2="common"),
                      n.chain=3, n.iter=500, n.burnin=200, pd="pv")
  testthat::expect_equal(all(c("beta.2", "d.1", "sd.beta.1", "totresdev") %in% mb.result$parameters.to.save), TRUE)


  mb.result <- mb.run(copdnet, fun=tpoly(degree = 4, pool.1 = "abs", method.1="common",
                                         pool.2="rel", method.2="common",
                                         pool.3="abs", method.3="random",
                                         pool.4="rel", method.4="common"),
                      n.chain=3, n.iter=500, n.burnin=200, pd="pv", rho="dunif(0,1)", covar = "varadj")
  testthat::expect_equal(all(c("beta.3", "beta.1", "d.2", "d.4",
                               "sd.beta.3", "rho",
                               "totresdev") %in% mb.result$parameters.to.save), TRUE)



  # Class effects
  mb.result <- mb.run(classnetwork, fun=tpoly(degree = 3, pool.1 = "rel", method.1="common",
                                              pool.2="abs", method.2="common",
                                              pool.3="rel", method.3="random"),
                      positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200, pd="pv",
                      class.effect = list("beta.3"="random"))
  testthat::expect_equal(all(c("D.3") %in% mb.result$parameters.to.save), TRUE)
  testthat::expect_equal(all(c("sd.D.3") %in% mb.result$parameters.to.save), TRUE)

  # UME
  mb.result <- mb.run(copdnet, fun=tpoly(degree = 1, pool.1 = "rel", method.1="common"),
                      positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                      UME=TRUE)
  testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("d.1", colnames(mb.result$BUGSoutput$sims.matrix))]),
                         4)

})





testthat::test_that("mb.run function (+ tuser()) works correctly", {
  skip_on_ci()
  skip_on_cran()

  testthat::expect_warning(mb.run(mb.network(alog_pcfb), pd="plugin",  n.chain=3, n.iter=500, n.burnin=200), "Plugin method only works")

  testthat::expect_error(mb.run(mb.network(copd), pd="plugin", rho=0.5, covar="AR1",  n.chain=3, n.iter=500, n.burnin=200), "pD cannot be calculated")

  alognet <- mb.network(alog_pcfb)
  expect_error(mb.run(alognet, pd="pd.kl", n.chain=3, n.iter=500, n.burnin=200), NA)

  # Class effects
  user.fun <- ~exp(beta.1*time + beta.2 + time)
  result <- mb.run(classnetwork, fun=tuser(fun=user.fun,
                                           pool.1="rel", method.1="random",
                                           pool.2="rel", method.2="common"),
                   class.effect=list("beta.2"="random"),
                   n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal(all(c("D.2", "sd.D.2") %in% result$parameters.to.save), TRUE)
  testthat::expect_equal(all(c("D.1") %in% result$parameters.to.save), FALSE)

  result <- mb.run(classnetwork, fun=tuser(fun=user.fun,
                                           pool.1="abs", method.1="random",
                                           pool.2="rel", method.2="common"),
                   class.effect=list("beta.2"="random"),
                   n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal(all(c("D.2", "sd.D.2") %in% result$parameters.to.save), TRUE)
  testthat::expect_equal(all(c("BETA.1") %in% result$parameters.to.save), FALSE)
  testthat::expect_equal(all(c("BETA.2") %in% result$parameters.to.save), FALSE)

  testthat::expect_error(mb.run(classnetwork, fun=tuser(fun=user.fun,
                                                        pool.1="abs", method.1="random",
                                                        pool.2="rel", method.2="common"),
                                class.effect=list("beta.1"="common"),
                                n.chain=3, n.iter=500, n.burnin=200), "Class effects can only be specified")

  # UME
  user.fun <- ~exp(beta.1*time)
  result <- mb.run(painnet, fun=tuser(fun=user.fun,
                                      pool.1="rel", method.1="random"),
                   UME=TRUE,
                   n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal("d.1[3,15]" %in% colnames(result$BUGSoutput$sims.matrix), TRUE)
})






test_that("mb.update function correctly", {
  skip_on_ci()
  skip_on_cran()

  result <- mb.run(copdnet, fun=tloglin(method.rate="random"),
                   UME=TRUE,
                   n.chain=3, n.iter=500, n.burnin=200)

  expect_error(mb.update(result, param="test"))

  update <- mb.update(result, param="resdev")
  expect_equal(sort(names(update)), sort(c("study", "arm", "mean", "fup")))

  update <- mb.update(result, param="theta")
  expect_equal(sort(names(update)), sort(c("study", "arm", "mean", "fup")))

  update <- mb.update(result, param="dev")
  expect_equal(sort(names(update)), sort(c("study", "arm", "mean", "fup")))

})
