testthat::context("Testing run.functions")
data <- mb.network(osteopain)

classnetwork <- mb.network(goutSUA_CFBcomb)


################ Testing mb.run wrapped functions ################

testthat::test_that("mb.exponential function works correctly", {
  mb.result <- mb.exponential(data, lambda=list(pool="rel", method="common"),
                                    positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  expect_equal(all(c("d.lambda", "totresdev") %in% mb.result$parameters.to.save), TRUE)

  # Class effects
  mb.result <- mb.exponential(classnetwork, lambda=list(pool="rel", method="common"),
                                    positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                                    class.effect = list("lambda"="random"))
  expect_equal(all(c("D.lambda", "sd.D.lambda") %in% mb.result$parameters.to.save), TRUE)

  mb.result <- mb.exponential(classnetwork, lambda=list(pool="arm", method="random"),
                                    positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                                    class.effect = list("lambda"="random"))
  testthat::expect_equal(all(c("BETA.lambda", "sd.BETA.lambda") %in% mb.result$parameters.to.save), TRUE)

  # UME
  mb.result <- mb.exponential(data, lambda=list(pool="rel", method="common"),
                                    positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                                    UME = "lambda")
  testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("d.lambda", colnames(mb.result$BUGSoutput$sims.matrix))]),
               (28+1)*(28/2)+1)
})

testthat::test_that("mb.emax function works correctly", {
  mb.result <- mb.emax(data,
                             emax=list(pool="rel", method="common"),
                             et50=list(pool="const", method="common"),
                             positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal(all(c("d.emax", "beta.et50", "totresdev") %in% mb.result$parameters.to.save), TRUE)

  mb.result <- mb.emax(data, emax=list(pool="rel", method="common"),
                             et50=list(pool="rel", method="random"),
                             positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal(all(c("d.emax", "d.et50", "sd.et50", "totresdev") %in% mb.result$parameters.to.save), TRUE)

  result <- mb.emax(data, rho="estimate", covar="CS", n.iter=500, n.burnin=200)
  test <- all(c("d.emax", "d.et50", "rho") %in% result$parameters.to.save)
  testthat::expect_equal(test, TRUE)

  # Class effects
  mb.result <- mb.emax(classnetwork, emax=list(pool="rel", method="common"),
                             et50=list(pool="arm", method="common"),
                             positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                             class.effect = list("emax"="common"))
  testthat::expect_equal(all(c("D.emax") %in% mb.result$parameters.to.save), TRUE)
  testthat::expect_equal(all(c("sd.D.emax") %in% mb.result$parameters.to.save), FALSE)

  mb.result <- mb.emax(classnetwork, emax=list(pool="rel", method="common"),
                             et50=list(pool="rel", method="common"),
                             positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                             class.effect = list("emax"="common", "et50"="random"))
  testthat::expect_equal(all(c("D.emax", "D.et50", "sd.D.et50") %in% mb.result$parameters.to.save), TRUE)
  testthat::expect_equal(all(c("sd.D.emax") %in% mb.result$parameters.to.save), FALSE)

  mb.result <- mb.emax(classnetwork, emax=list(pool="arm", method="random"),
                             et50=list(pool="arm", method="common"),
                             positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                             class.effect = list("emax"="common", "et50"="random"))
  testthat::expect_equal(all(c("BETA.emax", "BETA.et50", "sd.BETA.et50") %in% mb.result$parameters.to.save), TRUE)
  testthat::expect_equal(all(c("sd.BETA.emax") %in% mb.result$parameters.to.save), FALSE)

  testthat::expect_error(mb.emax(classnetwork, emax=list(pool="rel", method="common"),
                                    et50=list(pool="const", method="random"),
                          positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                          class.effect = list("emax"="common", "et50"="random")))


  # UME
  mb.result <- mb.emax(data, emax=list(pool="rel", method="common"),
                             et50=list(pool="rel", method="common"),
                             positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                             UME=TRUE)
  testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("d.emax", colnames(mb.result$BUGSoutput$sims.matrix))]),
               (28+1)*(28/2)+1)
  testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("d.et50", colnames(mb.result$BUGSoutput$sims.matrix))]),
               (28+1)*(28/2)+1)

  mb.result <- mb.emax(data, emax=list(pool="rel", method="common"),
                             et50=list(pool="rel", method="common"),
                             positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                             UME=c("et50"))
  testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("d.emax", colnames(mb.result$BUGSoutput$sims.matrix))]),
               29)
  testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("d.et50", colnames(mb.result$BUGSoutput$sims.matrix))]),
               (28+1)*(28/2)+1)
})

testthat::test_that("mb.emax.hill function works correctly", {
  mb.result <- mb.emax.hill(data, emax=list(pool="rel", method="common"),
                                  et50=list(pool="rel", method="common"),
                                  hill=list(pool="const", method="common"),
                                  positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal(all(c("d.emax", "d.et50", "beta.hill", "totresdev") %in% mb.result$parameters.to.save), TRUE)


  mb.result <- mb.emax.hill(mb.network(alog_pcfb), emax=list(pool="rel", method="common"),
                                et50=list(pool="rel", method="common"),
                                hill=list(pool="const", method="common"),
                                positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal(all(c("d.emax", "d.et50", "beta.hill", "totresdev") %in% mb.result$parameters.to.save), TRUE)
})

testthat::test_that("mb.fract.first function works correctly", {
  testthat::expect_error(mb.fract.first(data,
                                           slope=list(pool="rel", method="random"),
                                           power=list(pool="rel", method="common"),
                                 positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200))

  mb.result <- mb.fract.first(data, slope=list(pool="rel", method="random"),
                                    power=list(pool="const", method="random"),
                                    positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal(all(c("d.slope", "beta.power", "sd.slope", "totresdev", "sd.beta.power") %in% mb.result$parameters.to.save), TRUE)
})

testthat::test_that("mb.fract.second function works correctly", {
  mb.result <- mb.fract.second(mb.network(alog_pcfb), slope.1=list(pool="rel", method="common"),
                                     slope.2=list(pool="arm", method="common"),
                                     power.1=list(pool="const", method="common"),
                                     power.2=list(pool="const", method="random"),
                                     positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal(all(c("d.slope.1", "beta.slope.2", "beta.power.1", "beta.power.2", "sd.beta.power.2") %in% mb.result$parameters.to.save), TRUE)
})

testthat::test_that("mb.linear function works correctly", {
  mb.result <- mb.linear(data, slope=list(pool="rel", method="common"),
                               positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal(all(c("d.slope") %in% mb.result$parameters.to.save), TRUE)
})

testthat::test_that("mb.quadratic function works correctly", {
  mb.result <- mb.quadratic(data, beta.1=list(pool="rel", method="random"),
                                  beta.2=list(pool="rel", method="common"),
                                  positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal(all(c("d.1", "sd.1", "d.2") %in% mb.result$parameters.to.save), TRUE)

  mb.result <- mb.quadratic(mb.network(alog_pcfb), beta.1=list(pool="rel", method="random"),
                                  beta.2=list(pool="const", method="common"),
                                  positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal(all(c("d.1", "sd.1", "beta.2") %in% mb.result$parameters.to.save), TRUE)
})


testthat::test_that("mb.piecelinear function works correctly", {
  mb.result <- mb.piecelinear(data, slope.1=list(pool="rel", method="common"),
                                    slope.2=list(pool="rel", method="common"),
                                    knot=list(pool="const", method="random"),
                                    alpha="arm",
                                  positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal(all(c("d.slope.1", "d.slope.2", "beta.knot", "sd.beta.knot", "totresdev") %in% mb.result$parameters.to.save), TRUE)

  mb.result <- mb.piecelinear(data, slope.1=list(pool="rel", method="random"),
                                    slope.2=list(pool="rel", method="common"),
                                    knot=list(pool="const", method=1),
                                    alpha="study",
                                    positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200)
  testthat::expect_equal("beta.knot" %in% mb.result$parameters.to.save, FALSE)

  expect_error(mb.piecelinear(data, slope.1=list(pool="rel", method="random"),
                              slope.2=list(pool="rel", method="common"),
                              knot=list(pool="const", method=1),
                              alpha="study",
                              positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                              UME=c("slope.2")), NA)

  expect_error(mb.piecelinear(data, slope.1=list(pool="rel", method="random"),
                              slope.2=list(pool="rel", method="common"),
                              knot=list(pool="const", method=1),
                              alpha="study",
                              positive.scale=TRUE,  n.chain=3, n.iter=500, n.burnin=200,
                              UME=c("notaparameter")), "do not match time-course")
})


testthat::test_that("mb.run function works correctly", {
    testthat::expect_warning(mb.run(mb.network(alog_pcfb), pd="plugin",  n.chain=3, n.iter=500, n.burnin=200), "Plugin method only works")

    testthat::expect_error(mb.run(data, pd="plugin", rho=0.5, covar="AR1",  n.chain=3, n.iter=500, n.burnin=200), "pD cannot be calculated")

    alognet <- mb.network(alog_pcfb)
    expect_error(mb.run(alognet, pd="pd.kl", n.chain=3, n.iter=500, n.burnin=200), NA)

    result <- mb.run(data, fun="linear", beta.1=list(pool="rel", method="random"),
                        rho=0.5, covar="CS",
                        n.chain=3, n.iter=500, n.burnin=200)
    test <- all(c("d.1", "sd.1", "rho") %in% result$parameters.to.save)
    testthat::expect_equal(test, TRUE)

    # Class effects
    user.fun <- ~exp(alpha + beta.1*time + beta.2 + time)
    result <- mb.run(classnetwork, fun="user", user.fun=user.fun,
                        beta.1=list(pool="rel", method="random"),
                        beta.2=list(pool="rel", method="common"),
                        class.effect=list("beta.2"="random"),
                        n.chain=3, n.iter=500, n.burnin=200)
    testthat::expect_equal(all(c("D.2", "sd.D.2") %in% result$parameters.to.save), TRUE)
    testthat::expect_equal(all(c("D.1") %in% result$parameters.to.save), FALSE)

    user.fun <- ~exp(alpha + beta.1*time + beta.2 + time)
    result <- mb.run(classnetwork, fun="user", user.fun=user.fun,
                        beta.1=list(pool="arm", method="random"),
                        beta.2=list(pool="arm", method="common"),
                        class.effect=list("beta.2"="random"),
                        n.chain=3, n.iter=500, n.burnin=200)
    testthat::expect_equal(all(c("BETA.2", "sd.BETA.2") %in% result$parameters.to.save), TRUE)
    testthat::expect_equal(all(c("BETA.1") %in% result$parameters.to.save), FALSE)

    testthat::expect_error(mb.run(classnetwork, fun="user", user.fun=user.fun,
                           beta.1=list(pool="rel", method="random"),
                           class.effect=list("beta.2"="common"),
                           n.chain=3, n.iter=500, n.burnin=200))

    # UME
    user.fun <- ~exp(alpha + beta.1*time)
    result <- mb.run(data, fun="user", user.fun=user.fun,
                        beta.1=list(pool="rel", method="random"),
                        UME=TRUE,
                        n.chain=3, n.iter=500, n.burnin=200)
    testthat::expect_equal("d.1[3,15]" %in% colnames(result$BUGSoutput$sims.matrix), TRUE)
  })





test_that("mb.update function correctly", {

  result <- mb.run(data, fun="exponential",
                     beta.1=list(pool="rel", method="random"),
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
