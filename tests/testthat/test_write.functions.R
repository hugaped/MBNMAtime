testthat::context("Testing write.functions")
network <- mb.network(osteopain)



testthat::test_that("testing prior writing functions", {
  skip_on_ci()
  skip_on_cran()

  emax1 <- mb.run(network,
                  fun=temax(pool.emax="rel", method.emax="random",
                            pool.et50="rel", method.et50="common"),
                  n.chain=3, n.iter=200, n.burnin=100)

  ############### Testing prior writing functions ###############

  testthat::test_that("get.prior", {
    testthat::expect_equal(sort(names(get.prior(emax1$model.arg$jagscode))), sort(c("muinv.R", "alpha", "sd.emax", "inv.R")))
    testthat::expect_equal(class(get.prior(emax1$model.arg$jagscode)), "list")
    testthat::expect_equal(class(get.prior(emax1$model.arg$jagscode)[[1]]), "character")

    priors <- get.prior(emax1$model.arg$jagscode)
    for (i in seq_along(priors)) {
      expect_equal(grepl("d[a-z]+\\(.+\\)", priors[i]), TRUE)
    }

    testthat::expect_error(get.prior(5))
  })



  testthat::test_that("replace.prior", {
    priors <- list("alpha"="dnorm(-1, 0.01)",
                   "sd.emax"="dnorm(0,0.5) T(0,)")

    testthat::expect_equal(grep("T\\(0,\\)", MBNMAtime:::replace.prior(priors, mbnma=emax1)) %in% grep("sd.emax", emax1$model.arg$jagscode),
                           TRUE)

    priors <- list("banana"="dnorm(-1,0.01)")

    testthat::expect_error(MBNMAtime:::replace.prior(priors, mbnma=emax1))

  })

})




################### Testing mb.write ################


testthat::test_that("test.mb.write", {
  funlist <- list(tpoly(degree=1, pool.1 = "rel", method.1="common"), tloglin(pool.rate="rel", method.rate="common"))

  jags <- mb.write(fun=tpoly(degree=1, pool.1 = "rel", method.1="common"), intercept=TRUE, positive.scale = TRUE)
  testthat::expect_equal(length(grep("alpha\\[i", jags))>0, TRUE)
  testthat::expect_equal(length(grep("exp\\(alpha\\[i", jags))>0, TRUE)

  jags <- mb.write(fun=tpoly(degree=1, pool.1 = "rel", method.1="common"), intercept=FALSE, positive.scale = FALSE)
  testthat::expect_equal(length(grep("alpha\\[i", jags))>0, FALSE)

  jags <- mb.write(fun=tloglin(pool.rate="rel", method.rate="random"))
  testthat::expect_equal(length(grep("rate", jags))>0, TRUE)
  testthat::expect_equal(length(grep("sd\\.rate ~ dnorm", jags))>0, TRUE)
  testthat::expect_equal(length(grep("dev", jags))>0, TRUE)

  jags <- mb.write(fun=tloglin(pool.rate="rel", method.rate="random"), rho="dunif(0,1)", covar="CS")
  testthat::expect_equal(length(grep("rho ~ dunif\\(0\\,1\\)", jags))>0, TRUE)
  testthat::expect_equal(length(grep("cor\\[i,r,c\\]", jags))>0, FALSE)
  testthat::expect_equal(length(grep("dev", jags))>0, FALSE)

  jags <- mb.write(fun=tpoly(degree=2, pool.1 = "rel", method.1="common",
                             pool.2="abs", method.2="random"),
                   rho=0.1)
  testthat::expect_equal(length(grep("rho <- 0.1", jags))>0, TRUE)

  jags <- mb.write(fun=tpoly(degree=2, pool.1 = "rel", method.1="common",
                             pool.2="abs", method.2="random"),
                   rho="dunif(-1,1)", covar="AR1")
  testthat::expect_equal(length(grep("rho ~ dunif\\(-1,1", jags))>0, TRUE)
  testthat::expect_equal(length(grep("cor\\[i,r,c\\]", jags))>0, TRUE)
  testthat::expect_equal(length(grep("dev", jags))>0, FALSE)

  testthat::expect_error(mb.write(fun=tpoly(degree=2, pool.1 = "rel", method.1="common",
                                            pool.2="abs", method.2="random"),
                                  rho=5, covar="AR1"), "cannot be outside the bounds")

  testthat::expect_error(mb.write(fun=tloglin(pool.rate="rel", method.rate="random"), class.effect="alpha"))
  testthat::expect_error(mb.write(fun=tloglin(pool.rate="rel", method.rate="random"), class.effect=list(fake="common"), "list element names in `class.effect`"))


  jags <- mb.write(fun=tfpoly(degree=2, pool.1="abs", method.1="common",
                              pool.2="rel", method.2="random", method.power1="common", method.power2="random"),
                   rho="dunif(-1,1)", covar="AR1")
  testthat::expect_equal(length(grep("beta\\.1", jags))>0, TRUE)
  testthat::expect_equal(length(grep("d\\.1", jags))>0, FALSE)
  testthat::expect_equal(length(grep("d\\.2", jags))>0, TRUE)
  testthat::expect_equal(length(grep("power1", jags))>0, TRUE)
  testthat::expect_equal(length(grep("power2", jags))>0, TRUE)
  testthat::expect_equal(length(grep("sd\\.1", jags))>0, FALSE)
  testthat::expect_equal(length(grep("sd\\.beta\\.2", jags))>0, TRUE)
  testthat::expect_equal(length(grep("sd\\.power1", jags))>0, FALSE)
  testthat::expect_equal(length(grep("sd\\.power2", jags))>0, TRUE)
  testthat::expect_equal(length(grep("rho ~ dunif\\(-1,1", jags))>0, TRUE)
  testthat::expect_equal(length(grep("cor\\[i,r,c\\]", jags))>0, TRUE)
  testthat::expect_equal(length(grep("dev", jags))>0, FALSE)


  # Class effects
  jags <- mb.write(fun=tspline(type="rcs", knots=5, pool.1 = "rel", method.1="common",
                               pool.2="rel", method.2="random", pool.3="abs", method.3="common", pool.4="abs", method.4="random"),
                   class.effect=list(beta.1="common", beta.2="random"), rho="dunif(0,1)", covar="varadj")
  testthat::expect_equal(length(grep("D\\.1", jags))>0, TRUE)
  testthat::expect_equal(length(grep("D\\.2", jags))>0, TRUE)
  testthat::expect_equal(length(grep("D\\.3", jags))>0, FALSE)
  testthat::expect_equal(length(grep("sd.D\\.1", jags))>0, FALSE)
  testthat::expect_equal(length(grep("sd.D\\.2", jags))>0, TRUE)
  testthat::expect_equal(length(grep("beta\\.3", jags))>0, TRUE)
  testthat::expect_equal(length(grep("beta\\.4", jags))>0, TRUE)
  testthat::expect_equal(length(grep("sd.beta\\.4", jags))>0, TRUE)
  testthat::expect_equal(length(grep("sd.beta\\.3", jags))>0, FALSE)


  # UME
  jags <- mb.write(fun=tspline(type="bs", degree=3, knots=1, pool.1 = "rel", method.1="common",
                               pool.2="abs", method.2="random", pool.3="rel", method.3="common", pool.4="rel", method.4="random"),
                   UME=c("beta.1", "beta.4"))
  testthat::expect_equal(length(grep("d\\.1\\[1\\,", jags))>0, TRUE)
  testthat::expect_equal(length(grep("d\\.3\\[1\\,", jags))>0, FALSE)
  testthat::expect_equal(length(grep("beta\\.2 ~", jags))>0, TRUE)
  testthat::expect_equal(length(grep("d\\.4\\[1\\,", jags))>0, TRUE)
  testthat::expect_equal(length(grep("sd.beta\\.2", jags))>0, TRUE)
  testthat::expect_equal(length(grep("sd.beta\\.1", jags))>0, FALSE)
  testthat::expect_equal(length(grep("sd.beta\\.3", jags))>0, FALSE)
  testthat::expect_equal(length(grep("sd\\.beta\\.4", jags))>0, TRUE)

  testthat::expect_error(mb.write(fun=tspline(type="bs", degree=3, knots=1, pool.1 = "rel", method.1="common",
                                              pool.2="abs", method.2="random", pool.3="rel", method.3="common", pool.4="rel", method.4="random"),
                                  UME=c("beta.1", "beta.2", "beta.3")), "can only be specified for time-course")


  testthat::expect_silent(mb.write())

})






