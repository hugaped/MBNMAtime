testthat::context("Testing write.functions")
network <- mb.network(osteopain)

testthat::test_that("write functions pass correctly:", {

  testthat::expect_equal(1,1) # Avoids empty tests

  datalist <- list(osteopain=osteopain, copd=copd, goutSUA_CFBcomb=goutSUA_CFBcomb,
                   hyalarthritis=hyalarthritis, diabetes=diabetes, alog_pcfb=alog_pcfb)

  skip_on_appveyor()
  skip_on_ci()
  skip_on_cran()

  for (i in seq_along(datalist)) {

    print(names(datalist)[i])
    network <- mb.network(datalist[[i]])


    testthat::test_that("testing prior writing functions", {

      testthat::expect_equal(1,1) # Avoids empty tests

      emax1 <- suppressWarnings(mb.run(network,
                                       fun=temax(pool.emax="rel", method.emax="random",
                                                 pool.et50="rel", method.et50="common"),
                                       corparam = TRUE,
                                       n.chain=3, n.iter=200, n.burnin=100))

      emax2 <- suppressWarnings(mb.run(network,
                                       fun=temax(pool.emax="rel", method.emax="random",
                                                 pool.et50="rel", method.et50="common"),
                                       corparam = FALSE,
                                       n.chain=3, n.iter=200, n.burnin=100))

      ############### Testing prior writing functions ###############

      testthat::test_that(paste0("get.prior for ", names(datalist)[i]), {

        vars <- c("mu.z", "sd.emax", "z", "rhoparam")
        vars2 <- c("mu.1", "mu.2", "sd.emax", "et50", "emax")

        if (!names(datalist)[i] %in% c("goutSUA_CFBcomb", "hyalarthritis", "diabetes", "copd", "alog_pcfb")) {
          vars <- append(vars, "alpha")
          vars2 <- append(vars2, "alpha")
        }
        testthat::expect_equal(sort(names(get.prior(emax1$model.arg$jagscode))), sort(vars))

        testthat::expect_equal(class(get.prior(emax1$model.arg$jagscode)), "list")
        testthat::expect_equal(class(get.prior(emax1$model.arg$jagscode)[[1]]), "character")

        testthat::expect_equal(sort(names(get.prior(emax2$model.arg$jagscode))), sort(vars2))

        priors <- get.prior(emax1$model.arg$jagscode)
        for (i in seq_along(priors)) {
          expect_equal(grepl("d[a-z]+\\(.+\\)", priors[i]), TRUE)
        }

        testthat::expect_error(get.prior(5))
      })



      testthat::test_that("replace.prior", {
        priors <- list("z"="dnorm(-1, 0.01)",
                       "sd.emax"="dnorm(0,0.5) T(0,)")

        testthat::expect_equal(grep("T\\(0,\\)", MBNMAtime:::replace.prior(priors, mbnma=emax1)) %in% grep("sd.emax", emax1$model.arg$jagscode),
                               TRUE)

        priors <- list("banana"="dnorm(-1,0.01)")

        testthat::expect_error(MBNMAtime:::replace.prior(priors, mbnma=emax1))

      })

    })

  }


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
    testthat::expect_equal(length(grep("rho ~ dunif\\(-?[0-9],1", jags))>0, TRUE)
    testthat::expect_equal(length(grep("cor\\[i,r,c\\]", jags))>0, TRUE)
    testthat::expect_equal(length(grep("dev", jags))>0, FALSE)

    testthat::expect_error(mb.write(fun=tpoly(degree=2, pool.1 = "rel", method.1="common",
                                              pool.2="abs", method.2="random"),
                                    rho=5, covar="AR1"), "cannot be outside the bounds")

    testthat::expect_error(mb.write(fun=tloglin(pool.rate="rel", method.rate="random"), class.effect="alpha"))
    testthat::expect_error(mb.write(fun=tloglin(pool.rate="rel", method.rate="random"), class.effect=list(fake="common"), "list element names in `class.effect`"))


    jags <- mb.write(fun=tfpoly(degree=2, pool.1="abs", method.1="common",
                                pool.2="rel", method.2="random", method.power1=1, method.power2=-1),
                     rho="dunif(-1,1)", covar="AR1")
    testthat::expect_equal(length(grep("beta\\.1", jags))>0, TRUE)
    testthat::expect_equal(length(grep("d\\.1", jags))>0, FALSE)
    testthat::expect_equal(length(grep("d\\.2", jags))>0, TRUE)
    testthat::expect_equal(length(grep("power1", jags))>0, TRUE)
    testthat::expect_equal(length(grep("power2", jags))>0, TRUE)
    testthat::expect_equal(length(grep("sd\\.1", jags))>0, FALSE)
    testthat::expect_equal(length(grep("sd\\.beta\\.2", jags))>0, TRUE)
    testthat::expect_equal(length(grep("sd\\.power1", jags))>0, FALSE)
    testthat::expect_equal(length(grep("sd\\.power2", jags))>0, FALSE)
    testthat::expect_equal(length(grep("rho ~ dunif\\(-?[0-9],1", jags))>0, TRUE)
    testthat::expect_equal(length(grep("cor\\[i,r,c\\]", jags))>0, TRUE)
    testthat::expect_equal(length(grep("dev", jags))>0, FALSE)


    # Class effects
    jags <- mb.write(fun=tspline(type="ns", knots=3, pool.1 = "rel", method.1="common",
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

})

