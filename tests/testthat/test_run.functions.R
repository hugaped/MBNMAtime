testthat::context("Testing run.functions")

datalist <- list(osteopain=osteopain, copd=copd, goutSUA_CFBcomb=goutSUA_CFBcomb,
                 hyalarthritis=hyalarthritis, diabetes=diabetes, alog_pcfb=alog_pcfb)

################ Testing mb.run wrapped functions ################

testthat::test_that("run.functions tests:", {

  testthat::expect_equal(1,1) # Avoids empty tests

  skip_on_appveyor()
  skip_on_ci()
  skip_on_cran()

  n.iter <- 500
  n.burnin <- 200

  for (i in seq_along(datalist)) {

    print(names(datalist)[i])

    network <- mb.network(datalist[[i]])


    testthat::test_that(paste0(names(datalist)[i], ": exponential time-course function works correctly"), {

      # SUPPRESSES WARNINGS FOR VERSION 0.2.2 - REMOVE AFTER THIS AND TEST WITHOUT TO ENSURE WARNINGS IDENTIFIED
      suppressWarnings({

        if (!names(datalist)[i] %in% c("goutSUA_CFBcomb", "hyalarthritis", "diabetes", "alog_pcfb")) {
          mb.result <- mb.run(network, fun=titp(pool.emax="rel", method.emax="common"),
                              positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin)
          expect_equal(all(c("emax", "totresdev") %in% mb.result$parameters.to.save), TRUE)
        }

        if ("n" %in% names(datalist[[i]])) {

          if (any(is.na(datalist[[i]]$n))) {

            expect_error(mb.run(network, link="smd", fun=titp(pool.emax="rel", method.emax="random"),
                                positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin),
                         "Missing values in n")

          } else {
            mb.result <- mb.run(network, link="smd", fun=titp(pool.emax="rel", method.emax="random"),
                                positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin)
            expect_equal(all(c("emax", "sd.emax", "totresdev") %in% mb.result$parameters.to.save), TRUE)
          }

        } else {

          expect_error(mb.run(network, link="smd", fun=titp(pool.emax="rel", method.emax="random"),
                              positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin))
        }


        # Class effects
        if ("class" %in% names(datalist[[i]])) {

          mb.result <- mb.run(network, fun=titp(pool.emax="rel", method.emax="common"),
                              positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin,
                              class.effect = list("emax"="random"))
          expect_equal(all(c("emax", "emax", "sd.EMAX") %in% mb.result$parameters.to.save), TRUE)

          mb.result <- mb.run(network, fun=titp(pool.emax="rel", method.emax="random"),
                              positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin,
                              class.effect = list("emax"="random"))
          testthat::expect_equal(all(c("EMAX", "sd.EMAX", "sd.emax") %in% mb.result$parameters.to.save), TRUE)

        }

        # UME
        mb.result <- mb.run(network, link="log", fun=titp(pool.emax="rel", method.emax="common"),
                            positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin,
                            UME = "emax")
        testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("emax", colnames(mb.result$BUGSoutput$sims.matrix))]),
                               ncol(combn(network$treatments,2))+1)

      })
    })



    testthat::test_that(paste0(names(datalist)[i], ": itp time-course function works correctly"), {

      mb.result <- suppressWarnings(mb.run(network, fun=titp(pool.emax="rel", method.emax="common",
                                                             pool.rate="abs", method.rate="common"),
                                           positive.scale=TRUE,  n.iter=n.iter, n.burnin=n.burnin, pd="plugin"))
      testthat::expect_equal(all(c("emax", "rate", "totresdev") %in% mb.result$parameters.to.save), TRUE)

      mb.result <- mb.run(network, fun=titp(pool.emax="rel", method.emax="common",
                                            pool.rate="rel", method.rate="random"),
                          positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="popt")
      testthat::expect_equal(all(c("emax", "rate", "sd.rate", "totresdev") %in% mb.result$parameters.to.save), TRUE)


      result <- mb.run(network, fun=titp(),
                       rho="dunif(0,1)", covar="CS", n.iter=n.iter, n.burnin=n.burnin, pd="pv")
      test <- all(c("emax", "rate", "rho") %in% result$parameters.to.save)
      testthat::expect_equal(test, TRUE)

      # Class effects
      if ("class" %in% names(datalist[[i]])) {
        mb.result <- mb.run(network, fun=titp(pool.emax="rel", method.emax="common",
                                              pool.rate="abs", method.rate="common"),
                            positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin,
                            class.effect = list("emax"="common"))
        testthat::expect_equal(all(c("EMAX") %in% mb.result$parameters.to.save), TRUE)
        testthat::expect_equal(all(c("sd.EMAX") %in% mb.result$parameters.to.save), FALSE)

        testthat::expect_error(mb.run(network, fun=titp(pool.emax="rel", method.emax="common",
                                                             pool.rate="abs", method.rate="common"),
                                      class.effect = list("rate"="common")),
                               "Class effects can only"
        )
      }

      # UME
      mb.result <- mb.run(network, fun=titp(pool.emax="rel", method.emax="common",
                                            pool.rate="rel", method.rate="common"),
                          positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin,
                          UME=TRUE)
      testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("emax", colnames(mb.result$BUGSoutput$sims.matrix))]),
                             ncol(combn(network$treatments,2))+1)
      testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("rate", colnames(mb.result$BUGSoutput$sims.matrix))]),
                             ncol(combn(network$treatments,2))+1)

    })


    testthat::test_that(paste0(names(datalist)[i], ": emax time-course function works correctly"), {

      # SUPPRESSES WARNINGS FOR VERSION 0.2.2 - REMOVE AFTER THIS AND TEST WITHOUT TO ENSURE WARNINGS IDENTIFIED
      suppressWarnings({

        mb.result <- suppressWarnings(mb.run(network, fun=temax(pool.emax="rel", method.emax="common",
                                                                pool.et50="abs", method.et50="common"),
                                             positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="plugin"))
        testthat::expect_equal(all(c("emax", "et50", "totresdev") %in% mb.result$parameters.to.save), TRUE)

        mb.result <- mb.run(network, fun=temax(pool.emax="rel", method.emax="common",
                                               pool.et50="rel", method.et50="random"),
                            positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="popt")
        testthat::expect_equal(all(c("emax", "et50", "sd.et50", "totresdev") %in% mb.result$parameters.to.save), TRUE)


        result <- mb.run(network, fun=temax(),
                         rho="dunif(0,1)", covar="CS", n.iter=n.iter, n.burnin=n.burnin, pd="pv")
        test <- all(c("emax", "et50", "rho") %in% result$parameters.to.save)
        testthat::expect_equal(test, TRUE)

        # Class effects
        if ("class" %in% names(datalist[[i]])) {

          mb.result <- mb.run(network, fun=temax(pool.emax="rel", method.emax="common",
                                                 pool.et50="abs", method.et50="common"),
                              positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin,
                              class.effect = list("emax"="common"))
          testthat::expect_equal(all(c("EMAX") %in% mb.result$parameters.to.save), TRUE)
          testthat::expect_equal(all(c("sd.EMAX") %in% mb.result$parameters.to.save), FALSE)

          testthat::expect_error(mb.run(network, fun=temax(pool.emax="rel", method.emax="common",
                                                           pool.et50="abs", method.et50="common"),
                                        class.effect = list("et50"="common")),
                                 "Class effects can only"
          )

        }


        # UME
        mb.result <- mb.run(network, fun=temax(pool.emax="rel", method.emax="common",
                                               pool.et50="rel", method.et50="common"),
                            positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin,
                            UME=TRUE)
        testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("emax", colnames(mb.result$BUGSoutput$sims.matrix))]),
                               ncol(combn(network$treatments,2))+1)
        testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("et50", colnames(mb.result$BUGSoutput$sims.matrix))]),
                               ncol(combn(network$treatments,2))+1)


        expect_error(mb.run(network, fun=temax(pool.emax="rel", method.emax="common",
                                               pool.et50="rel", method.et50="common"),
                            corparam = TRUE,
                            positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin,
                            UME=TRUE), NA)

        # Include a Hill parameter
        mb.result <- mb.run(network, fun=temax(pool.emax="rel", method.emax="common",
                                               pool.et50="rel", method.et50="common",
                                               pool.hill = "abs", method.hill = 2),
                            n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="pv")
        testthat::expect_equal(all(c("hill") %in% mb.result$parameters.to.save), TRUE)

      })
    })








    testthat::test_that(paste0(names(datalist)[i], ": polynomial time-course function works correctly"), {

      mb.result <- mb.run(network, fun=tpoly(degree = 2, pool.1 = "rel", method.1="common",
                                             pool.2="abs", method.2="common"),
                          n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="pv")
      testthat::expect_equal(all(c("beta.2", "d.1", "totresdev") %in% mb.result$parameters.to.save), TRUE)


      mb.result <- mb.run(network, fun=tpoly(degree = 4, pool.1 = "rel", method.1="common",
                                             pool.2="rel", method.2="common",
                                             pool.3="abs", method.3="random",
                                             pool.4="rel", method.4="random"),
                          n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="pv", rho=0.8)
      testthat::expect_equal(all(c("beta.3", "d.1", "d.2", "d.4",
                                   "sd.beta.3", "sd.beta.4",
                                   "totresdev", "rho") %in% mb.result$parameters.to.save), TRUE)


      if (names(datalist)[i] %in% c("copd", "goutSUA_CFBcomb")) {
        mb.result <- mb.run(network, fun=tpoly(degree = 4, pool.1 = "rel", method.1="common",
                                               pool.2="rel", method.2="common",
                                               pool.3="abs", method.3="random",
                                               pool.4="rel", method.4="random"),
                            corparam = FALSE,
                            n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="pv", rho=0.8)
        testthat::expect_equal(all(c("beta.3", "d.1", "d.2", "d.4",
                                     "sd.beta.3", "sd.beta.4",
                                     "totresdev", "rho") %in% mb.result$parameters.to.save), TRUE)
      } else {

        mb.result <- mb.run(network, fun=tpoly(degree = 4, pool.1 = "rel", method.1="common",
                                               pool.2="rel", method.2="common",
                                               pool.3="abs", method.3="random",
                                               pool.4="rel", method.4="random"),
                            corparam = TRUE,
                            n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="pv", rho=0.8)
        testthat::expect_equal(all(c("beta.3", "d.1", "d.2", "d.4",
                                     "sd.beta.3", "sd.beta.4",
                                     "totresdev", "rho", "rhoparam") %in% mb.result$parameters.to.save), TRUE)
      }





      # Class effects
      if ("class" %in% names(datalist[[i]])) {
        mb.result <- mb.run(network, fun=tpoly(degree = 3, pool.1 = "rel", method.1="common",
                                               pool.2="abs", method.2="common",
                                               pool.3="rel", method.3="random"),
                            positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="pv",
                            class.effect = list("beta.3"="random"))
        testthat::expect_equal(all(c("D.3") %in% mb.result$parameters.to.save), TRUE)
        testthat::expect_equal(all(c("sd.D.3") %in% mb.result$parameters.to.save), TRUE)
      }

      # UME
      mb.result <- mb.run(network, fun=tpoly(degree = 1, pool.1 = "rel", method.1="common"),
                          positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin,
                          UME=TRUE)
      testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("d.1", colnames(mb.result$BUGSoutput$sims.matrix))]),
                             ncol(combn(network$treatments,2))+1)

    })





    testthat::test_that(paste0(names(datalist)[i], ": Fractional polynomial time-course function works correctly"), {

      mb.result <- mb.run(network, fun=tfpoly(degree = 2, pool.1 = "rel", method.1="random",
                                              pool.2="abs", method.2="common"),
                          n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="pv")
      testthat::expect_equal(all(c("beta.2", "d.1", "sd.beta.1", "totresdev") %in% mb.result$parameters.to.save), TRUE)


      mb.result <- mb.run(network, fun=tpoly(degree = 4, pool.1 = "abs", method.1="common",
                                             pool.2="rel", method.2="common",
                                             pool.3="abs", method.3="random",
                                             pool.4="rel", method.4="common"),
                          n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="pv", rho="dunif(0,1)", covar = "varadj")
      testthat::expect_equal(all(c("beta.3", "beta.1", "d.2", "d.4",
                                   "sd.beta.3", "rho",
                                   "totresdev") %in% mb.result$parameters.to.save), TRUE)



      # Class effects
      if ("class" %in% names(datalist[[i]])) {
        mb.result <- mb.run(network, fun=tpoly(degree = 3, pool.1 = "rel", method.1="common",
                                               pool.2="abs", method.2="common",
                                               pool.3="rel", method.3="random"),
                            positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="pv",
                            class.effect = list("beta.3"="random"))
        testthat::expect_equal(all(c("D.3") %in% mb.result$parameters.to.save), TRUE)
        testthat::expect_equal(all(c("sd.D.3") %in% mb.result$parameters.to.save), TRUE)

        mb.result <- mb.run(network, fun=tpoly(degree = 3, pool.1 = "rel", method.1="common",
                                               pool.2="abs", method.2="common",
                                               pool.3="rel", method.3="random"),
                            corparam = TRUE,
                            positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin, pd="pv",
                            class.effect = list("beta.3"="random"))
      }


      # UME
      mb.result <- mb.run(network, fun=tpoly(degree = 1, pool.1 = "rel", method.1="common"),
                          positive.scale=TRUE,  n.chain=3, n.iter=n.iter, n.burnin=n.burnin,
                          UME=TRUE)
      testthat::expect_equal(ncol(mb.result$BUGSoutput$sims.matrix[,grepl("d.1", colnames(mb.result$BUGSoutput$sims.matrix))]),
                             ncol(combn(network$treatments,2))+1)

    })





    testthat::test_that(paste0(names(datalist)[i], ": mb.run function (+ tuser()) works correctly"), {

      testthat::expect_warning(mb.run(network, pd="plugin",  n.chain=3, n.iter=n.iter, n.burnin=n.burnin), "Plugin method only works")

      testthat::expect_error(mb.run(network, pd="plugin", rho=0.5, covar="AR1",  n.chain=3, n.iter=n.iter, n.burnin=n.burnin), "pD cannot be calculated")

      expect_error(mb.run(network, pd="pd.kl", n.chain=3, n.iter=n.iter, n.burnin=n.burnin), NA)

      # Class effects
      user.fun <- ~exp(beta.1*time + beta.2 + time)

      if ("class" %in% names(datalist[[i]])) {
        result <- mb.run(network, fun=tuser(fun=user.fun,
                                                 pool.1="rel", method.1="random",
                                                 pool.2="rel", method.2="common"),
                         class.effect=list("beta.2"="random"),
                         n.chain=3, n.iter=n.iter, n.burnin=n.burnin)
        testthat::expect_equal(all(c("D.2", "sd.D.2") %in% result$parameters.to.save), TRUE)
        testthat::expect_equal(all(c("D.1") %in% result$parameters.to.save), FALSE)

        result <- mb.run(network, fun=tuser(fun=user.fun,
                                                 pool.1="abs", method.1="random",
                                                 pool.2="rel", method.2="common"),
                         class.effect=list("beta.2"="random"),
                         n.chain=3, n.iter=n.iter, n.burnin=n.burnin)
        testthat::expect_equal(all(c("D.2", "sd.D.2") %in% result$parameters.to.save), TRUE)
        testthat::expect_equal(all(c("BETA.1") %in% result$parameters.to.save), FALSE)
        testthat::expect_equal(all(c("BETA.2") %in% result$parameters.to.save), FALSE)

        testthat::expect_error(mb.run(network, fun=tuser(fun=user.fun,
                                                              pool.1="abs", method.1="random",
                                                              pool.2="rel", method.2="common"),
                                      class.effect=list("beta.1"="common"),
                                      n.chain=3, n.iter=n.iter, n.burnin=n.burnin), "Class effects can only be specified")

      }


      # UME
      user.fun <- ~exp(beta.1*time)
      result <- mb.run(network, fun=tuser(fun=user.fun,
                                          pool.1="rel", method.1="random"),
                       UME=TRUE,
                       n.chain=3, n.iter=n.iter, n.burnin=n.burnin)

      if (length(network$treatments)>3) {
        testthat::expect_equal(paste0("d.1[2,", length(network$treatments)-1, "]") %in% colnames(result$BUGSoutput$sims.matrix), TRUE)
      }
    })






    test_that(paste0(names(datalist)[i], ": mb.update function correctly"), {

      result <- mb.run(network, fun=tloglin(method.rate="random"),
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


  }


})
