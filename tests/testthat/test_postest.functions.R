testthat::context("Testing predict.functions")

datalist <- list(osteopain=osteopain, copd=copd, goutSUA_CFBcomb=goutSUA_CFBcomb,
                 hyalarthritis=hyalarthritis, diabetes=diabetes, alog_pcfb=alog_pcfb)

seed <- 890421

testthat::test_that("post-estimation tests pass correctly", {

  skip_on_appveyor()
  skip_on_ci()
  skip_on_cran()

  testthat::expect_equal(1,1) # Avoids empty test

  for (i in seq_along(datalist)) {

    print(names(datalist)[i])
    network <- mb.network(datalist[[i]])

    testthat::test_that(paste0(names(datalist)[i], ": predict.functions tests pass correctly"), {

      testthat::expect_equal(1,1) # Avoids empty test

      loglin <- mb.run(network, fun=tloglin(pool.rate="rel", method.rate="common"), jags.seed=seed)


      # SUPPRESSES WARNINGS FOR VERSION 0.2.2 - REMOVE AFTER THIS AND TEST WITHOUT TO ENSURE WARNINGS IDENTIFIED
      suppressWarnings({

        emax <- mb.run(network, fun=temax(pool.emax="rel", method.emax="random",
                                          pool.et50="abs", method.et50="common"), pd="pv", jags.seed=seed)

      })


      bs <- mb.run(network, fun=tspline(type="bs", degree=1, knots=3,
                                        pool.1="rel", method.1="common",
                                        pool.2="abs", method.2="random",
                                        pool.3 = "rel", method.3="random"), pd="pv", jags.seed=seed)


      ls <- mb.run(network, fun=tspline(type="ls", knots = 25/250),
                   rho="dunif(0,1)", covar="varadj", pd="pv", jags.seed=seed)

      loglin.ar1 <- mb.run(network, fun=tloglin(pool.rate="rel", method.rate="common"), covar="AR1",
                           rho="dunif(0,1)", n.iter=1500, pd="pv", jags.seed=seed)

      resdev <- mb.run(network, fun=tpoly(degree=1), parameters.to.save = "resdev", n.iter=1000, pd="pv", jags.seed=seed)


      model.list <- list(loglin, emax, bs, ls, loglin.ar1)



      testthat::test_that(paste0(names(datalist)[i], ": predict.mbnma functions correctly"), {



        treats.list <- list(network$treatments[2:3],
                            NULL)
        if (length(network$treatments)>5) {
          treats.list[[length(treats.list)+1]] <- c(1,3,4)
          treats.list[[length(treats.list)+1]] <- network$treatments[2:4]
        }
        treats.list[[length(treats.list)+1]] <- NULL
        treats.list[[length(treats.list)+1]] <- NULL
        treats.list[[length(treats.list)+1]] <- network$treatments[1:3]
        treats.list[[length(treats.list)+1]] <- network$treatments[c(1,3)]



        ref.resp.list <- list(network$data.ab[network$data.ab$treatment==1,],
                              network$data.ab[network$data.ab$treatment==2,],
                              list(beta.1=~rnorm(n,-0.1,0.1), beta.3=~rnorm(n, 0.2, 0.01), beta.4=~rnorm(n,0,0)),
                              list(beta.1=-1, beta.2=0.1),
                              NULL,
                              NULL
        )

        times.list <- list(c(0:10),
                           c(1,10:20),
                           seq(0, max(bs$model.arg$jagsdata$time, na.rm=TRUE), length.out=20),
                           c(2,4,6,7,10),
                           seq(0, max(ls$model.arg$jagsdata$time, na.rm=TRUE), length.out=20),
                           c(0:20)
        )

        E0.list <- list(7,
                        ~rnorm(n, 7,2),
                        0,
                        ~rnorm(n,5,5),
                        10,
                        0)
        synth.list <- rep(c("common", "random"),3)

        for (k in 1:6) {
          print(paste0("modellist: ", k))

          if (length(model.list)<k) {
            mbnma <- model.list[[k-1]]
          } else {
            mbnma <- model.list[[k]]
          }
          E0 <- E0.list[[k]]
          ref.resp <- ref.resp.list[[k]]

          if (length(mbnma$model.arg$class.effect)>0) {
            treats <- network$classes[1:2]
          } else {
            if (length(treats.list)<k) {
              treats <- NULL
            } else {
              treats <- treats.list[[k]]
            }
          }

          times <- times.list[[k]]
          synth <- synth.list[[k]]

          # Tests using ref.resp
          pred <- suppressWarnings(predict(mbnma, times=times,
                                           level=ifelse(length(mbnma$model.arg$class.effect)>0,
                                                        "class", "treatment"),
                                           E0=E0, treats=treats,
                                           ref.resp=ref.resp, synth=synth))


          if (!is.null(treats)) {
            testthat::expect_equal(length(pred$pred.mat), length(treats))
          } else {
            testthat::expect_equal(length(pred$pred.mat), length(mbnma$network$treatments))
          }

          if (is.numeric(treats)) {
            testthat::expect_equal(names(pred$pred.mat), mbnma$network$treatments[treats])
          } else if (is.character(treats)) {
            testthat::expect_equal(names(pred$pred.mat), treats)
          }

          testthat::expect_identical(names(pred), c("summary", "pred.mat", "network", "times", "link", "lim"))
          testthat::expect_equal(nrow(pred$pred.mat[[1]]), mbnma$BUGSoutput$n.sims)
          testthat::expect_equal(nrow(pred$summary[[1]]), length(times))
          testthat::expect_equal(identical(pred$summary[[1]]$time, times), TRUE)
        }

        # Tests of class models
        if ("class" %in% names(datalist[[i]])) {

          class <- mb.run(network, fun=tfpoly(degree=2, method.1="common", pool.1="rel",
                                              method.2="random", pool.2="rel"),
                          class.effect = list(beta.1="random"), pd="pv", jags.seed=seed
          )

          expect_error(predict(class, level="class"), "all relative effects must be modelled with class")

          class2 <- mb.run(network, fun=tpoly(degree=2, method.1="common", pool.1="rel",
                                              method.2="random", pool.2="rel"),
                           class.effect = list(beta.1="random", beta.2="common"), pd="pv", jags.seed=seed
          )

          pred <- predict(class2, level="class")
          testthat::expect_equal(names(pred$pred.mat), class2$network$classes)

          pred <- predict(class2, level="class", treats = class2$network$classes[2:3])
          testthat::expect_equal(names(pred$pred.mat), class2$network$classes[2:3])

          pred <- predict(class2, level="class", treats = c(1,3,5))
          testthat::expect_equal(names(pred$pred.mat), class2$network$classes[c(1,3,5)])

        }


        # Tests using ref.resp
        ref.resp <- list("emax"=-1)
        testthat::expect_error(predict(emax,
                                       E0=7,
                                       ref.resp=ref.resp),
                               NA)

        ref.resp <- list("d.emax"=-1) # incorrect prior name ("d.emax" rather than "emax")
        testthat::expect_error(predict(emax,
                                       E0=7,
                                       ref.resp=ref.resp), "Named elements of")

        ref.resp <- list("rate"=~rnorm(n, -1,1))
        testthat::expect_error(predict(loglin,
                                       E0=E0,
                                       ref.resp=ref.resp), NA)

        ref.resp <- list("rate"="rnorm(n, -1,1)", "beta.2"="rnorm(n, 1, 0.1)") # beta.2 is not a relative effect in quad
        testthat::expect_error(predict(loglin,
                                       E0=E0,
                                       ref.resp=ref.resp), "Named elements of")


        # Error due to wrong parameters being saved from model
        testthat::expect_error(predict(resdev), "Parameters required for estimation of time-course")


        # Expect no error even if ref.resp is NULL
        testthat::expect_error(predict(bs, ref.resp=NULL), NA)

        # Expect no error even when only a single time point is predicted
        testthat::expect_error(predict(bs, times=2), NA)

      })



      testthat::test_that(paste0(names(datalist)[i], ": ref.synth functions correctly"), {
        ref.resp <- datalist[[i]][datalist[[i]]$treatment==datalist[[i]]$treatment[1], ]

        testthat::expect_error(suppressWarnings(ref.synth(ref.resp, emax, synth="random",
                                                          n.burnin=100, n.iter=200)), NA)

      })




      testthat::test_that(paste0(names(datalist)[i], ": get.model.vals functions correctly"), {

        vals <- get.model.vals(loglin)
        testthat::expect_equal(names(vals), c("alpha", "d.1", "timecourse", "time.params"))

        vals <- get.model.vals(emax)
        testthat::expect_equal(names(vals), c("alpha", "d.1", "beta.2", "timecourse", "time.params"))

        # Can add more here if problems
      })
    })
  }

})


