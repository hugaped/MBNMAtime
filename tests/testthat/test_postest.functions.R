testthat::context("Testing predict.functions")
painnet <- mb.network(osteopain)
alognet <- mb.network(alog_pcfb)
copdnet <- mb.network(copd)
goutnet <- mb.network(goutSUA_CFBcomb)
obesenet <- mb.network(obesityBW_CFB)


testthat::test_that("predict.functions tests pass correctly", {
  skip_on_ci()
  skip_on_cran()

  loglin <- mb.run(painnet, fun=tloglin(pool.rate="rel", method.rate="common"))

  emax <- mb.run(alognet, fun=temax(pool.emax="rel", method.emax="random",
                                    pool.et50="abs", method.et50="common"), pd="pv")




  ################### Testing add_index ################

  testthat::test_that("predict.mbnma functions correctly", {
    skip_on_ci()
    skip_on_cran()

    # Run models
    bs <- mb.run(copdnet, fun=tspline(type="bs", degree=1, knots=3,
                                      pool.1="rel", method.1="common",
                                      pool.2="abs", method.2="random",
                                      pool.3 = "rel", method.3="random"), pd="pv")

    class <- mb.run(goutnet, fun=tfpoly(degree=2, method.1="common", pool.1="rel",
                                        method.2="random", pool.2="rel"),
                    class.effect = list(beta.1="random"), pd="pv"
    )

    ls <- mb.run(obesenet, fun=tspline(type="ls", knots = 25/250),
                 rho="dunif(0,1)", covar="varadj", pd="pv")

    loglin.ar1 <- mb.run(alognet, fun=tloglin(pool.rate="rel", method.rate="common"), covar="AR1",
                         rho="dunif(0,1)", n.iter=1500, pd="pv")

    resdev <- mb.run(alognet, fun=tpoly(degree=1), parameters.to.save = "resdev", n.iter=1000, pd="pv")


    model.list <- list(loglin, emax, bs, class, ls, loglin.ar1)
    treats.list <- list(c(1,5,8,15), c("alog_50", "alog_25"), NULL, c(3,5,7), NULL, NULL)
    ref.resp.list <- list(painnet$data.ab[painnet$data.ab$treatment==1,],
                          alognet$data.ab[alognet$data.ab$treatment==2,],
                          list(beta.1=~rnorm(n,-0.1,0.1), beta.3=~rnorm(n, 0.2, 0.01), beta.4=~rnorm(n,0,0)),
                          list(beta.1=-1, beta.2=0.1),
                          NULL,
                          NULL
                          )
    times.list <- list(c(0:10), c(1,10:20), seq(0, max(bs$model.arg$jagsdata$time, na.rm=TRUE), length.out=20),
                       c(2,4,6,7,10), seq(0, max(ls$model.arg$jagsdata$time, na.rm=TRUE), length.out=20),
                       c(0:20)
                       )
    E0.list <- list(7, ~rnorm(n, 7,2), 0, ~rnorm(n,5,5), 10, 0)
    synth.list <- rep(c("common", "random"),3)

    for (i in 1:6) {
      print(paste0("modellist: ", i))
      mbnma <- model.list[[i]]
      E0 <- E0.list[[i]]
      ref.resp <- ref.resp.list[[i]]
      treats <- treats.list[[i]]
      times <- times.list[[i]]
      synth <- synth.list[[i]]

      # Tests using ref.resp
      pred <- suppressWarnings(predict(mbnma, times=times,
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

      testthat::expect_identical(names(pred), c("summary", "pred.mat", "network", "times", "link"))
      testthat::expect_equal(nrow(pred$pred.mat[[1]]), mbnma$BUGSoutput$n.sims)
      testthat::expect_equal(nrow(pred$summary[[1]]), length(times))
      testthat::expect_equal(identical(pred$summary[[1]]$time, times), TRUE)
    }

    # Tests of class models
    expect_error(predict(class, level="class"), "all relative effects must be modelled with class")

    class2 <- mb.run(goutnet, fun=tpoly(degree=2, method.1="common", pool.1="rel",
                                        method.2="random", pool.2="rel"),
                    class.effect = list(beta.1="random", beta.2="common"), pd="pv"
    )

    pred <- predict(class2, level="class")
    testthat::expect_equal(names(pred$pred.mat), class2$network$classes)

    pred <- predict(class2, level="class", treats = c("Allo", "Arha"))
    testthat::expect_equal(names(pred$pred.mat), c("Allo", "Arha"))

    pred <- predict(class2, level="class", treats = c(1,3,5))
    testthat::expect_equal(names(pred$pred.mat), class2$network$classes[c(1,3,5)])



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



  testthat::test_that("ref.synth functions correctly", {
    ref.resp <- osteopain[osteopain$treatname=="Placebo_0",]

    testthat::expect_error(ref.synth(ref.resp, emax, synth="random",
                                     n.burnin=100, n.iter=200), NA)

  })




  testthat::test_that("get.model.vals functions correctly", {

    vals <- get.model.vals(loglin)
    testthat::expect_equal(names(vals), c("alpha", "d.1", "timecourse", "time.params"))

    vals <- get.model.vals(emax)
    testthat::expect_equal(names(vals), c("alpha", "d.1", "beta.2", "timecourse", "time.params"))

    # Can add more here if problems
  })
})



