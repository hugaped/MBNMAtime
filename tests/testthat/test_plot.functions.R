testthat::context("Testing plot.functions")


datalist <- list(osteopain=osteopain, copd=copd, goutSUA_CFBcomb=goutSUA_CFBcomb,
                 hyalarthritis=hyalarthritis, diabetes=diabetes, alog_pcfb=alog_pcfb)

testthat::test_that("plot function tests pass correctly", {

  skip_on_appveyor()
  skip_on_ci()
  skip_on_cran()

  n.iter <- 200
  n.burnin <- 100
  n.thin <- 1
  seed <- 890421

  testthat::expect_equal(1,1) # Avoids empty test

  for (dat in seq_along(datalist)) {

    print(names(datalist)[dat])

    network <- mb.network(datalist[[dat]])

    # SUPPRESSES WARNINGS FOR VERSION 0.2.2 - REMOVE AFTER THIS AND TEST WITHOUT TO ENSURE WARNINGS IDENTIFIED
    suppressWarnings({
    emax <- mb.run(network, fun=temax(pool.emax="rel", method.emax="random",
                                      pool.et50="rel", method.et50="common"),
                   n.chain=3, n.iter=1200, n.burnin=800, jags.seed=seed)

    bs <- mb.run(network, fun=tspline(type="bs", knots=2,
                                      pool.1="abs", method.1="random",
                                      pool.2="rel", method.2="common",
                                      pool.3="rel", method.3="common"
    ), omega=matrix(c(10,0,0,10), nrow=2),
    n.chain=3, n.iter=1200, n.burnin=800, intercept = FALSE, jags.seed=seed)

    loglin <- mb.run(network, fun=tloglin(pool.rate="rel", method.rate="random"),
                     n.chain=3, n.iter=1200, n.burnin=800, intercept = FALSE, jags.seed=seed)


    if ("class" %in% names(datalist[[dat]])) {

      emax.class.random <- suppressWarnings(
        mb.run(network,
               fun=temax(pool.emax="rel", method.emax="common",
                         pool.et50="rel", method.et50="common"),
               positive.scale=TRUE,
               n.chain=3, n.iter=1200, n.burnin=800, jags.seed=seed,
               class.effect=list("et50"="random")
        )
      )
    }

    })


    testthat::test_that(paste0("plot(mb.predict) functions correctly for ", names(datalist)[dat]), {

      if (length(unique(datalist[[dat]]$treatment))>6) {
        treats <- c(1,3,4,6)
      } else {
        treats <- c(1,3)
      }


      data.ab <- network[["data.ab"]]
      ref.estimate <- network$data.ab[network$data.ab$treatment==1,]
      times <- c(0:15)
      baseline <- 7

      pred.emax <- suppressWarnings(
        predict(emax, times=times, E0=baseline, treats=treats,
                ref.resp=ref.estimate)
      )

      pred.bs <- suppressWarnings(
        predict(bs, E0=0, ref.resp=ref.estimate)
      )

      pred.loglin <- predict(loglin, E0=0,
                             ref.resp=list(rate=0))


      predicts <- list(pred.emax, pred.bs, pred.loglin)

      for (i in seq_along(predicts)) {
        print(i)

        predict <- predicts[[i]]


        #### Silent runs ####

        testthat::expect_silent(plot(predict, disp.obs = FALSE, overlay.ref=FALSE,
                                     max.col.scale = NULL))

        testthat::expect_silent(plot(predict, disp.obs = FALSE,
                                     overlay.ref=FALSE, max.col.scale = "jazz"))

        # overlay.ref
        testthat::expect_message(plot(predict, disp.obs = FALSE, overlay.ref=TRUE,
                                      max.col.scale = NULL),
                                 "Reference treatment in plots is")

        # disp.obs
        testthat::expect_silent(plot(predict, disp.obs = TRUE, overlay.ref=FALSE,
                                     max.col.scale = NULL))

        testthat::expect_silent(plot(predict, disp.obs = TRUE, overlay.ref=FALSE,
                                     max.col.scale = 100))

        testthat::expect_silent(plot(predict, disp.obs = TRUE, overlay.ref=FALSE,
                                     col = "green"))

        # both
        testthat::expect_message(plot(predict, disp.obs = TRUE, overlay.ref=TRUE))


        #### error runs ####

        testthat::expect_error(plot(predict, disp.obs = TRUE, overlay.ref=FALSE,
                                    max.col.scale = 0))

        testthat::expect_error(plot(predict, disp.obs = TRUE, overlay.ref=FALSE,
                                    max.col.scale = "jazz"))

      }

      # Change treats to not include 1
      if (length(unique(datalist[[dat]]$treatment))>6) {
        treats <- c(2,3,4,6)
      } else {
        treats <- c(2,3)
      }

      pred <- suppressWarnings(
        predict(object=emax, E0=7, treats=treats,
                ref.resp=ref.estimate)
      )


      # Test overlay.nma

      expect_error(plot(pred.emax, overlay.nma=quantile(network$data.ab$time)[2:3], NA))

      expect_error(plot(pred.bs, method="common", overlay.nma = quantile(network$data.ab$time)[c(1,3)],
                        disp.obs = TRUE), NA)

      expect_error(plot(pred.loglin, overlay.nma=quantile(network$data.ab$time)[c(1:3)]), NA)

      expect_warning(plot(pred.loglin, overlay.nma=c(1.099,1.1)), "No time-point in predicted values is between")

      expect_error(plot(pred.bs, overlay.ref=FALSE, overlay.nma=quantile(network$data.ab$time)[2:3]), "must be TRUE")
      expect_error(plot(pred.emax, overlay.nma=c(0,200)), NA)

      expect_error(plot(pred, overlay.nma=quantile(network$data.ab$time)[c(1:2)]), "Reference treatment")

    })




    testthat::test_that(paste0("plot(mbnma) functions correctly for ", names(datalist)[dat]), {
      testthat::expect_silent(plot(emax, params=NULL))
      testthat::expect_silent(plot(bs, params=NULL))
      testthat::expect_silent(plot(loglin, params=NULL))

      testthat::expect_silent(plot(bs, params=c("beta.2", "beta.3")))

      forest <- plot(emax, params=NULL)
      testthat::expect_equal(length(unique(forest$data$timeparam)), 2)

      forest <- plot(bs, params="beta.2")
      testthat::expect_equal(length(unique(forest$data$timeparam)), 1)

      testthat::expect_error(plot(bs, params="beta.1"), "does not vary")

      forest <- plot(loglin, params=NULL)
      testthat::expect_equal(length(unique(forest$data$timeparam)), 1)

      if ("class" %in% names(datalist[[dat]])) {
        testthat::expect_silent(plot(emax.class.random, params=NULL))

        expect_error(plot(emax.class.random, params="ET50"), NA)
        expect_error(plot(emax.class.random, params="emax"), NA)

        forest <- plot(emax.class.random, params=NULL)
        testthat::expect_equal(length(unique(forest$data$timeparam)), 3)
      }

      testthat::expect_error(plot(bs, params=c("badgers")), "Must contain elements of set")
      testthat::expect_error(plot(emax, params=1))
    })


    test_that(paste0("devplot functions correctly for ", names(datalist)[dat]), {

      # Warnings suppressed due to infinite deviance values from very poorly fitting models!
      g <- suppressWarnings(devplot(emax, dev.type="resdev", n.iter=100))
      expect_identical(class(g$graph), c("gg", "ggplot"))
      expect_identical(names(g), c("graph", "dev.data"))
      expect_equal(any(g$dev.data$mean<0), FALSE)


      if (!names(datalist)[dat] %in% "hyalarthritis") {

        # Warnings suppressed due to infinite deviance values from very poorly fitting models!
        g <- suppressWarnings(devplot(bs, dev.type="dev", n.iter=100))
        expect_equal(any(g$dev.data$mean<0), TRUE)

        g <- suppressWarnings(devplot(loglin, dev.type="dev", n.iter=100))
        expect_equal(any(g$dev.data$mean<0), TRUE)
      }


      # Warnings suppressed due to infinite deviance values from very poorly fitting models!
      g <- suppressWarnings(devplot(bs, plot.type="box", facet=FALSE, xaxis="fup", n.iter=100))
      expect_equal(!("facet" %in% names(g$dev.data)), TRUE)
      expect_equal(is.integer(g$dev.data$fup), TRUE)

      expect_error(suppressWarnings(devplot(emax, dev.type="theta")))
      expect_error(devplot(emax1, plot.type="badger"))
      expect_error(devplot(emax2, xaxis="lemming"))

      if ("class" %in% names(datalist[[dat]])) {
        g <- suppressWarnings(devplot(emax.class.random))
        expect_identical(class(g$graph), c("gg", "ggplot"))
        expect_identical(names(g), c("graph", "dev.data"))

        expect_error(devplot(emax.class.random, facet=NULL))
      }
    })


    test_that(paste0("fitplot functions correctly for ", names(datalist)[dat]), {

      # Suppresses warning: In unique.default(x) : reached elapsed time limit
      g <- suppressWarnings(fitplot(emax, n.iter=100))
      expect_identical(class(g), c("gg", "ggplot"))

      g <- suppressWarnings(fitplot(bs, n.iter=100))
      expect_identical(class(g), c("gg", "ggplot"))

      g <- suppressWarnings(fitplot(loglin, n.iter=100))
      expect_identical(class(g), c("gg", "ggplot"))

      expect_error(suppressWarnings(fitplot(loglin, disp.obs=FALSE, n.iter=100)), NA)
      expect_error(fitplot(emax, disp.obs=FALSE, treat.labs=paste0("badger", 1:length(network$treatments)), n.iter=100), NA)
      expect_error(fitplot(emax, disp.obs=FALSE, treat.labs=c(1:100), n.iter=100), "treat.labs must be the same length")

    })


    test_that(paste0("plot.mb.rank functions correctly for ", names(datalist)[dat]), {

      ranks <- rank(emax, param=c("auc"), direction=-1, n.iter=10)
      g <- plot(ranks)
      expect_s3_class(g, "ggplot")
      expect_silent(plot(ranks))

      ranks <- rank(emax, param=c("emax"), direction=-1)
      g <- plot(ranks)
      expect_equal(length(levels(g$data$treat)), length(network$treatments))

      ranks <- rank(bs, param=c("d.2"), direction=-1, n.iter=10)
      expect_silent(plot(ranks))

      if ("class" %in% names(datalist[[dat]])) {
        ranks <- rank(emax.class.random, param=c("ET50"), direction=-1, level="class")
        g <- plot(ranks)
        expect_equal(length(levels(g$data$treat)), length(network$classes))

        expect_silent(plot(ranks, treat.labs = paste0("badger", 1:length(network$classes))))
        expect_error(plot(ranks, treat.labs = c(network$classes, "extra")), "must be the same length")
      }
    })


    # Only test for datasets in which inconsistency can be assessed

    if (names(datalist)[dat] %in% c("osteopain", "goutSUA_CFBcomb", "hyalarthritis")) {

      test_that(paste0("plot.mb.nodesplit functions correctly for ", names(datalist)[dat]), {

        network <- mb.network(datalist[[dat]])
        comp <- mb.nodesplit.comparisons(network)

        # SUPPRESSES WARNINGS FOR VERSION 0.2.2 - REMOVE AFTER THIS AND TEST WITHOUT TO ENSURE WARNINGS IDENTIFIED
        suppressWarnings(
        nodesplit <- mb.nodesplit(network, comparisons=comp,
                                  nodesplit.parameters="all",
                                  fun=temax(pool.emax="rel", pool.et50="rel",
                                            method.emax="common", method.et50="common"),
                                  positive.scale=TRUE, intercept=TRUE,
                                  class.effect=list(),
                                  n.iter=200, n.burnin=100, n.thin=1, n.chain=2, jags.seed=seed)
        )

        g <- plot.nodesplit(nodesplit)
        expect_equal(length(g), 2)
        expect_equal(length(plot.nodesplit(nodesplit, params = "emax")), 2)
        expect_equal(length(plot.nodesplit(nodesplit, params = "et50")), 2)
        expect_equal(length(plot.nodesplit(nodesplit, plot.type = "forest")), 1)
        expect_equal(length(plot.nodesplit(nodesplit, plot.type = "density")), 1)

        expect_error(plot.nodesplit(nodesplit, params="badgers"), "not a time-course parameter")
        expect_error(plot.nodesplit(nodesplit, plot.type="badgers"))

        g <- plot.nodesplit(nodesplit, plot.type=NULL)
        testthat::expect_equal(length(g), 2)

      })
    }



    testthat::test_that(paste0("plot(mb.network) functions correctly for ", names(datalist)[dat]), {
      testthat::expect_silent(plot(network, layout=igraph::in_circle(),
                                   edge.scale=1, label.distance=0))

      testthat::expect_silent(plot(network, layout=igraph::as_star(),
                                   edge.scale=1, label.distance=0))

      testthat::expect_silent(plot(network, layout = igraph::with_fr(),
                                   edge.scale=10, label.distance=0))

      testthat::expect_silent(plot(network, layout = igraph::with_fr(),
                                   edge.scale=0.5, label.distance=10))

      testthat::expect_silent(plot(network, layout=igraph::in_circle(),
                                   edge.scale=0.5, label.distance=-10))


      if ("class" %in% names(datalist[[dat]])) {
        testthat::expect_silent(plot(network, layout=igraph::as_star(),
                                      edge.scale=1, label.distance=0))
        testthat::expect_silent(plot(network, layout=igraph::in_circle(),
                                     level="class", remove.loops=TRUE))
        testthat::expect_silent(plot(network, layout=igraph::in_circle(),
                                      level="treatment"))
        testthat::expect_error(plot(network, layout=igraph::in_circle(),
                                    level="apple"))
      }

    })



    test_that(paste0("timeplot functions correctly for ", names(datalist)[dat]), {

      expect_silent(timeplot(network))

      expect_silent(timeplot(network, plotby="rel"))


      if ("class" %in% names(datalist[[dat]])) {

        expect_error(timeplot(network, level="treatment"), NA)
        expect_error(timeplot(network, level="class"), NA)
        expect_silent(timeplot(network, plotby="rel", level="class"))

      } else {
        expect_error(timeplot(network, level="class"), "cannot be set to class")
      }

    })



    test_that(paste0("binplot functions correctly for ", names(datalist)[dat]), {

      network <- suppressWarnings(mb.network(datalist[[dat]]))

      expect_error(binplot(network), NA)

      if (names(datalist)[dat]=="diabetes") {
        expect_error(binplot(network, overlay.nma = c(0,5)), "No NMA can be performed between")
      } else {
        expect_message(binplot(network, overlay.nma = c(0,5,5.001)), "not possible between")
      }


      expect_error(binplot(network, overlay.nma=10), "Must have length >= 2")

    })

  }


  testthat::test_that("alpha.scale functions correctly", {
    n.cut <- 10
    cols <- MBNMAtime:::alpha.scale(n.cut, col="blue")
    testthat::expect_equal(length(cols), n.cut+1)

    testthat::expect_silent(MBNMAtime:::alpha.scale(n.cut, col="green"))
    testthat::expect_silent(MBNMAtime:::alpha.scale(n.cut, col="red"))
    testthat::expect_silent(MBNMAtime:::alpha.scale(n.cut, col=c(10,10,50)))

    testthat::expect_error(MBNMAtime:::alpha.scale(n.cut, col="cheesecake"))
    testthat::expect_error(MBNMAtime:::alpha.scale(n.cut, col=c(10,10)))
    testthat::expect_error(MBNMAtime:::alpha.scale(n.cut, col=c(-10,10,50)))
    testthat::expect_error(MBNMAtime:::alpha.scale(n.cut, col=c(10,10,500)))

  })
})


