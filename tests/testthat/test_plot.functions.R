testthat::context("Testing plot.functions")

painnet <- mb.network(osteopain)
copdnet <- mb.network(copd)
alognet <- mb.network(alog_pcfb)

# Class data
classdata <- osteopain
classdata$class[classdata$treatment=="Pl_0"] <- 1
classdata$class[classdata$treatment!="Pl_0"] <- 2
classnetwork <- mb.network(classdata)



testthat::test_that("plot.functions with mbnma models correctly", {
  skip_on_ci()
  skip_on_cran()

  emax <- mb.run(painnet, fun=temax(pool.emax="rel", method.emax="random",
                                    pool.et50="rel", method.et50="common"),
                                    n.chain=3, n.iter=1200, n.burnin=800)

  bs <- mb.run(copdnet, fun=tspline(type="bs", knots=2,
                                    pool.1="abs", method.1="random",
                                    pool.2="rel", method.2="common",
                                    pool.3="rel", method.3="common"
                                    ), omega=matrix(c(10,0,0,10), nrow=2),
               n.chain=3, n.iter=1200, n.burnin=800, intercept = FALSE)

  loglin <- mb.run(alognet, fun=tloglin(pool.rate="rel", method.rate="random"),
                   n.chain=3, n.iter=1200, n.burnin=800, intercept = FALSE)


  emax.class.random <- suppressWarnings(
    mb.run(classnetwork,
           fun=temax(pool.emax="rel", method.emax="common",
                     pool.et50="rel", method.et50="common"),
           positive.scale=TRUE,
           n.chain=3, n.iter=1200, n.burnin=800,
           class.effect=list("et50"="random")
           )
    )




  testthat::test_that("plot(mb.predict) functions correctly", {
    skip_on_ci()
    skip_on_cran()

    treats <- c(1,5,8,10)
    data.ab <- painnet[["data.ab"]]
    ref.estimate <- data.ab[data.ab$treatname=="Placebo_0",]
    times <- c(0:15)
    baseline <- 7

    pred.emax <- suppressWarnings(
      predict(emax, times=times, E0=baseline, treats=treats,
              ref.resp=ref.estimate)
    )

    pred.bs <- suppressWarnings(
      predict(bs, E0=0, ref.resp=copdnet$data.ab[copdnet$data.ab$treatment==1,])
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
    predict <- suppressWarnings(
      predict(object=emax, E0=7, treats=c(2:5),
              ref.resp=ref.estimate)
    )


    # Test overlay.nma
    expect_error(plot(pred.emax, overlay.nma=c(1,5), NA))

    expect_error(plot(pred.bs, method="common", overlay.nma = c(10,20), disp.obs = TRUE), NA)

    expect_error(plot(pred.loglin, overlay.nma=c(10,30)), NA)

    expect_error(plot(pred.loglin, overlay.nma=c(1,1.5)), "Network reference treatment")

    expect_error(plot(pred.bs, overlay.ref=FALSE, overlay.nma=c(10,20)), "must be TRUE")
    expect_error(plot(pred.emax, overlay.nma=c(0,20)), "Assertion on")
    expect_error(plot(pred.emax, overlay.nma=c(1,200)), "Assertion on")


    testthat::expect_error(plot(predict, disp.obs = FALSE, overlay.ref=TRUE,
                                max.col.scale = NULL), "must be included")
  })






  testthat::test_that("plot(mbnma) functions correctly", {
    testthat::expect_silent(plot(emax, params=NULL))
    testthat::expect_silent(plot(bs, params=NULL))
    testthat::expect_silent(plot(emax.class.random, params=NULL))
    testthat::expect_silent(plot(loglin, params=NULL))

    testthat::expect_silent(plot(bs, params=c("beta.2", "beta.3")))

    forest <- plot(emax, params=NULL)
    testthat::expect_equal(length(unique(forest$data$timeparam)), 2)

    forest <- plot(bs, params="beta.2")
    testthat::expect_equal(length(unique(forest$data$timeparam)), 1)

    testthat::expect_error(plot(bs, params="beta.1"), "does not vary")

    forest <- plot(loglin, params=NULL)
    testthat::expect_equal(length(unique(forest$data$timeparam)), 1)

    forest <- plot(emax.class.random, params=NULL)

    testthat::expect_error(plot(bs, params=c("badgers")), "Must contain elements of set")
    testthat::expect_error(plot(emax, params=1))
  })




  test_that("devplot functions correctly", {
    skip_on_ci()
    skip_on_cran()

    g <- devplot(emax, dev.type="resdev", n.iter=100)
    expect_identical(class(g$graph), c("gg", "ggplot"))
    expect_identical(names(g), c("graph", "dev.data"))

    g <- devplot(bs, dev.type="dev", n.iter=100)
    expect_equal(any(g$dev.data$mean<0), TRUE)

    g <- devplot(loglin, dev.type="dev", n.iter=100)
    expect_equal(any(g$dev.data$mean<0), TRUE)

    g <- devplot(bs, plot.type="box", facet=FALSE, xaxis="fup", n.iter=100)
    expect_equal(!("facet" %in% names(g$dev.data)), TRUE)
    expect_equal(is.integer(g$dev.data$fup), TRUE)

    g <- devplot(emax.class.random)
    expect_identical(class(g$graph), c("gg", "ggplot"))
    expect_identical(names(g), c("graph", "dev.data"))

    expect_error(devplot(emax, dev.type="theta"))
    expect_error(devplot(emax.class.random, facet=NULL))
    expect_error(devplot(emax1, plot.type="badger"))
    expect_error(devplot(emax2, xaxis="lemming"))
  })



  test_that("fitplot functions correctly", {
    skip_on_ci()
    skip_on_cran()

    g <- fitplot(emax, n.iter=100)
    expect_identical(class(g), c("gg", "ggplot"))

    g <- fitplot(bs, n.iter=100)
    expect_identical(class(g), c("gg", "ggplot"))

    expect_error(fitplot(loglin, disp.obs=FALSE, n.iter=100), NA)
    expect_error(fitplot(emax, disp.obs=FALSE, treat.labs=c(1:29), n.iter=100), NA)
    expect_error(fitplot(emax, disp.obs=FALSE, treat.labs=c(1:10), n.iter=100), "treat.labs must be the same length")

  })





  test_that("plot.mb.rank functions correctly", {
    skip_on_ci()
    skip_on_cran()

    ranks <- rank(emax, params=c("auc", "emax"), direction=-1, n.iter=10)
    g <- plot(ranks)
    expect_equal(length(g), 2)
    expect_identical(class(g[[1]]), c("gg", "ggplot"))
    expect_silent(plot(ranks))

    ranks <- rank(emax.class.random, params=c("emax"), direction=-1)
    g <- plot(ranks)
    expect_equal(length(levels(g$emax$data$treat)), 29)

    ranks <- rank(emax.class.random, params=c("ET50"), direction=-1, level="class")
    g <- plot(ranks)
    expect_equal(length(levels(g$ET50$data$treat)), 2)

    expect_silent(plot(ranks, treat.labs = c("placebo", "active")))
    expect_error(plot(ranks, treat.labs = c("placebo", "active", "extra")), "must be the same length")

    ranks <- rank(bs, params=c("auc", "d.2"), direction=-1, n.iter=10)
    expect_silent(plot(ranks))
  })
})


############################################################


test_that("plot.mb.nodesplit functions correctly", {
  skip_on_ci()
  skip_on_cran()

  network <- mb.network(osteopain, reference = "Pl_0")
  comp <- mb.nodesplit.comparisons(network)
  nodesplit <- mb.nodesplit(network, comparisons=comp,
                            nodesplit.parameters="all",
                            fun=temax(pool.emax="rel", pool.et50="rel",
                                      method.emax="common", method.et50="common"),
                            positive.scale=TRUE, intercept=TRUE,
                            class.effect=list(),
                            parallel=TRUE,
                            n.iter=200, n.burnin=100, n.thin=1, n.chain=2)

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




testthat::test_that("plot(mb.network) functions correctly", {
  testthat::expect_silent(plot(painnet, layout=igraph::in_circle(),
                               edge.scale=1, label.distance=0))

  testthat::expect_silent(plot(alognet, layout=igraph::as_star(),
                               edge.scale=1, label.distance=0))

  testthat::expect_silent(plot(copdnet, layout = igraph::with_fr(),
                               edge.scale=10, label.distance=0))

  testthat::expect_silent(plot(painnet, layout = igraph::with_fr(),
                               edge.scale=0.5, label.distance=10))

  testthat::expect_silent(plot(alognet, layout=igraph::in_circle(),
                               edge.scale=0.5, label.distance=-10))

  testthat::expect_error(plot.mb.network(painnet[["data.abd"]], layout=igraph::in_circle(),
                                         edge.scale=0.5, label.distance=-10))



  network.gout <- mb.network(goutSUA_CFB)
  testthat::expect_warning(plot(network.gout, layout=igraph::as_star(),
                                edge.scale=1, label.distance=0))
  testthat::expect_silent(plot(network.gout, layout=igraph::in_circle(),
                               level="class", remove.loops=TRUE))
  testthat::expect_warning(plot(network.gout, layout=igraph::in_circle(),
                                level="treatment"))
  testthat::expect_error(plot(network.gout, layout=igraph::in_circle(),
                              level="apple"))


})





test_that("timeplot functions correctly", {

  expect_silent(timeplot(painnet))
  expect_silent(timeplot(alognet))

  expect_silent(timeplot(painnet, plotby="rel"))

  expect_error(timeplot(painnet, level="class"), "cannot be set to class")

  goutnet <- mb.network(goutSUA_CFBcomb)
  expect_error(timeplot(goutnet, level="treatment"), NA)
  expect_error(timeplot(goutnet, level="class"), NA)
  expect_silent(timeplot(classnetwork, plotby="rel", level="class"))

})
