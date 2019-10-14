testthat::context("Testing plot.functions")
network <- mb.network(osteopain)

# emax1 <- mb.emax(network,
#                     emax=list(pool="rel", method="random"),
#                     et50=list(pool="rel", method="common"),
#                     positive.scale=TRUE,
#                     n.chain=3, n.iter=1200, n.burnin=800)
#
# emax2 <- mb.emax(network,
#                     emax=list(pool="rel", method="common"),
#                     et50=list(pool="const", method="common"),
#                     positive.scale=TRUE,
#                     n.chain=3, n.iter=1200, n.burnin=800)
#
# emax.hill <- mb.emax.hill(network,
#                              emax=list(pool="arm", method="common"),
#                              et50=list(pool="arm", method="random"),
#                              hill=list(pool="const", method=-0.5),
#                              positive.scale=TRUE, n.chain=3, n.iter=1200, n.burnin=800)



# Class data
classdata <- osteopain
classdata$class[classdata$treatment=="Pl_0"] <- 1
classdata$class[classdata$treatment!="Pl_0"] <- 2
classnetwork <- mb.network(classdata)

# emax.class.random <- mb.emax(classnetwork,
#                                 emax=list(pool="rel", method="common"),
#                                 et50=list(pool="rel", method="common"),
#                                 positive.scale=TRUE,
#                                 n.chain=3, n.iter=1200, n.burnin=800,
#                                 class.effect=list("et50"="random"))
#
# emax.class.2 <- mb.emax(classnetwork,
#                                 emax=list(pool="arm", method="common"),
#                                 et50=list(pool="arm", method="random"),
#                                 positive.scale=TRUE,
#                                 n.chain=3, n.iter=1200, n.burnin=800,
#                                 class.effect=list("emax"="common", "et50"="random"))
#
#
# treats <- c(1:10)
# data.ab <- network[["data.ab"]]
# ref.estimate <- data.ab[data.ab$treatname=="Placebo_0",]
# times <- c(0:15)
# baseline <- 7
#
# predict1 <- predict(emax1, times=times, baseline=baseline, treats=treats,
#                          ref.estimate=ref.estimate)
# predict2 <- predict(emax2, times=times, baseline=baseline, treats=treats,
#                           ref.estimate=ref.estimate)
# predict3 <- predict(emax.hill, times=times, baseline=baseline, treats=treats,
#                           ref.estimate=ref.estimate)
#
# predicts <- list(predict1, predict2, predict3)
#
# testthat::test_that("plot(mb.predict) functions correctly", {
#
#   for (i in seq_along(predicts)) {
#     print(i)
#
#     predict <- predicts[[i]]
#
#
#     #### Silent runs ####
#
#     testthat::expect_silent(plot(predict, disp.obs = FALSE, overlay.ref=FALSE,
#                    max.col.scale = NULL))
#
#     testthat::expect_silent(plot(predict, disp.obs = FALSE,
#                                  overlay.ref=FALSE, max.col.scale = "jazz"))
#
#     # overlay.ref
#     testthat::expect_message(plot(predict, disp.obs = FALSE, overlay.ref=TRUE,
#                                   max.col.scale = NULL),
#                    "Reference treatment in plots is")
#
#     # disp.obs
#     testthat::expect_silent(plot(predict, disp.obs = TRUE, overlay.ref=FALSE,
#                                  max.col.scale = NULL))
#
#     testthat::expect_silent(plot(predict, disp.obs = TRUE, overlay.ref=FALSE,
#                                  max.col.scale = 100))
#
#     testthat::expect_silent(plot(predict, disp.obs = TRUE, overlay.ref=FALSE,
#                                  col = "green"))
#
#     # both
#     testthat::expect_message(plot(predict, disp.obs = TRUE, overlay.ref=TRUE))
#
#
#     #### error runs ####
#
#     testthat::expect_error(plot(predict, disp.obs = TRUE, overlay.ref=FALSE,
#                                  max.col.scale = 0))
#
#     testthat::expect_error(plot(predict, disp.obs = TRUE, overlay.ref=FALSE,
#                                 max.col.scale = "jazz"))
#
#     # Change treats to not include 1
#     predict <- predict(object=emax2, times=times, baseline=baseline, treats=c(2:5),
#                              ref.estimate=ref.estimate)
#
#     testthat::expect_error(plot(predict, disp.obs = FALSE, overlay.ref=TRUE,
#                                  max.col.scale = NULL))
#     }
# })


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
  testthat::expect_silent(plot(network, layout_in_circle = TRUE,
                 edge.scale=1, label.distance=0))

  testthat::expect_silent(plot(network, layout_in_circle = FALSE,
                               edge.scale=1, label.distance=0))

  testthat::expect_silent(plot(network, layout_in_circle = FALSE,
                               edge.scale=10, label.distance=0))

  testthat::expect_silent(plot(network, layout_in_circle = FALSE,
                               edge.scale=0.5, label.distance=10))

  testthat::expect_silent(plot(network, layout_in_circle = TRUE,
                               edge.scale=0.5, label.distance=-10))

  testthat::expect_error(plot.mb.network(network[["data"]], layout_in_circle = TRUE,
                               edge.scale=0.5, label.distance=-10))

  # network2 <- network
  # network2[["data.ab"]] <- network2[["data.ab"]][-c(1:10),]
  #
  # testthat::expect_warning(plot(network2, layout_in_circle = FALSE,
  #                              edge.scale=1, label.distance=0))

  network.gout <- mb.network(goutSUA_CFB)
  testthat::expect_warning(plot(network.gout, layout_in_circle = FALSE,
                               edge.scale=1, label.distance=0))
  testthat::expect_silent(plot(network.gout, layout_in_circle = TRUE,
                               level="class", remove.loops=TRUE))
  testthat::expect_warning(plot(network.gout, layout_in_circle = TRUE,
                              level="treatment"))
  testthat::expect_error(plot(network.gout, layout_in_circle = TRUE,
                               level="apple"))


})




# testthat::test_that("plot(mbnma) functions correctly", {
#   testthat::expect_silent(plot(emax1, params=NULL))
#   testthat::expect_silent(plot(emax2, params=NULL))
#   testthat::expect_silent(plot(emax.class.random, params=NULL))
#   testthat::expect_silent(plot(emax.class.2, params=NULL))
#
#   forest <- plot(emax1, params=NULL)
#   testthat::expect_equal(length(unique(forest$data$timeparam)), 2)
#
#   forest <- plot(emax1, params="d.et50")
#   testthat::expect_equal(length(unique(forest$data$timeparam)), 1)
#
#   forest <- plot(emax2, params=NULL)
#   testthat::expect_equal(length(unique(forest$data$timeparam)), 1)
#
#   forest <- plot(emax.class.random, params=NULL)
#
#   testthat::expect_error(plot(emax2, params=c("d.emax", "d.et50")))
#   testthat::expect_error(plot(emax2, params="badgers"))
#   testthat::expect_error(plot(emax2, params="beta.et50"))
#   testthat::expect_error(plot(emax2, params=1))
#   testthat::expect_error(plot(network, params=NULL), NA)
# })




# test_that("devplot functions correctly", {
#   g <- devplot(emax1, dev.type="resdev", n.iter=100)
#   expect_identical(class(g$graph), c("gg", "ggplot"))
#   expect_identical(names(g), c("graph", "dev.data"))
#
#   g <- devplot(emax1, dev.type="dev", n.iter=100)
#   expect_equal(any(g$dev.data$mean<0), TRUE)
#
#   g <- devplot(emax2, plot.type="box", facet=FALSE, xaxis="fup", n.iter=100)
#   expect_equal(!("facet" %in% names(g$dev.data)), TRUE)
#   expect_equal(is.integer(g$dev.data$fup), TRUE)
#
#   expect_error(devplot(emax2, dev.type="theta"))
#   expect_error(devplot(emax2, facet=NULL))
#   expect_error(devplot(emax1, plot.type="badger"))
#   expect_error(devplot(emax2, xaxis="lemming"))
# })
#
#
#
# test_that("fitplot functions correctly", {
#   g <- fitplot(emax1, n.iter=100)
#   expect_identical(class(g), c("gg", "ggplot"))
#
#   expect_error(fitplot(emax.hill, disp.obs=FALSE, n.iter=100), NA)
#   expect_error(fitplot(emax.hill, disp.obs=FALSE, treat.labs=c(1:29), n.iter=100), NA)
#   expect_error(fitplot(emax.hill, disp.obs=FALSE, treat.labs=c(1:10), n.iter=100), "treat.labs must be the same length")
#
# })





# test_that("plot.mb.rank functions correctly", {
#
#   ranks <- rank(emax1, params=c("auc", "d.emax"), direction=-1, n.iter=10)
#   g <- plot(ranks)
#   expect_equal(length(g), 2)
#   expect_identical(class(g[[1]]), c("gg", "ggplot"))
#   expect_silent(plot(ranks))
#
#   ranks <- rank(emax.class.random, params=c("d.emax"), direction=-1)
#   g <- plot(ranks)
#   expect_equal(length(levels(g$d.emax$data$treat)), 29)
#
#   ranks <- rank(emax.class.random, params=c("D.et50"), direction=-1, level="class")
#   g <- plot(ranks)
#   expect_equal(length(levels(g$D.et50$data$treat)), 2)
#
#   expect_silent(plot(ranks, treat.labs = c("placebo", "active")))
#   expect_error(plot(ranks, treat.labs = c("placebo", "active", "extra")), "must be the same length")
#
#   ranks <- rank(emax.hill, params=c("auc", "beta.et50"), direction=-1, n.iter=10)
#   expect_silent(plot(ranks))
# })
#
#
#
# test_that("timeplot functions correctly", {
#
#   alognet <- mb.network(alog_pcfb)
#
#   expect_silent(timeplot(network))
#   expect_output(timeplot(alognet), "Absence of observations")
#
#   expect_error(timeplot(network, level="class"), "cannot be set to class")
#
#   expect_silent(timeplot(classnetwork, level="treatment"))
#   expect_silent(timeplot(classnetwork, level="class"))
#
# })
#
#
#
# test_that("plot.mb.nodesplit functions correctly", {
#
#   network <- mb.network(osteopain, reference = "Pl_0")
#   comp <- mb.nodesplit.comparisons(network)
#   nodesplit <- mb.nodesplit(network, comparisons=comp,
#                             nodesplit.parameters="all",
#                             fun="emax", user.fun=NULL,
#                             alpha="study",
#                             beta.1=list(pool="rel", method="common"),
#                             beta.2=list(pool="rel", method="common"),
#                             beta.3=NULL, beta.4=NULL,
#                             positive.scale=TRUE, intercept=TRUE,
#                             class.effect=list(),
#                             parallel=TRUE,
#                             n.iter=200, n.burnin=100, n.thin=1, n.chain=2)
#
#   g <- plot(nodesplit)
#   expect_equal(length(g), 4)
#   expect_equal(length(plot(nodesplit, params = "beta.1")), 2)
#   expect_equal(length(plot(nodesplit, params = "beta.2")), 2)
#   expect_equal(length(plot(nodesplit, plot.type = "forest")), 2)
#   expect_equal(length(plot(nodesplit, plot.type = "density")), 2)
#
#   expect_error(plot(nodesplit, params="badgers"))
#   expect_error(plot(nodesplit, plot.type="badgers"))
#
#   g <- plot(nodesplit, plot.type=NULL)
#   testthat::expect_equal(length(g), 4)
#
# })
