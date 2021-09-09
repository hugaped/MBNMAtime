testthat::context("Testing nodesplit.functions")

#### Pain data ####

network <- mb.network(osteopain)
last.data <- get.latest.time(network)


# Create dataset with mixed up treatment codes
testdata <- osteopain

testdata$treatment <- factor(as.character(testdata$treatment),
                             levels = unique(as.character(testdata$treatment)))

testnetwork <- mb.network(testdata)

last.test <- get.latest.time(testnetwork)


#### Gout data ####
network.gout <- mb.network(goutSUA_CFB)
last.gout <- get.latest.time(network.gout)


#### Alogliptin data ####
alognet <- mb.network(alog_pcfb)
last.alog <- get.latest.time(alognet)



################ Testing inconsistency.loops ################


testthat::test_that("test.inconsistency.loops", {
  testthat::expect_equal(nrow(inconsistency.loops(last.data)), nrow(inconsistency.loops(last.test)))
  testthat::expect_equal(nrow(inconsistency.loops(last.data)), 2)

  alogloop <- inconsistency.loops(last.alog)
  expect_equal(nrow(alogloop), 4)
  expect_identical(alogloop$path, unique(alogloop$path))

  goutloop <- inconsistency.loops(last.gout)
  expect_equal(nrow(goutloop), 6)
  expect_identical(goutloop$path, unique(goutloop$path))
})




testthat::test_that("mb.nodesplit.comparisons is working", {
  network <- mb.network(osteopain, reference = "Pl_0")
  comp <- mb.nodesplit.comparisons(network)
  testthat::expect_equal(nrow(comp), 2)
  testthat::expect_equal(names(comp), c("t1", "t2", "path"))
  checkmate::expect_class(comp$t1, "numeric")
  checkmate::expect_class(comp$t2, "numeric")
  checkmate::expect_class(comp$path, "factor")
  testthat::expect_equal(sort(as.matrix(comp[,1:2])[1,]), as.matrix(comp[,1:2])[1,])

  network <- mb.network(osteopain, reference = "Ce_200")
  comp <- mb.nodesplit.comparisons(network)
  testthat::expect_equal(nrow(comp), 6)
  testthat::expect_equal(names(comp), c("t1", "t2", "path"))
  checkmate::expect_class(comp$t1, "numeric")
  checkmate::expect_class(comp$t2, "numeric")
  checkmate::expect_class(comp$path, "factor")
  testthat::expect_equal(sort(as.matrix(comp[,1:2])[1,]), as.matrix(comp[,1:2])[1,])
})





testthat::test_that("mb.nodesplit is working", {
  skip_on_ci()
  skip_on_cran()

  # Emax time-course
  network <- mb.network(osteopain, reference = "Pl_0")
  comp <- mb.nodesplit.comparisons(network)
  nodesplit <- mb.nodesplit(network, comparisons=comp,
                              nodesplit.parameters="all",
                            fun=temax(pool.emax="rel", method.emax="common",
                                      pool.et50="rel", method.et50="common"),
                              positive.scale=TRUE, intercept=TRUE,
                              class.effect=list(),
                              n.iter=200, n.burnin=100, n.thin=1, n.chain=2)

  testthat::expect_equal(nrow(comp), length(nodesplit))
  testthat::expect_equal(names(nodesplit)[1], "Ro_25 vs Ce_200")
  checkmate::expect_list(nodesplit[[2]], len=2)
  checkmate::expect_list(nodesplit[[1]][[1]], len=6)
  checkmate::expect_list(nodesplit[[2]][[2]], len=6)
  checkmate::expect_character(nodesplit[[1]][[2]]$parameter, len=1)
  testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
  checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
  testthat::expect_equal(names(nodesplit[[2]][[2]]$quantiles), c("difference", "direct", "indirect", "nma"))

  network <- mb.network(osteopain, reference = "Ce_200")
  comp <- mb.nodesplit.comparisons(network)
  nodesplit <- mb.nodesplit(network, comparisons=comp,
                               nodesplit.parameters="all",
                            fun=temax(pool.emax="rel", method.emax="common",
                                      pool.et50="abs", method.et50="common"),
                               positive.scale=TRUE, intercept=TRUE,
                               class.effect=list(),
                               parallel=TRUE,
                               n.iter=200, n.burnin=100, n.thin=1, n.chain=2)

  testthat::expect_equal(nrow(comp), length(nodesplit))
  testthat::expect_equal(names(nodesplit)[5], "Na_1000 vs Pl_0")
  checkmate::expect_list(nodesplit[[2]], len=1) # length= n parameters
  checkmate::expect_list(nodesplit[[3]][[1]], len=6)
  checkmate::expect_list(nodesplit[[4]][[1]], len=6)
  checkmate::expect_character(nodesplit[[5]][[1]]$parameter, len=1)
  testthat::expect_equal(names(nodesplit[[6]][[1]]$`overlap matrix`), c("direct", "indirect"))
  checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
  testthat::expect_equal(names(nodesplit[[2]][[1]]$quantiles), c("difference", "direct", "indirect", "nma"))



  # Piecewise linear time-course
  network <- mb.network(osteopain, reference = "Pl_0")
  comp <- mb.nodesplit.comparisons(network)
  nodesplit <- mb.nodesplit(network, comparisons=comp,
                               nodesplit.parameters="all",
                            fun=tspline(type="ls", knots = 0.1),
                               positive.scale=TRUE, intercept=TRUE,
                               class.effect=list(),
                               parallel=TRUE,
                               n.iter=200, n.burnin=100, n.thin=1, n.chain=2)

  testthat::expect_equal(nrow(comp), length(nodesplit))
  checkmate::expect_list(nodesplit[[2]], len=2)
  checkmate::expect_list(nodesplit[[1]][[1]], len=6)
  checkmate::expect_list(nodesplit[[2]][[2]], len=6)
  checkmate::expect_character(nodesplit[[1]][[2]]$parameter, len=1)
  testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
  checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
  testthat::expect_equal(names(nodesplit[[2]][[2]]$quantiles), c("difference", "direct", "indirect", "nma"))



  # Gout data (bspline)
  network <- mb.network(goutSUA_CFB, reference = "Plac")
  comp <- mb.nodesplit.comparisons(network)
  nodesplit <- mb.nodesplit(network, comparisons=comp,
                               nodesplit.parameters="all",
                               fun=tspline(type="bs", knots=2,
                                           pool.2="abs", method.2="random"),
                               positive.scale=TRUE, intercept=TRUE,
                               class.effect=list(), omega=matrix(c(10,0,0,10), nrow=2),
                               parallel=TRUE,
                               n.iter=200, n.burnin=100, n.thin=1, n.chain=2)

  testthat::expect_equal(nrow(comp), length(nodesplit))
  testthat::expect_equal(names(nodesplit)[1], "RDEA_600 vs RDEA_400")
  checkmate::expect_list(nodesplit[[2]], len=2)
  testthat::expect_equal(names(nodesplit[[2]]), c("beta.1", "beta.3"))
  checkmate::expect_list(nodesplit[[1]][[1]], len=6)
  checkmate::expect_list(nodesplit[[2]][[2]], len=6)
  checkmate::expect_character(nodesplit[[1]][[2]]$parameter, len=1)
  testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
  checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
  testthat::expect_equal(names(nodesplit[[2]][[2]]$quantiles), c("difference", "direct", "indirect", "nma"))



  # Alogliptin dataset
  network <- mb.network(alog_pcfb)
  comp <- mb.nodesplit.comparisons(network)
  nodesplit <- mb.nodesplit(network, comparisons=comp,
                            nodesplit.parameters="all",
                            fun=tloglin(pool.rate="rel", method.rate="common"),
                            positive.scale=TRUE, intercept=TRUE,
                            class.effect=list(),
                            parallel=TRUE,
                            n.iter=200, n.burnin=100, n.thin=1, n.chain=2)

  testthat::expect_equal(nrow(comp), length(nodesplit))
  testthat::expect_equal(names(nodesplit[[2]]), c("rate"))
  checkmate::expect_list(nodesplit[[2]], len=1)
  checkmate::expect_list(nodesplit[[3]][[1]], len=6)
  checkmate::expect_character(nodesplit[[4]][[1]]$parameter, len=1)
  testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
  checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
  testthat::expect_equal(names(nodesplit[[2]][[1]]$quantiles), c("difference", "direct", "indirect", "nma"))


  expect_error(mb.nodesplit(network, comparisons=comp, nodesplit.parameters="all", fun=tloglin(pool.rate="abs", method.rate="common")),
                            "Parameter specified for nodesplit.parameters")
})
