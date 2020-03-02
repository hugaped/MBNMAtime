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
  testthat::expect_equal(sort(comp[1,1:2]), comp[1,1:2])

  network <- mb.network(osteopain, reference = "Ce_200")
  comp <- mb.nodesplit.comparisons(network)
  testthat::expect_equal(nrow(comp), 6)
  testthat::expect_equal(names(comp), c("t1", "t2", "path"))
  checkmate::expect_class(comp$t1, "numeric")
  checkmate::expect_class(comp$t2, "numeric")
  checkmate::expect_class(comp$path, "factor")
  testthat::expect_equal(sort(comp[1,1:2]), comp[1,1:2])
})





testthat::test_that("mb.nodesplit is working", {

  # Emax time-course
  network <- mb.network(osteopain, reference = "Pl_0")
  comp <- mb.nodesplit.comparisons(network)
  nodesplit <- mb.nodesplit(network, comparisons=comp,
                              nodesplit.parameters="all",
                              fun="emax", user.fun=NULL,
                              alpha="study",
                              beta.1=list(pool="rel", method="common"),
                              beta.2=list(pool="rel", method="common"),
                              beta.3=NULL, beta.4=NULL,
                              positive.scale=TRUE, intercept=TRUE,
                              class.effect=list(),
                              parallel=TRUE,
                              n.iter=200, n.burnin=100, n.thin=1, n.chain=2)

  testthat::expect_equal(nrow(comp), length(nodesplit))
  testthat::expect_equal(names(nodesplit)[1], paste(c("d", comp[1,1:2]), collapse="."))
  checkmate::expect_list(nodesplit[[2]], len=2)
  checkmate::expect_list(nodesplit[[1]][[1]], len=9)
  checkmate::expect_list(nodesplit[[2]][[2]], len=9)
  checkmate::expect_character(nodesplit[[1]][[2]]$parameter, len=1)
  testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
  checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
  testthat::expect_equal(names(nodesplit[[2]][[2]]$quantiles), c("difference", "direct", "indirect", "nma"))
  checkmate::expect_class(nodesplit[[1]][[2]]$forest.plot, c("gg", "ggplot"))
  checkmate::expect_class(nodesplit[[2]][[1]]$density.plot, c("gg", "ggplot"))
  checkmate::expect_class(nodesplit[[1]][[1]]$direct, c("rjags", "mbnma"))
  checkmate::expect_class(nodesplit[[2]][[2]]$indirect, c("rjags", "mbnma"))


  network <- mb.network(osteopain, reference = "Ce_200")
  comp <- mb.nodesplit.comparisons(network)
  nodesplit <- mb.nodesplit(network, comparisons=comp,
                               nodesplit.parameters="all",
                               fun="emax", user.fun=NULL,
                               alpha="study",
                               beta.1=list(pool="rel", method="common"),
                               beta.2=list(pool="const", method="common"),
                               beta.3=NULL, beta.4=NULL,
                               positive.scale=TRUE, intercept=TRUE,
                               class.effect=list(),
                               parallel=TRUE,
                               n.iter=200, n.burnin=100, n.thin=1, n.chain=2)

  testthat::expect_equal(nrow(comp), length(nodesplit))
  testthat::expect_equal(names(nodesplit)[5], paste(c("d", comp[5,1:2]), collapse="."))
  checkmate::expect_list(nodesplit[[2]], len=1) # length= n parameters
  checkmate::expect_list(nodesplit[[3]][[1]], len=9)
  checkmate::expect_list(nodesplit[[4]][[1]], len=9)
  checkmate::expect_character(nodesplit[[5]][[1]]$parameter, len=1)
  testthat::expect_equal(names(nodesplit[[6]][[1]]$`overlap matrix`), c("direct", "indirect"))
  checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
  testthat::expect_equal(names(nodesplit[[2]][[1]]$quantiles), c("difference", "direct", "indirect", "nma"))
  checkmate::expect_class(nodesplit[[3]][[1]]$forest.plot, c("gg", "ggplot"))
  checkmate::expect_class(nodesplit[[4]][[1]]$density.plot, c("gg", "ggplot"))
  checkmate::expect_class(nodesplit[[5]][[1]]$direct, c("rjags", "mbnma"))
  checkmate::expect_class(nodesplit[[1]][[1]]$indirect, c("rjags", "mbnma"))




  # Piecewise linear time-course
  network <- mb.network(osteopain, reference = "Pl_0")
  comp <- mb.nodesplit.comparisons(network)
  nodesplit <- mb.nodesplit(network, comparisons=comp,
                               nodesplit.parameters="all",
                               fun="piecelinear", user.fun=NULL,
                               alpha="study",
                               beta.1=list(pool="rel", method="common"),
                               beta.2=list(pool="rel", method="random"),
                               beta.3=list(pool="arm", method="common"),
                               beta.4=NULL,
                               positive.scale=TRUE, intercept=TRUE,
                               class.effect=list(),
                               parallel=TRUE,
                               n.iter=200, n.burnin=100, n.thin=1, n.chain=2)

  testthat::expect_equal(nrow(comp), length(nodesplit))
  testthat::expect_equal(names(nodesplit)[1], paste(c("d", comp[1,1:2]), collapse="."))
  checkmate::expect_list(nodesplit[[2]], len=2)
  checkmate::expect_list(nodesplit[[1]][[1]], len=9)
  checkmate::expect_list(nodesplit[[2]][[2]], len=9)
  checkmate::expect_character(nodesplit[[1]][[2]]$parameter, len=1)
  testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
  checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
  testthat::expect_equal(names(nodesplit[[2]][[2]]$quantiles), c("difference", "direct", "indirect", "nma"))
  checkmate::expect_class(nodesplit[[1]][[2]]$forest.plot, c("gg", "ggplot"))
  checkmate::expect_class(nodesplit[[2]][[1]]$density.plot, c("gg", "ggplot"))
  checkmate::expect_class(nodesplit[[1]][[1]]$direct, c("rjags", "mbnma"))
  checkmate::expect_class(nodesplit[[2]][[2]]$indirect, c("rjags", "mbnma"))

  testthat::expect_equal("beta.3" %in% nodesplit[[1]][[2]]$direct[["parameters.to.save"]], TRUE)
  testthat::expect_equal("beta.3" %in% nodesplit[[1]][[2]]$indirect[["parameters.to.save"]], TRUE)



  # Gout data (emax)
  network <- mb.network(goutSUA_CFB, reference = "Plac")
  comp <- mb.nodesplit.comparisons(network)
  nodesplit <- mb.nodesplit(network, comparisons=comp,
                               nodesplit.parameters="all",
                               fun="emax", user.fun=NULL,
                               alpha="study", beta.1=list(pool="rel", method="common"),
                               beta.2=list(pool="rel", method="common"),
                               beta.3=NULL, beta.4=NULL,
                               positive.scale=TRUE, intercept=TRUE,
                               class.effect=list(),
                               parallel=TRUE,
                               n.iter=200, n.burnin=100, n.thin=1, n.chain=2)

  testthat::expect_equal(nrow(comp), length(nodesplit))
  testthat::expect_equal(names(nodesplit)[1], paste(c("d", comp[1,1:2]), collapse="."))
  checkmate::expect_list(nodesplit[[2]], len=2)
  checkmate::expect_list(nodesplit[[1]][[1]], len=9)
  checkmate::expect_list(nodesplit[[2]][[2]], len=9)
  checkmate::expect_character(nodesplit[[1]][[2]]$parameter, len=1)
  testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
  checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
  testthat::expect_equal(names(nodesplit[[2]][[2]]$quantiles), c("difference", "direct", "indirect", "nma"))
  checkmate::expect_class(nodesplit[[1]][[2]]$forest.plot, c("gg", "ggplot"))
  checkmate::expect_class(nodesplit[[2]][[1]]$density.plot, c("gg", "ggplot"))
  checkmate::expect_class(nodesplit[[1]][[1]]$direct, c("rjags", "mbnma"))
  checkmate::expect_class(nodesplit[[2]][[2]]$indirect, c("rjags", "mbnma"))




  # Alogliptin dataset
  network <- mb.network(alog_pcfb)
  comp <- mb.nodesplit.comparisons(network)
  nodesplit <- mb.nodesplit(network, comparisons=comp,
                            nodesplit.parameters="all",
                            fun="emax", user.fun=NULL,
                            alpha="study", beta.1=list(pool="rel", method="common"),
                            beta.2=list(pool="arm", method="common"),
                            beta.3=NULL, beta.4=NULL,
                            positive.scale=TRUE, intercept=TRUE,
                            class.effect=list(),
                            parallel=TRUE,
                            n.iter=200, n.burnin=100, n.thin=1, n.chain=2)

  testthat::expect_equal(nrow(comp), length(nodesplit))
  testthat::expect_equal(names(nodesplit)[1], paste(c("d", comp[1,1:2]), collapse="."))
  checkmate::expect_list(nodesplit[[2]], len=1)
  checkmate::expect_list(nodesplit[[3]][[1]], len=9)
  checkmate::expect_character(nodesplit[[4]][[1]]$parameter, len=1)
  testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
  checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
  testthat::expect_equal(names(nodesplit[[2]][[1]]$quantiles), c("difference", "direct", "indirect", "nma"))
  checkmate::expect_class(nodesplit[[3]][[1]]$forest.plot, c("gg", "ggplot"))
  checkmate::expect_class(nodesplit[[4]][[1]]$density.plot, c("gg", "ggplot"))
  checkmate::expect_class(nodesplit[[1]][[1]]$direct, c("rjags", "mbnma"))
  checkmate::expect_class(nodesplit[[2]][[1]]$indirect, c("rjags", "mbnma"))

})
