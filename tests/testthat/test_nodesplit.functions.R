testthat::context("Testing nodesplit.functions")

datalist <- list(osteopain=osteopain, goutSUA_CFBcomb=goutSUA_CFBcomb,
                 hyalarthritis=hyalarthritis)

testthat::test_that("Node-split tests pass correctly", {

  skip_on_appveyor()
  skip_on_ci()
  skip_on_cran()

  n.iter <- 200
  n.burnin <- 100
  n.thin <- 1

  seed <- 890421

  testthat::expect_equal(1,1) # Avoids empty tests

  for (i in seq_along(datalist)) {

    print(names(datalist)[i])
    network <- mb.network(datalist[[i]])
    last.data <- get.latest.time(network)$data.ab

    # Create dataset with mixed up treatment codes
    testdata <- datalist[[i]]
    testdata$treatment <- factor(as.character(testdata$treatment),
                                 levels = unique(as.character(testdata$treatment)))

    testnetwork <- mb.network(testdata)
    last.test <- get.latest.time(testnetwork)$data.ab

    testthat::test_that(paste0(names(datalist)[i], ": test.inconsistency.loops"), {

      if (names(datalist)[i]=="osteopain") {
        testthat::expect_equal(nrow(inconsistency.loops(last.data)), nrow(inconsistency.loops(last.test)))
        testthat::expect_equal(nrow(inconsistency.loops(last.data)), 2)

      } else if (names(datalist)[i]=="alog_pcfb") {

        loop <- inconsistency.loops(last.data)
        expect_equal(nrow(loop), 4)
        expect_identical(loop$path, unique(loop$path))

      } else if (names(datalist)[i]=="goutSUA_CFBcomb") {

        loop <- inconsistency.loops(last.data)
        expect_equal(nrow(loop), 7)
        expect_identical(loop$path, unique(loop$path))
      }

      testthat::expect_error(inconsistency.loops(last.data), NA)

    })



    testthat::test_that(paste0(names(datalist)[i], ": mb.nodesplit.comparisons is working"), {

      # Expect error
      compnet <- mb.network(alog_pcfb, reference =
                              ifelse(is.numeric(alog_pcfb$treatment),
                                     as.numeric(mb.network(alog_pcfb)$treatments[3]),
                                     mb.network(alog_pcfb)$treatments[3]))
      expect_error(mb.nodesplit.comparisons(compnet), "arising from independent sources")


      # Tests should pass
      comp <- mb.nodesplit.comparisons(network)

      if (names(datalist)[i]=="osteopain") {
        testthat::expect_equal(nrow(comp), 2)
      }

      testthat::expect_equal(names(comp), c("t1", "t2", "path"))
      checkmate::expect_class(comp$t1, "numeric")
      checkmate::expect_class(comp$t2, "numeric")
      checkmate::expect_class(comp$path, "factor")
      testthat::expect_equal(sort(as.matrix(comp[,1:2])[1,]), as.matrix(comp[,1:2])[1,])

      compnet <- mb.network(datalist[[i]], reference =
                              ifelse(is.numeric(datalist[[i]]$treatment),
                                                as.numeric(network$treatments[3]),
                                                network$treatments[3]))

      comp <- mb.nodesplit.comparisons(compnet)

      if (names(datalist)[i]=="osteopain") {
        testthat::expect_equal(nrow(comp), 6)
      }

      testthat::expect_equal(names(comp), c("t1", "t2", "path"))
      checkmate::expect_class(comp$t1, "numeric")
      checkmate::expect_class(comp$t2, "numeric")
      checkmate::expect_class(comp$path, "factor")
      testthat::expect_equal(sort(as.matrix(comp[,1:2])[1,]), as.matrix(comp[,1:2])[1,])

    })





    ###### mb.nodesplit ######

    testthat::test_that(paste("mb.nodesplit is working for:", names(datalist)[i]), {

      copdnet <- mb.network(copd)
      expect_error(mb.nodesplit.comparisons(copdnet), "No closed loops of treatments")

      # Emax time-course
      comp <- mb.nodesplit.comparisons(network)[1:2,]


      # REMOVE SUPPRESSWARNINGS FROM VERSION 0.2.3 ONWARNS
      suppressWarnings(
      nodesplit <- mb.nodesplit(network, comparisons=comp,
                                nodesplit.parameters="all",
                                fun=temax(pool.emax="rel", method.emax="common",
                                          pool.et50="rel", method.et50="common"),
                                positive.scale=TRUE, intercept=TRUE,
                                class.effect=list(),
                                n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, jags.seed=seed)
      )


      testthat::expect_equal(nrow(comp), length(nodesplit))
      testthat::expect_equal(any(sapply(network$treatments, function(x) {grepl(x, names(nodesplit))})),
                             TRUE)
      checkmate::expect_list(nodesplit[[2]], len=2)
      checkmate::expect_list(nodesplit[[1]][[1]], len=6)
      checkmate::expect_list(nodesplit[[2]][[2]], len=6)
      checkmate::expect_character(nodesplit[[1]][[2]]$parameter, len=1)
      testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
      checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
      testthat::expect_equal(names(nodesplit[[2]][[2]]$quantiles), c("difference", "direct", "indirect", "nma"))

      net2 <- mb.network(datalist[[i]], reference=network$treatments[3])
      comp <- mb.nodesplit.comparisons(net2)[1:2,]

      # REMOVE SUPPRESSWARNINGS FROM VERSION 0.2.3 ONWARNS
      #suppressWarnings(
      nodesplit <- mb.nodesplit(net2, comparisons=comp,
                                nodesplit.parameters="all",
                                fun=temax(pool.emax="rel", method.emax="common",
                                          pool.et50="abs", method.et50="common"),
                                positive.scale=TRUE, intercept=TRUE,
                                class.effect=list(),
                                parallel=TRUE,
                                n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, jags.seed=seed)
      #)

      testthat::expect_equal(nrow(comp), length(nodesplit))
      testthat::expect_equal(any(sapply(net2$treatments, function(x) {grepl(x, names(nodesplit))})),
                             TRUE)
      checkmate::expect_list(nodesplit[[2]], len=1) # length= n parameters
      checkmate::expect_list(nodesplit[[2]][[1]], len=6)
      checkmate::expect_character(nodesplit[[1]][[1]]$parameter, len=1)
      testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
      checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
      testthat::expect_equal(names(nodesplit[[2]][[1]]$quantiles), c("difference", "direct", "indirect", "nma"))



      # Piecewise linear time-course
      net2 <- mb.network(datalist[[i]], reference=network$treatments[2])
      comp <- mb.nodesplit.comparisons(net2)[1:2,]

      maxtime <- max(net2$data.ab$time, na.rm=TRUE)
      knots <- stats::quantile(0:maxtime, probs = c(0.1))
      names(knots) <- NULL
      nodesplit <- mb.nodesplit(net2, comparisons=comp,
                                nodesplit.parameters="all",
                                fun=tspline(type="ls", knots = knots),
                                positive.scale=TRUE, intercept=TRUE,
                                class.effect=list(),
                                parallel=TRUE,
                                n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, jags.seed=seed)

      testthat::expect_equal(nrow(comp), length(nodesplit))
      checkmate::expect_list(nodesplit[[2]], len=2)
      checkmate::expect_list(nodesplit[[1]][[1]], len=6)
      checkmate::expect_list(nodesplit[[2]][[2]], len=6)
      checkmate::expect_character(nodesplit[[1]][[2]]$parameter, len=1)
      testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
      checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
      testthat::expect_equal(names(nodesplit[[2]][[2]]$quantiles), c("difference", "direct", "indirect", "nma"))



      # Gout data (bspline)
      net2 <- mb.network(datalist[[i]])
      comp <- mb.nodesplit.comparisons(net2)[1:2,]

      nodesplit <- mb.nodesplit(net2, comparisons=comp,
                                nodesplit.parameters="all",
                                fun=tspline(type="bs", nknots=2,
                                            pool.2="abs", method.2="random"),
                                positive.scale=TRUE, intercept=TRUE, corparam=TRUE,
                                class.effect=list(),
                                parallel=TRUE,
                                n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, jags.seed=seed)

      testthat::expect_equal(nrow(comp), length(nodesplit))
      testthat::expect_equal(any(sapply(net2$treatments, function(x) {grepl(x, names(nodesplit))})),
                             TRUE)
      checkmate::expect_list(nodesplit[[2]], len=2)
      testthat::expect_equal(names(nodesplit[[2]]), c("beta.1", "beta.3"))
      checkmate::expect_list(nodesplit[[1]][[1]], len=6)
      checkmate::expect_list(nodesplit[[2]][[2]], len=6)
      checkmate::expect_character(nodesplit[[1]][[2]]$parameter, len=1)
      testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
      checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
      testthat::expect_equal(names(nodesplit[[2]][[2]]$quantiles), c("difference", "direct", "indirect", "nma"))



      # titp

      if (!names(datalist)[i] %in% c("goutSUA_CFBcomb", "hyalarthritis")) {
        if (names(datalist)[i]=="osteopain") {
          positive <- TRUE
          intercept <- TRUE
        } else {
          positive <- FALSE
          intercept <- FALSE
        }

        # REMOVE SUPPRESSWARNINGS FROM VERSION 0.2.3 ONWARNS
        #suppressWarnings(
          nodesplit <- mb.nodesplit(network, comparisons=comp,
                                    nodesplit.parameters="all",
                                    fun=titp(),
                                    positive.scale=positive, intercept=intercept,
                                    class.effect=list(),
                                    parallel=TRUE,
                                    n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, jags.seed=seed)
        #)


        testthat::expect_equal(nrow(comp), length(nodesplit))
        testthat::expect_equal(names(nodesplit[[2]]), c("emax", "rate"))
        checkmate::expect_list(nodesplit[[2]], len=2)
        checkmate::expect_list(nodesplit[[2]][[1]], len=6)
        checkmate::expect_character(nodesplit[[2]][[1]]$parameter, len=1)
        testthat::expect_equal(names(nodesplit[[2]][[1]]$`overlap matrix`), c("direct", "indirect"))
        checkmate::expect_list(nodesplit[[1]][[1]]$quantiles, len=4, unique=TRUE)
        testthat::expect_equal(names(nodesplit[[2]][[1]]$quantiles), c("difference", "direct", "indirect", "nma"))


        expect_error(mb.nodesplit(network, comparisons=comp, nodesplit.parameters="all", fun=tloglin(pool.rate="abs", method.rate="common")),
                     "Parameter specified for nodesplit.parameters")
      }
    })

  }

})
