testthat::context("Testing rank.functions")

datalist <- list(osteopain=osteopain, copd=copd, goutSUA_CFBcomb=goutSUA_CFBcomb,
                 hyalarthritis=hyalarthritis, diabetes=diabetes, alog_pcfb=alog_pcfb)

testthat::test_that("rank.functions tests pass correctly", {

  testthat::expect_equal(1,1) # Avoids empty tests

  skip_on_ci()
  skip_on_cran()
  skip_on_appveyor()


  for (i in seq_along(datalist)) {

    print(names(datalist)[i])

    network <- mb.network(datalist[[i]])

    emax <- mb.run(network,
                   fun=temax(pool.emax="rel", method.emax="common",
                             pool.et50="rel", method.et50="random",
                             pool.hill="abs", method.hill=2),
                   pd="pv", n.iter=1000)

    if ("n" %in% names(network$data.ab) & !any(is.na(network$data.ab[["n"]]))) {
      bs <- mb.run(network,
                   fun=tspline(type = "bs", degree=2, knots = 2,
                               pool.2="abs", pool.3 = "abs", method.3="random"), pd="pv", link="smd")
    } else {
      bs <- mb.run(network,
                   fun=tspline(type = "bs", degree=2, knots = 2,
                               pool.2="abs", pool.3 = "abs", method.3="random"), pd="pv")
    }


    resdev <- mb.run(network, fun=tpoly(degree=1), parameters.to.save = "resdev", n.iter=1000, pd="pv")


    ############# Rank AUC ###########

    testthat::test_that(paste0(names(datalist)[i], ": rankauc functions correctly"), {

      model.list <- list(emax, bs)
      treats.list <- list(c(2:3),
                          c(1:3))
      int.list <- list(c(0,10), c(1,3))
      subs.list <- list(10, 40)
      dec.list <- list(TRUE, FALSE)
      for (i in seq_along(model.list)) {
        auc <- MBNMAtime:::rankauc(model.list[[i]], decreasing=dec.list[[i]], treats=model.list[[i]]$network$treatments[treats.list[[i]]],
                                   int.range=int.list[[i]], subdivisions=subs.list[[i]], n.iter=100)

        testthat::expect_equal(names(auc), c("summary", "prob.matrix", "rank.matrix", "auc.int"))
        testthat::expect_equal(nrow(auc[["summary"]]), length(treats.list[[i]]))
        testthat::expect_equal(nrow(auc[["prob.matrix"]]), ncol(auc[["prob.matrix"]]))
        testthat::expect_equal(nrow(auc[["prob.matrix"]]), length(treats.list[[i]]))
        testthat::expect_equal(nrow(auc[["rank.matrix"]]), 100)
        testthat::expect_equal(colnames(auc[["rank.matrix"]]), model.list[[i]]$network$treatments[treats.list[[i]]])
      }

      i <- 1
      testthat::expect_error(MBNMAtime:::rankauc(model.list[[i]], decreasing=5, treats=treats.list[[i]],
                                                 int.range=int.list[[i]], subdivisions=subs.list[[i]], n.iter=100))

      testthat::expect_error(MBNMAtime:::rankauc(model.list[[i]], decreasing=dec.list[[i]], treats=c("Placecbo", "Celebrex"),
                                                 int.range=int.list[[i]], subdivisions=subs.list[[i]], n.iter=100))

      testthat::expect_error(MBNMAtime:::rankauc(model.list[[i]], decreasing=dec.list[[i]], treats=treats.list[[i]],
                                                 int.range=c(1:10), subdivisions=subs.list[[i]], n.iter=100))

      i <- 2
      testthat::expect_error(MBNMAtime:::rankauc(model.list[[i]], decreasing=dec.list[[i]], treats=treats.list[[i]],
                                                 int.range=c(-5,5), subdivisions=subs.list[[i]], n.iter=100))

      testthat::expect_error(MBNMAtime:::rankauc(model.list[[i]], decreasing=dec.list[[i]], treats=treats.list[[i]],
                                                 subdivisions=subs.list[[i]], n.iter=100))

      testthat::expect_error(MBNMAtime:::rankauc(model.list[[i]], decreasing=dec.list[[i]], treats=treats.list[[i]],
                                                 int.range=int.list[[i]], subdivisions=-10, n.iter=100))

      # Error due to wrong parameters being saved from model
      testthat::expect_error(MBNMAtime:::rankauc(resdev, decreasing=dec.list[[i]],
                                                 treats=treats.list[[i]],
                                                 int.range=int.list[[i]],
                                                 subdivisions=subs.list[[i]], n.iter=100))

    })



    ############# rank.mbnma #############

    testthat::test_that(paste0(names(datalist)[i], ": rank.mbnma functions correctly"), {

      model.list <- list(emax, bs)
      treats.list <- list(c(1,2,3), network$treatments[c(1,3)])

      i <- 1

      rank <- rank(emax, params=c("emax", "et50"),
                   direction=-1, treats=treats.list[[i]])

      testthat::expect_equal(names(rank), c("emax", "et50"))
      testthat::expect_equal(names(rank$et50), c("summary", "prob.matrix", "rank.matrix", "cum.matrix", "lower_better"))
      testthat::expect_equal(nrow(rank$et50[["summary"]]), length(treats.list[[i]]))
      testthat::expect_equal(nrow(rank$emax[["prob.matrix"]]), ncol(rank$emax[["prob.matrix"]]))
      testthat::expect_equal(nrow(rank$emax[["prob.matrix"]]), length(treats.list[[i]]))
      testthat::expect_equal(nrow(rank$emax[["rank.matrix"]]), model.list[[i]]$BUGSoutput$n.sims)

      # Check that treatment codes can be character or numeric when estimating AUC
      expect_error(rank(emax, params="auc",
                        direction=-1, treats=treats.list[[1]], n.iter=100), NA)
      expect_error(rank(emax, params="auc",
                        direction=-1, treats=c("Badgers"), n.iter=100), "includes treatments/classes not included")

      if (is.numeric(treats.list[[i]])) {
        matchtreat <- emax$network$treatments[treats.list[[i]]]
      } else if (is.character(treats.list[[i]])) {
        matchtreat <- treats.list[[i]]
      }
      testthat::expect_equal(colnames(rank$et50[["rank.matrix"]]), matchtreat)




      i <- 2
      rank <- rank(bs, params=c("d.4", "d.1"),
                   direction=-1, treats=treats.list[[i]])

      testthat::expect_equal(sort(names(rank)), sort(c("d.4", "d.1")))
      testthat::expect_equal(names(rank$d.1), c("summary", "prob.matrix", "rank.matrix", "cum.matrix", "lower_better"))
      testthat::expect_equal(nrow(rank$d.1[["summary"]]), length(treats.list[[i]]))
      testthat::expect_equal(nrow(rank$d.1[["prob.matrix"]]), ncol(rank$d.1[["prob.matrix"]]))
      testthat::expect_equal(nrow(rank$d.4[["prob.matrix"]]), length(treats.list[[i]]))
      testthat::expect_equal(nrow(rank$d.4[["rank.matrix"]]), model.list[[i]]$BUGSoutput$n.sims)

      if (is.numeric(treats.list[[i]])) {
        matchtreat <- bs$network$treatments[treats.list[[i]]]
      } else if (is.character(treats.list[[i]])) {
        matchtreat <- treats.list[[i]]
      }
      testthat::expect_equal(colnames(rank$d.4[["rank.matrix"]]), matchtreat)
      testthat::expect_equal(colnames(rank$d.1[["rank.matrix"]]), matchtreat)


      expect_error(rank(bs, params=c("beta.2"),
                        direction=-1, treats=treats.list[[i]]), "does not vary by treatment")


      # Class effect models
      if ("classes" %in% names(network)) {
        fpoly <- mb.run(network, fun=tfpoly(degree=2),
                        class.effect = list("beta.2"="random"), pd="pv",
                        rho="dunif(0,1)", n.iter=1000)

        testthat::expect_error(rank(fpoly,
                                    direction=-1, param="D.2", treats=c("1","wer"), level="class"), "classes not included")

        testthat::expect_error(rank(fpoly,
                                    direction=-1, treats=fpoly$network$classes[c(2,3)], param="auc", level="class"), NA)

        testthat::expect_silent(rank(fpoly,
                                     direction=-1, treats=c(1,2), param="D.2", level="class"))

        expect_silent(rank(fpoly,
                           direction=-1, treats=fpoly$network$classes[c(2,3)], param="D.2", level="class"))
      }

      expect_error(rank(emax, direction=-1, level="class"), "not a class effect model")
    })




    ################ rank.mb.predict ###############

    testthat::test_that(paste0(names(datalist)[i], ": rank.mb.predict functions correctly"), {

      preds <- predict(emax, E0=7,
                       ref.resp=list(emax=~rnorm(n, -0.5, 0.05), et50=-0.2))

      ranks <- rank(preds, lower_better=TRUE, treat=emax$network$treatments[1:3])
      expect_equal(length(ranks[[1]]), 3)
      expect_equal(ranks[[1]]$summary$treatment, emax$network$treatments[1:3])
      expect_error(plot(ranks), NA)

      preds <- predict(bs)
      ranks <- rank(preds, lower_better=FALSE)
      expect_equal(length(ranks[[1]]), 3)
      expect_equal(ranks[[1]]$summary$treatment, bs$network$treatments)
      expect_error(plot(ranks), NA)

      expect_error(rank(preds, lower_better=TRUE, time=preds$times), "Must have length 1")

    })

  }

})


