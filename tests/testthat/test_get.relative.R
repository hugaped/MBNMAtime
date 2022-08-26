testthat::context("Testing get.relative")


testthat::test_that("get.relative tests pass correctly", {
  skip_on_ci()
  skip_on_cran()

  testthat::test_that("get.relative passes", {

    alognet <- mb.network(alog_pcfb)
    copdnet <- mb.network(copd)

    itp <- mb.run(copdnet, titp(), corparam = FALSE, n.iter=2000)
    loglin <- mb.run(alognet, tloglin(), n.iter=2000)

    expect_error(get.relative(mbnma=loglin, mbnma.add=itp, time=20),
                 "mbnma and mbnma.add must have a single treatment")


    copdnew <- copd
    levels(copdnew$treatment)[1] <- "placebo"
    copdnet <- mb.network(copdnew)


    itp <- mb.run(copdnet, titp(), corparam = FALSE)
    loglin <- mb.run(alognet, tloglin(), n.iter=2000)

    expect_error(get.relative(mbnma=loglin, mbnma.add=itp), NA)

    expect_error(get.relative(mbnma=loglin, mbnma.add=itp, time=200), NA)

    treats <- c("placebo", "alog_100", "Aclidinium")
    rels <- get.relative(mbnma=loglin, mbnma.add=itp, treats=treats)
    expect_equal(treats, rownames(rels$mean))

    treats <- c("alog_100", "placebo", "Aclidinium")
    rels <- get.relative(mbnma=loglin, mbnma.add=itp, treats=treats)
    expect_equal(c("placebo", "alog_100", "Aclidinium"), rownames(rels$mean))

    expect_error(get.relative(mbnma=loglin, mbnma.add=itp, treats=c("alog_100", "Aclidinium")),
                 "mbnma and mbnma.add must have a single treatment")


    # Test performing MBNMA with a different reference treatment for alog and check again
    alognet2 <- mb.network(alog_pcfb, reference="alog_100")
    loglin2 <- mb.run(alognet2, tloglin(), n.iter=2000)

    treats <- c("alog_100", "placebo", "Aclidinium")
    rels <- get.relative(mbnma=loglin2, mbnma.add=itp, treats=treats)
    expect_equal(rownames(rels$mean), c("placebo", "alog_100", "Aclidinium"))

  })

})
