testthat::context("Testing get.relative")

datalist <- list(osteopain=osteopain, copd=copd, goutSUA_CFBcomb=goutSUA_CFBcomb,
                 hyalarthritis=hyalarthritis, diabetes=diabetes, alog_pcfb=alog_pcfb)

n.iter <- 2000
seed <- 890421

for (i in 2:length(datalist)) {

  print(names(datalist)[i])

  network <- mb.network(datalist[[i]])

  testthat::test_that(paste0(names(datalist)[i], ": get.relative tests pass correctly"), {

    skip_on_ci()
    skip_on_cran()
    skip_on_appveyor()

    if (names(datalist)[i] %in% c("goutSUA_CFBcomb", "hyalarthritis", "alog_pcfb")) {
      itp <- mb.run(network, tpoly(degree=2), corparam = FALSE, n.iter=n.iter, jags.seed=seed)
    } else {
      itp <- mb.run(network, titp(), corparam = FALSE, n.iter=n.iter, jags.seed=seed)
    }

    loglin <- mb.run(mb.network(datalist[[i-1]]), tloglin(), n.iter=n.iter, jags.seed=seed)

    expect_error(get.relative(mbnma=loglin, mbnma.add=itp, time=20),
                 "mbnma and mbnma.add must have a single treatment")


    # Create new network with same treatment
    netnew <- datalist[[i-1]]
    if (class(datalist[[i-1]]$treatment) != class(datalist[[i]]$treatment)) {

      netnew <- datalist[[i-2]]

    }

    if (is.factor(netnew$treatment)) {
      levels(netnew$treatment)[1] <- itp$network$treatments[1]
    } else if (is.character(netnew$treatment)) {
      netnew$treatment[netnew$treatment==netnew$treatment[1]] <-
        itp$network$treatments[1]
    } else if (is.numeric(netnew$treatment)) {
      netnew$treatment[netnew$treatment==netnew$treatment[1]] <-
        as.numeric(itp$network$treatments[1])
    }


    netnew <- mb.network(netnew)

    loglin <- mb.run(netnew, tloglin(), n.iter=n.iter, jags.seed=seed)


    expect_error(get.relative(mbnma=loglin, mbnma.add=itp), NA)

    expect_error(get.relative(mbnma=loglin, mbnma.add=itp, time=200), NA)

    treats <- c(network$treatments[1:2], netnew$treatments[3])
    rels <- get.relative(mbnma=loglin, mbnma.add=itp, treats=treats)
    expect_equal(any(is.na(match(treats, rownames(rels$mean)))), FALSE)

    treats <- c(netnew$treatments[3], network$treatments[c(1,3)])
    rels <- get.relative(mbnma=loglin, mbnma.add=itp, treats=treats)
    expect_equal(c(network$treatments[1],
                   netnew$treatments[3],
                   network$treatments[3]),
                 rownames(rels$mean))

    expect_error(get.relative(mbnma=loglin, mbnma.add=itp,
                              treats=c(network$treatments[2], netnew$treatments[3])),
                 "mbnma and mbnma.add must have a single treatment")


    # Test performing MBNMA with a different reference treatment for alog and check again
    ref <- ifelse(is.numeric(datalist[[i]]$treatment),
                  as.numeric(network$treatments[3]),
                  network$treatments[3])

    netref <- mb.network(datalist[[i]], reference=ref)
    loglin2 <- mb.run(netref, tloglin(), n.iter=n.iter, jags.seed=seed)


    # Create new network with same treatment
    netnew <- datalist[[i-1]]
    if (class(datalist[[i-1]]$treatment) != class(datalist[[i]]$treatment)) {
      netnew <- datalist[[i-2]]
    }

    if (is.factor(netnew$treatment)) {
      levels(netnew$treatment)[1] <- loglin2$network$treatments[1]
    } else if (is.character(netnew$treatment)) {
      netnew$treatment[netnew$treatment==netnew$treatment[1]] <-
        loglin2$network$treatments[1]
    } else if (is.numeric(netnew$treatment)) {
      netnew$treatment[netnew$treatment==netnew$treatment[1]] <-
        as.numeric(loglin2$network$treatments[1])
    }

    netnew <- mb.network(netnew)

    if (names(datalist)[i-1] %in% c("diabetes")) {

      # WARNING CAN BE REMOVED AFTER v0.2.2
      itp2 <- suppressWarnings(mb.run(netnew, temax(), corparam = TRUE, n.iter=n.iter, jags.seed=seed))
    } else {
      itp2 <- mb.run(netnew, titp(), corparam = TRUE, n.iter=n.iter, jags.seed=seed)
    }


    treats <- c(netnew$treatments[3], netref$treatments[c(1,3)])
    rels <- get.relative(mbnma=loglin2, mbnma.add=itp2, treats=treats)
    expect_equal(any(is.na(match(treats, rownames(rels$mean)))), FALSE)

  })

}


