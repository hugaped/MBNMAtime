testthat::context("Testing write.functions")
network <- mb.network(osteopain)

emax1 <- mb.emax(network,
                    emax=list(pool="rel", method="random"),
                    et50=list(pool="rel", method="common"),
                    positive.scale=TRUE,
                    n.chain=3, n.iter=200, n.burnin=100)


################### Testing time.fun ################

testthat::test_that("time.fun functions correctly", {
  timecourse <- time.fun(fun="exponential", alpha="study")
  testthat::expect_equal(grepl("alpha\\[i\\]", timecourse$jagscode), TRUE)
  testthat::expect_equal(grepl("exp", timecourse$relationship), TRUE)

  timecourse <- time.fun(fun="user", user.fun="alpha+(1/beta.1*time) + (beta.2^beta.3)", beta.2="rel.common", beta.3=2, alpha="study")
  testthat::expect_equal(grepl("alpha\\[i\\]", timecourse$jagscode), TRUE)
  testthat::expect_equal(grepl("time\\[i,m\\]", timecourse$jagscode), TRUE)
  testthat::expect_equal(grepl("beta.3\\[", timecourse$jagscode), FALSE)

  # Log linear user.fun
  timecourse <- time.fun(fun="user", user.fun="exp(alpha+(beta.1*time))", beta.1="rel.common", alpha="study")
  testthat::expect_equal(timecourse$time.function, "user")
  testthat::expect_equal(timecourse$jagscode, "exp(alpha[i]+(beta.1[i,k]*time[i,m]))")
  #testthat::expect_equal(timecourse$jagscode, "exp(alpha[i]+(beta.1*time[i,m]))")

  timecourse <- time.fun(fun="emax", alpha="arm", beta.1="const.random", beta.2="rel.random")
  testthat::expect_equal(grepl("beta\\.1", timecourse$jagscode), TRUE)
  testthat::expect_equal(grepl("beta\\.2", timecourse$jagscode), TRUE)
  testthat::expect_equal(grepl("beta\\.3", timecourse$jagscode), FALSE)

  timecourse <- time.fun(fun="piecelinear", alpha="study", beta.1="rel.common", beta.2="rel.random", beta.3="const.random")
  testthat::expect_equal(grepl("time\\[i\\,m\\] < beta\\.3", timecourse$jagscode), TRUE)

  testthat::expect_message(time.fun(fun="emax.hill", user.fun="alpha+(1/beta.1*time) + (beta.2^beta.3)",
                          beta.1="rel.common", beta.2="rel.random", beta.3="const.common", alpha="study"))

  timecourse.fe <- time.fun(fun="fract.poly.first", beta.1="rel.common")
  timecourse.re <- time.fun(fun="fract.poly.first", beta.1="rel.random")
  testthat::expect_equal(timecourse.fe$jagscode, timecourse.re$jagscode)

  testthat::expect_length(time.fun(fun="emax.hill", alpha="arm", beta.1="rel.common", beta.2="rel.random", beta.3="const.random"), 3)

  testthat::expect_silent(time.fun())

})




################### Testing mb.write ################


testthat::test_that("test.mb.write", {
  funlist <- c("linear", "exponential")

  for (i in seq_along(funlist)) {
    print(funlist[i])
    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", intercept=TRUE, positive.scale = TRUE)
    testthat::expect_equal(grepl("alpha\\[i", jags), TRUE)
    testthat::expect_equal(grepl("exp\\(alpha\\[i", jags), TRUE)

    jags <- mb.write(fun=funlist[i], beta.1="rel.common", intercept=FALSE, positive.scale = FALSE)
    testthat::expect_equal(grepl("alpha\\[i", jags), FALSE)

    jags <- mb.write(fun=funlist[i], alpha="arm")
    testthat::expect_equal(grepl("alpha\\[i,k\\] ~ dnorm\\(0,0\\.0001\\)", jags), TRUE)

    jags <- mb.write(fun=funlist[i], alpha="study")
    testthat::expect_equal(grepl("alpha\\[i\\] ~ dnorm\\(0,0\\.0001\\)", jags), TRUE)

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.random")
    testthat::expect_equal(grepl("d\\.1", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.1 ~ dnorm", jags), TRUE)
    testthat::expect_equal(grepl("dev", jags), TRUE)

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", rho="estimate", covar="CS")
    testthat::expect_equal(grepl("rho ~ dunif", jags), TRUE)
    testthat::expect_equal(grepl("cor\\[i,r,c\\]", jags), FALSE)
    testthat::expect_equal(grepl("dev", jags), FALSE)

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", rho="estimate", covar="CS")
    testthat::expect_equal(grepl("rho ~ dunif", jags), TRUE)
    testthat::expect_equal(grepl("cor\\[i,r,c\\]", jags), FALSE)

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", rho="estimate", covar="CS")
    testthat::expect_equal(grepl("rho ~ dunif", jags), TRUE)
    testthat::expect_equal(grepl("cor\\[i,r,c\\]", jags), FALSE)

    testthat::expect_error(mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", rho=0.1))
    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", rho=0.1, covar="AR1")
    testthat::expect_equal(grepl("rho <- 0.1", jags), TRUE)

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", rho="estimate", covar="AR1")
    testthat::expect_equal(grepl("rho ~ dunif\\(-1,1", jags), TRUE)
    testthat::expect_equal(grepl("cor\\[i,r,c\\]", jags), TRUE)
    testthat::expect_equal(grepl("dev", jags), FALSE)

    testthat::expect_error(mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", rho=5, covar="CS"))

    testthat::expect_error(mb.write(fun=funlist[i], beta.1="rel.common", intercept=TRUE, positive.scale = TRUE, class.effect="alpha"))
    testthat::expect_error(mb.write(fun=funlist[i], beta.1="rel.common", intercept=TRUE, positive.scale = TRUE, class.effect=list("fake")))
    testthat::expect_error(mb.write(fun=funlist[i], beta.1="rel.common", intercept=TRUE, positive.scale = TRUE, class.effect=list("alpha"=FALSE)))

    jags <- mb.write(fun=funlist[i], beta.1="rel.common", intercept=TRUE, positive.scale = TRUE, class.effect=list("beta.1"="common"))
    testthat::expect_equal(grepl("# Priors on relative class effects", jags), TRUE)
    testthat::expect_equal(grepl("D\\.1\\[class\\[k\\]\\]\n", jags), TRUE)

    jags <- mb.write(fun=funlist[i], beta.1="rel.common", intercept=TRUE, positive.scale = TRUE, class.effect=list("beta.1"="random"))
    testthat::expect_equal(grepl("D\\.1\\[class\\[k\\]\\]\\, tau\\.D\\.1\\)", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.D\\.1", jags), TRUE)

    jags <- mb.write(fun=funlist[i], beta.1="arm.common", intercept=TRUE, positive.scale = TRUE, class.effect=list("beta.1"="common"))
    testthat::expect_equal(grepl("# Priors on absolute class effects", jags), TRUE)
    testthat::expect_equal(grepl("BETA\\.1\\[class\\[k\\]\\]\n", jags), TRUE)

    jags <- mb.write(fun=funlist[i], beta.1="arm.common", intercept=TRUE, positive.scale = TRUE, class.effect=list("beta.1"="random"))
    testthat::expect_equal(grepl("BETA\\.1\\[class\\[k\\]\\]\\, tau\\.BETA\\.1\\)", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.BETA\\.1", jags), TRUE)

    testthat::expect_error(mb.write(fun=funlist[i], beta.1="rel.common", class.effect=list("beta.1"="random"), UME=c("beta.1")))
  }


  funlist <- c("emax", "quadratic")
  for (i in seq_along(funlist)) {
    print(funlist[i])

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", beta.2="rel.random")
    testthat::expect_equal(grepl("d\\.1", jags), TRUE)
    testthat::expect_equal(grepl("d\\.2", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.2 ~ dnorm", jags), TRUE)

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="arm.common", beta.2="arm.random")
    testthat::expect_equal(grepl("beta\\.1", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.2", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.beta\\.2 ~ dnorm", jags), TRUE)

    jags <- mb.write(fun=funlist[i], alpha="arm", beta.1="rel.common", beta.2="rel.random",
                        rho="estimate", covar="AR1")
    testthat::expect_equal(grepl("d\\.1", jags), TRUE)
    testthat::expect_equal(grepl("d\\.2", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.2 ~ dnorm", jags), TRUE)
    testthat::expect_equal(grepl("rho ~ dunif\\(-1,1", jags), TRUE)
    testthat::expect_equal(grepl("cor\\[i,r,c\\]", jags), TRUE)
    testthat::expect_equal(grepl("dev", jags), FALSE)

    jags <- mb.write(fun=funlist[i], beta.1="rel.common", beta.2="rel.random", intercept=TRUE, positive.scale = TRUE,
                        class.effect=list("beta.1"="common", "beta.2"="random"), rho="estimate", covar="CS")
    testthat::expect_equal(grepl("# Priors on relative class effects", jags), TRUE)
    testthat::expect_equal(grepl("D\\.1\\[class\\[k\\]\\]\n", jags), TRUE)
    testthat::expect_equal(grepl("D\\.2\\[class\\[k\\]\\]\\, tau\\.D\\.2\\)", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.D\\.2", jags), TRUE)

    jags <- mb.write(fun=funlist[i], beta.1="arm.common", beta.2="arm.random", intercept=TRUE, positive.scale = TRUE,
                        class.effect=list("beta.1"="common", "beta.2"="random"), rho="estimate", covar="CS")
    testthat::expect_equal(grepl("# Priors on absolute class effects", jags), TRUE)
    testthat::expect_equal(grepl("BETA\\.1\\[class\\[k\\]\\]\n", jags), TRUE)
    testthat::expect_equal(grepl("BETA\\.2\\[class\\[k\\]\\]\\, tau\\.BETA\\.2\\)", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.BETA\\.2", jags), TRUE)

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", beta.2="rel.random", UME=c("beta.1"))
    testthat::expect_equal(grepl("d\\.1\\[1\\]", jags), FALSE)
    testthat::expect_equal(grepl("d\\.2\\[1\\]", jags), TRUE)
    testthat::expect_equal(grepl("\\[treat\\[i,1\\],treat\\[i,k\\]\\]", jags), TRUE)
    testthat::expect_equal(grepl("UME priors", jags), TRUE)
  }

  funlist <- c("fract.poly.first")
  for (i in seq_along(funlist)) {
    print(funlist[i])

    testthat::expect_error(mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", beta.3="rel.random"))
    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", beta.3="const.random")
    testthat::expect_equal(grepl("d\\.1", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.3", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.beta\\.3", jags), TRUE)

    jags <- mb.write(fun=funlist[i], alpha="arm", beta.1="rel.common", beta.3="const.random",
                        rho="estimate", covar="AR1")
    testthat::expect_equal(grepl("d\\.1", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.3", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.beta\\.3", jags), TRUE)
    testthat::expect_equal(grepl("rho ~ dunif\\(-1,1", jags), TRUE)
    testthat::expect_equal(grepl("cor\\[i,r,c\\]", jags), TRUE)
    testthat::expect_equal(grepl("dev", jags), FALSE)

    jags <- mb.write(fun=funlist[i], beta.1="rel.common", beta.3=2, intercept=TRUE, positive.scale = TRUE,
                        class.effect=list("beta.1"="common"), rho="estimate", covar="AR1")
    testthat::expect_equal(grepl("# Priors on relative class effects", jags), TRUE)
    testthat::expect_equal(grepl("D\\.1\\[class\\[k\\]\\]\n", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.3", jags), TRUE)

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", beta.3=-1.5, UME=c("beta.1"))
    testthat::expect_equal(grepl("d\\.1\\[1\\]", jags), FALSE)
    testthat::expect_equal(grepl("\\[treat\\[i,1\\],treat\\[i,k\\]\\]", jags), TRUE)
    testthat::expect_equal(grepl("UME priors", jags), TRUE)
  }

  funlist <- c("emax.hill", "piecelinear")
  for (i in seq_along(funlist)) {
    print(funlist[i])

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", beta.2="rel.random", beta.3=5)
    testthat::expect_equal(grepl("d\\.1", jags), TRUE)
    testthat::expect_equal(grepl("d\\.2", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.3", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.2 ~ dnorm", jags), TRUE)

    if (funlist[i]!="piecelinear") {
      jags <- mb.write(fun=funlist[i], alpha="arm", beta.1="rel.common", beta.2="rel.random", beta.3="rel.common",
                          rho="estimate", covar="AR1")
      testthat::expect_equal(grepl("d\\.1", jags), TRUE)
      testthat::expect_equal(grepl("d\\.2", jags), TRUE)
      testthat::expect_equal(grepl("sd\\.2 ~ dnorm", jags), TRUE)
      testthat::expect_equal(grepl("d\\.3", jags), TRUE)
      testthat::expect_equal(grepl("rho ~ dunif\\(-1,1", jags), TRUE)
      testthat::expect_equal(grepl("cor\\[i,r,c\\]", jags), TRUE)
      testthat::expect_equal(grepl("dev", jags), FALSE)
    }

    jags <- mb.write(fun=funlist[i], beta.1="rel.common", beta.2="rel.random", intercept=TRUE, beta.3="const.common", positive.scale = TRUE,
                        class.effect=list("beta.1"="common", "beta.2"="random"), rho="estimate", covar="CS")
    testthat::expect_equal(grepl("# Priors on relative class effects", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.3", jags), TRUE)
    testthat::expect_equal(grepl("D\\.1\\[class\\[k\\]\\]\n", jags), TRUE)
    testthat::expect_equal(grepl("D\\.2\\[class\\[k\\]\\]\\, tau\\.D\\.2\\)", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.D\\.2", jags), TRUE)

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", beta.2="rel.random", beta.3="const.random", UME=c("beta.1"))
    testthat::expect_equal(grepl("d\\.1\\[1\\]", jags), FALSE)
    testthat::expect_equal(grepl("d\\.2\\[1\\]", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.beta\\.3", jags), TRUE)
    testthat::expect_equal(grepl("\\[treat\\[i,1\\],treat\\[i,k\\]\\]", jags), TRUE)
    testthat::expect_equal(grepl("UME priors", jags), TRUE)

    if (funlist[i]=="piecelinear") {
      testthat::expect_equal(grepl("beta\\.3 ~ dunif\\(0, maxtime\\)", jags), TRUE)

      jags <- mb.write(fun=funlist[i], alpha="arm",
                          beta.1="const.common", beta.2="rel.common", beta.3="const.random")
      testthat::expect_equal(grepl("\\(s\\.beta\\.1\\[i\\,k\\] \\* s\\.beta\\.3\\[i\\,k\\]\\)", jags), TRUE)
      testthat::expect_equal(grepl("sd\\.beta\\.3", jags), TRUE)

      jags <- mb.write(fun=funlist[i], alpha="arm",
                          beta.1="rel.common", beta.2="rel.common", beta.3="const.common")
      testthat::expect_equal(grepl("\\(beta\\.1\\[i\\,k\\] \\* s\\.beta\\.3\\[i\\,k\\]\\)", jags), TRUE)
      testthat::expect_equal(grepl("sd\\.beta\\.3", jags), FALSE)

    }

    testthat::expect_error(mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", beta.2="rel.random", UME=c("beta.1")))
  }

  funlist <- c("fract.poly.second")
  for (i in seq_along(funlist)) {
    print(funlist[i])

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", beta.2="rel.random", beta.3=5, beta.4=1)
    testthat::expect_equal(grepl("d\\.1", jags), TRUE)
    testthat::expect_equal(grepl("d\\.2", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.2 ~ dnorm", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.3", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.4", jags), TRUE)

    jags <- mb.write(fun=funlist[i], alpha="arm", beta.1="rel.common", beta.2="rel.random",
                        beta.3="const.random", beta.4="const.common",
                        rho="estimate", covar="AR1")
    testthat::expect_equal(grepl("d\\.1", jags), TRUE)
    testthat::expect_equal(grepl("d\\.2", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.2 ~ dnorm", jags), TRUE)
    testthat::expect_equal(grepl("rho ~ dunif\\(-1,1", jags), TRUE)
    testthat::expect_equal(grepl("cor\\[i,r,c\\]", jags), TRUE)
    testthat::expect_equal(grepl("dev", jags), FALSE)
    testthat::expect_equal(grepl("beta\\.3", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.beta\\.3", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.4", jags), TRUE)

    jags <- mb.write(fun=funlist[i], beta.1="rel.common", beta.2="rel.random", beta.3="const.common", beta.4="const.random",
                        intercept=TRUE, positive.scale = TRUE,
                        class.effect=list("beta.1"="common", "beta.2"="random"), rho="estimate", covar="CS")
    testthat::expect_equal(grepl("# Priors on relative class effects", jags), TRUE)
    testthat::expect_equal(grepl("D\\.1\\[class\\[k\\]\\]\n", jags), TRUE)
    testthat::expect_equal(grepl("D\\.2\\[class\\[k\\]\\]\\, tau\\.D\\.2\\)", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.D\\.2", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.3", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.4", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.beta\\.4", jags), TRUE)

    jags <- mb.write(fun=funlist[i], alpha="study", beta.1="rel.common", beta.2="rel.random",
                        beta.3="const.random", beta.4="const.common", UME=c("beta.1"))
    testthat::expect_equal(grepl("d\\.1\\[1\\]", jags), FALSE)
    testthat::expect_equal(grepl("d\\.2\\[1\\]", jags), TRUE)
    testthat::expect_equal(grepl("\\[treat\\[i,1\\],treat\\[i,k\\]\\]", jags), TRUE)
    testthat::expect_equal(grepl("UME priors", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.3", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.beta\\.3", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.4", jags), TRUE)

  }



  funlist <- c("user")
  for (i in seq_along(funlist)) {
    print(funlist[i])

    user.fun <- "exp(alpha + beta.1*time)"
    jags <- mb.write(fun=funlist[i], user.fun=user.fun, alpha="study", beta.1="rel.random")
    testthat::expect_equal(grepl("d\\.1", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.1 ~ dnorm", jags), TRUE)
    testthat::expect_equal(grepl("alpha\\[i\\]", jags), TRUE)

    user.fun <- "alpha + ((beta.1^beta.2) * time) + (beta.3 * time)"
    jags <- mb.write(fun=funlist[i], user.fun=user.fun,
                        alpha="arm", beta.1="rel.common", beta.2="const.random",
                        beta.3="rel.random",
                        rho="estimate", covar="AR1")
    testthat::expect_equal(grepl("d\\.1", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.2", jags), TRUE)
    testthat::expect_equal(grepl("d\\.3", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.3 ~ dnorm", jags), TRUE)
    testthat::expect_equal(grepl("rho ~ dunif\\(-1,1", jags), TRUE)
    testthat::expect_equal(grepl("cor\\[i,r,c\\]", jags), TRUE)
    testthat::expect_equal(grepl("dev", jags), FALSE)

    jags <- mb.write(fun=funlist[i], user.fun=user.fun,
                        alpha="arm", beta.1="rel.common", beta.2="const.common",
                        beta.3="rel.random",
                        intercept=TRUE, positive.scale = TRUE,
                        class.effect=list("beta.1"="common", "beta.3"="random"),
                        rho="estimate", covar="CS")
    testthat::expect_equal(grepl("# Priors on relative class effects", jags), TRUE)
    testthat::expect_equal(grepl("D\\.1\\[class\\[k\\]\\]\n", jags), TRUE)
    testthat::expect_equal(grepl("D\\.3\\[class\\[k\\]\\]\\, tau\\.D\\.3\\)", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.D\\.3", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.2", jags), TRUE)

    jags <- mb.write(fun=funlist[i], user.fun=user.fun,
                        alpha="arm", beta.1="rel.common", beta.2="const.random",
                        beta.3="rel.random", UME=TRUE)
    testthat::expect_equal(grepl("d\\.1\\[1\\]", jags), FALSE)
    testthat::expect_equal(grepl("d\\.3\\[c,k\\]", jags), TRUE)
    testthat::expect_equal(grepl("\\[treat\\[i,1\\],treat\\[i,k\\]\\]", jags), TRUE)
    testthat::expect_equal(grepl("UME priors", jags), TRUE)
    testthat::expect_equal(grepl("beta\\.2", jags), TRUE)
    testthat::expect_equal(grepl("sd\\.beta\\.2", jags), TRUE)

    jags <- mb.write(fun=funlist[i], user.fun=user.fun,
                        alpha="arm", beta.1="rel.common", beta.2="const.common",
                        beta.3="rel.random", UME="beta.3")
    testthat::expect_equal(grepl("d\\.1\\[1\\]", jags), TRUE)
    testthat::expect_equal(grepl("d\\.3\\[c,k\\]", jags), TRUE)


    # Errors for user.fun
    testthat::expect_error(mb.write(fun=funlist[i], beta.1="rel.common",
                             beta.2="const.common", beta.3="rel.random"))
    testthat::expect_error(mb.write(fun=funlist[i], user.fun=user.fun,
                             alpha="arm", beta.1="rel.common", beta.2="const.common",
                             beta.3="rel.random",
                             class.effect=list("beta.1"="common", "beta.2"="random")))
    testthat::expect_error(mb.write(fun=funlist[i], user.fun=user.fun,
                             alpha="arm", beta.1="rel.common", beta.2="const.common",
                             beta.3="rel.random",
                             intercept=FALSE, positive.scale = FALSE,
                             UME=c("beta.1", "beta.2")))

  }


  # expect error
  testthat::expect_error(mb.write(fun=funlist[i], alpha="golf", beta.1="rel.common", intercept=TRUE, positive.scale = TRUE))
  testthat::expect_error(mb.write(fun="linear", alpha="study", beta.1="rel.common", beta.2="rel.common", intercept=TRUE, positive.scale = TRUE))
  testthat::expect_error(mb.write(fun="quadratic", alpha="study", beta.1="rel.common", beta.3="const.common", intercept=TRUE, positive.scale = TRUE))
  testthat::expect_error(mb.write(fun="quadratic", alpha="study", beta.1="rel.common", beta.2="rel.random", beta.3="const.common", intercept=TRUE, positive.scale = TRUE))
  testthat::expect_error(mb.write(fun="piecelinear", alpha="study", beta.1="rel.common", beta.2="rel.random", beta.3="test", intercept=TRUE, positive.scale = TRUE))

  #testthat::expect_error(mb.write(fun="emax", beta.1="rel.common", beta.2="rel.random", intercept=TRUE, positive.scale = TRUE,
  #                         class.effect=list("beta.1"="common", "beta.2"="random", "beta.3"="random"))
  #)

  testthat::expect_silent(mb.write())
})






############### Testing prior writing functions ###############

testthat::test_that("get.prior", {
  testthat::expect_equal(sort(names(get.prior(emax1$model.arg$jagscode))), sort(c("mu.2", "mu.1", "alpha", "sd.emax", "inv.R")))
  testthat::expect_equal(class(get.prior(emax1$model.arg$jagscode)), "list")
  testthat::expect_equal(class(get.prior(emax1$model.arg$jagscode)[[1]]), "character")

  priors <- get.prior(emax1$model.arg$jagscode)
  for (i in seq_along(priors)) {
    expect_equal(grepl("d[a-z]+\\(.+\\)", priors[i]), TRUE)
  }

  testthat::expect_error(get.prior(5))
})



testthat::test_that("replace.prior", {
  priors <- list("alpha"="dnorm(-1, 0.01)",
                 "sd.emax"="dnorm(0,0.5) T(0,)",
                 "mu.1"="dnorm(10,0.1)")

  testthat::expect_equal(grepl("T\\(0,\\)", replace.prior(priors, mbnma=emax1)), TRUE)
  testthat::expect_equal(grepl("mu\\.1\\[i\\] ~ dnorm\\(10,0\\.1\\)", replace.prior(priors, mbnma=emax1)), TRUE)

  testthat::expect_equal(grepl("T\\(0,\\)", replace.prior(priors, model=emax1$model.arg$jagscode)), TRUE)
  testthat::expect_equal(grepl("mu\\.1\\[i\\] ~ dnorm\\(10,0\\.1\\)", replace.prior(priors, model=emax1$model.arg$jagscode)), TRUE)

  priors <- list("banana"="dnorm(-1,0.01)")

  testthat::expect_error(replace.prior(priors, mbnma=emax1))

})
