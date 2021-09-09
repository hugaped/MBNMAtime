testthat::context("Testing prepare.functions")

datalist <- list(osteopain, copd, goutSUA_CFBcomb)

################### Testing add_index ################

testthat::test_that("add_index functions correctly", {
  for (i in seq_along(datalist)) {
    testthat::expect_equal(is.character(add_index(datalist[[i]])[["treatments"]]), TRUE)
    testthat::expect_equal(is.numeric(add_index(datalist[[i]])[["data.ab"]]$treatment), TRUE)

    testthat::expect_error(mb.network(datalist[[i]], ref="test"))

    testthat::expect_error(mb.network(datalist[[i]], ref=40))
    testthat::expect_message(mb.network(datalist[[i]], ref=1))
    testthat::expect_message(mb.network(datalist[[i]], ref=NULL))

    testthat::expect_error(mb.network(datalist[[i]], ref="Notarealtreatment"))
  }

  testthat::expect_silent(mb.network(osteopain, ref="Pl_0"))

  network <- mb.network(osteopain, ref=3)
  network <- mb.network(osteopain, ref="Ce_200")
  testthat::expect_equal(network$treatments[1], "Ce_200")

  network <- mb.network(osteopain, ref="Et_90")
  testthat::expect_equal(network$treatments[1], "Et_90")

  #### Test character data frames ####
  testdata <- osteopain
  testdata$treatment <- as.character(testdata$treatment)

  testthat::expect_silent(mb.network(testdata, ref="Lu_200"))


  #### Test that if you remove a treatment but there is still a factor for it, labels are correct
  painnet <- mb.network(osteopain)
  expect_equal("Ox_44" %in% painnet$treatments, TRUE)

  dropID <- unique(osteopain$studyID[osteopain$treatment=="Ox_44"])
  noox <- osteopain[osteopain$studyID!=dropID,]
  painnet <- mb.network(noox)
  expect_equal("Ox_44" %in% painnet$treatments, FALSE)
})



test_that("mb.validate.data functions correctly", {

  # Checks data.ab has required column names
  df.err <- osteopain
  names(df.err)[1] <- "study"
  expect_error(mb.validate.data(df.err))

  # Checks there are no NAs
  df.err <- alog_pcfb
  df.err$y[1] <- NA
  expect_error(mb.validate.data(df.err), "y")

  # Checks that all SEs are positive
  df.err <- alog_pcfb
  df.err$se[1] <- -1
  expect_error(mb.validate.data(df.err), "All SEs must be >0")

  # Checks that studies have baseline measurement (unless change from baseline data is being used)
  df.err <- alog_pcfb
  expect_warning(mb.validate.data(df.err, CFB=FALSE))
  expect_silent(mb.validate.data(df.err, CFB = TRUE))

  df.err <- osteopain
  expect_silent(mb.validate.data(df.err))

  df.err <- osteopain[!(osteopain$studyID=="Baerwald 2010" & osteopain$time==0),]
  expect_silent(mb.validate.data(df.err, CFB=TRUE))
  expect_warning(mb.validate.data(df.err, CFB=FALSE))

  # Check that studies have >1 time point (unless CFB=TRUE)
  df.err <- osteopain[!(osteopain$studyID=="Baerwald 2010" & osteopain$time>0),]
  expect_warning(mb.validate.data(df.err, CFB=FALSE), "single time point")
  expect_silent(mb.validate.data(df.err, CFB=TRUE))

  # Checks that class codes are consistent within each treatment
  df.err <- goutSUA_CFBcomb
  df.err$class <- as.character(df.err$class)
  df.err$class[df.err$studyID==1201 & df.err$treatment=="Allo_289"] <- "Test"
  expect_error(mb.validate.data(df.err, CFB=TRUE), "Class codes are different")

  # Checks that studies have more than one arm if single.arm==FALSE
  df.err <- alog_pcfb
  df.err <- df.err[!(df.err$studyID==1 & df.err$treatment!="placebo"),]
  expect_error(mb.validate.data(df.err), "single study arm")

  # Checks that arms are balanced and treatment codes consistent at each time point within study
  df.err <- alog_pcfb
  df.err$treatment[df.err$studyID==1 & df.err$time==2][1] <- "alog_6.25"
  expect_error(mb.validate.data(df.err), "consistent treatment codes")

  expect_silent(mb.validate.data(osteopain))
  expect_silent(mb.validate.data(goutSUA_CFBcomb, CFB=TRUE))
  expect_silent(mb.validate.data(obesityBW_CFB, CFB=TRUE))
  expect_silent(mb.validate.data(alog_pcfb, CFB=TRUE))
})




test_that("mb.network functions correctly", {
  expect_message(mb.network(osteopain), "Reference treatment")
  expect_silent(mb.network(osteopain, reference = "Pl_0"))

  expect_error(mb.network(osteopain, reference = "NOTATREAT"), "Reference treatment specified is not")

  expect_silent(mb.network(goutSUA_CFBcomb, reference="RDEA_100", description="TEST"))

  expect_message(mb.network(alog_pcfb), "Reference treatment")

  expect_silent(mb.network(obesityBW_CFB, reference = "orli"))
})





test_that("genspline functions correctly", {
  x <- c(0:100)

  # B-splines
  bs <- genspline(x, spline="bs", knots = 1, degree=1)
  expect_equal(c("matrix", "array"), class(bs))
  expect_equal(ncol(bs), 2)

  bs <- genspline(x, spline="bs", knots = 1, degree=2)
  expect_equal(ncol(bs), 3)

  bs <- genspline(x, spline="bs", knots = 2, degree=2)
  expect_equal(ncol(bs), 4)

  bs <- genspline(x, spline="bs", knots = c(0.1,0.8), degree=2)
  expect_equal(ncol(bs), 4)

  expect_error(genspline(x, spline="bs", knots = -1, degree=2))
  expect_error(genspline(x, spline="bs", knots = c(1,4)), "probs")
  expect_error(genspline(x, spline="bs", knots = 1, degree=-2))

  # Piecewise linear splines
  ls <- genspline(x, spline="ls", knots = 1)
  expect_equal(c("matrix", "array"), class(ls))
  expect_equal(ncol(ls), 2)

  ls <- genspline(x, spline="ls", knots = 3)
  expect_equal(ncol(ls), 4)

  ls <- genspline(x, spline="ls", knots = c(0.1,0.8))
  expect_equal(ncol(ls), 3)

  expect_error(genspline(x, spline="ls", knots = -1))
  expect_error(genspline(x, spline="ls", knots = c(1,4)), "probs")
  expect_error(genspline(x, spline="ls", knots = 5), "unlikely to be able to support it")

  # Natural cubic splines
  ns <- genspline(x, spline="ns", knots = 1, degree=1)
  expect_equal(c("matrix", "array"), class(ns))
  expect_equal(ncol(ns), 2)

  ns <- genspline(x, spline="ns", knots = 1, degree=2)
  expect_equal(ncol(ns), 2) # degree doesn't matter for ns

  ns <- genspline(x, spline="ns", knots = 3)
  expect_equal(ncol(ns), 4)

  ns <- genspline(x, spline="ns", knots = c(0.1,0.8))
  expect_equal(ncol(ns), 3)

  expect_error(genspline(x, spline="ns", knots = -1))
  expect_error(genspline(x, spline="ns", knots = c(1,4)), "probs")
  expect_error(genspline(x, spline="ns", knots = 5), "unlikely to be able to support it")

  # Restricted cubic splines
  rcs <- genspline(x, spline="rcs", knots = 3)
  expect_equal(c("matrix", "array"), class(rcs))
  expect_equal(ncol(rcs), 2)

  rcs <- genspline(x, spline="rcs", knots = 3, degree=2)
  expect_equal(ncol(rcs), 2) # degree doesn't matter for ns

  rcs <- genspline(x, spline="rcs", knots = 4)
  expect_equal(ncol(rcs), 3)

  rcs <- genspline(x, spline="rcs", knots = c(0.1,0.8,0.9))
  expect_equal(ncol(rcs), 2)

  expect_error(genspline(x, spline="rcs", knots = 1, degree=1))
  expect_error(genspline(x, spline="rcs", knots = -1))
  expect_error(genspline(x, spline="rcs", knots = c(1,4,7)), "probs")
  expect_error(genspline(x, spline="rcs", knots = 7), "unlikely to be able to support it")

})
