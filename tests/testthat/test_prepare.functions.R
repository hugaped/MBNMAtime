testthat::context("Testing prepare.functions")
data <- osteopain

################### Testing add_index ################

testthat::test_that("add_index functions correctly", {
  testthat::expect_equal(is.character(add_index(data)[["treatments"]]), TRUE)
  testthat::expect_equal(is.numeric(add_index(data)[["data.ab"]]$treatment), TRUE)

  testthat::expect_error(mb.network(data, ref="test"))

  testthat::expect_error(mb.network(data, ref=40))

  testthat::expect_silent(mb.network(data, ref="Pl_0"))
  testthat::expect_message(mb.network(data, ref=1))

  testthat::expect_message(mb.network(data, ref=NULL))

  network <- mb.network(data, ref=3)
  network <- mb.network(data, ref="Ce_200")
  testthat::expect_equal(network$treatments[1], "Ce_200")

  network <- mb.network(data, ref="Et_90")
  testthat::expect_equal(network$treatments[1], "Et_90")

  #### Test character data frames ####
  testdata <- data
  testdata$treatment <- as.character(testdata$treatment)

  testthat::expect_error(mb.network(testdata, ref=1))
  testthat::expect_error(mb.network(testdata, ref="Notarealtreatment"))
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




test_that("jagstonetwork functions correctly", {
  network <- mb.network(osteopain)
  emax <- mb.emax(network, n.iter=500)

  df  <- jagstonetwork(emax)
  expect_identical(sort(names(df)), sort(c("studyID", "treatment", "time", "y", "se", "arm", "fupcount", "fups", "narm")))
  expect_equal(nrow(df)>0, TRUE)

  network <- mb.network(goutSUA_CFBcomb)
  emax <- mb.emax(network, class.effect = list(emax="random"), n.iter=500)
  df  <- jagstonetwork(emax)
  expect_equal("class" %in% names(df), TRUE)
  expect_equal(nrow(df)>0, TRUE)
  expect_equal(any(is.na(df$class)), FALSE)
})
