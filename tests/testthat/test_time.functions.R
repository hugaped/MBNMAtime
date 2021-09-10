testthat::context("Testing time.functions")

testthat::test_that("tloglin functions correctly", {
  timefun <- tloglin(pool.rate="rel", method.rate="common")
  expect_equal(timefun$nparam, 1)
  expect_equal(timefun$apool, c("rate"="rel"))
  expect_equal(timefun$amethod, c("rate"="common"))
  expect_equal(timefun$name, "loglin")

  timefun <- tloglin(pool.rate="abs", method.rate="random")
  expect_equal(timefun$nparam, 1)
  expect_equal(timefun$apool, c("rate"="abs"))
  expect_equal(timefun$amethod, c("rate"="random"))
  expect_equal(timefun$name, "loglin")
})


testthat::test_that("temax functions correctly", {
  timefun <- temax(pool.emax="rel", method.emax="common", pool.et50="rel", method.et50="common")
  expect_equal(timefun$nparam, 2)

  timefun <- temax(pool.emax="rel", method.emax="common", pool.et50="rel", method.et50="common",
                  pool.hill="rel", method.hill="common")
  expect_equal(timefun$nparam, 3)

  expect_message(temax(pool.emax="rel", method.emax="common", pool.et50="rel", method.et50="common"), "et50")
  expect_message(temax(pool.emax="rel", method.emax="common", pool.et50="rel", method.et50="common",
                       pool.hill="abs", method.hill="random"), "hill")

  timefun <- temax(pool.emax="abs", method.emax="random", pool.et50="rel", method.et50="random",
                   pool.hill="abs", method.hill="common")
  expect_equal(timefun$apool, c(emax="abs", et50="rel", hill="abs"))
  expect_equal(timefun$amethod, c(emax="random", et50="random", hill="common"))

})
