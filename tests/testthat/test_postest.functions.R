testthat::context("Testing predict.functions")
network <- mb.network(osteopain)

# exponential <- mb.exponential(network, lambda=list(pool="rel", method="common"),
#                                  positive.scale=TRUE,
#                                  n.chain=3, n.iter=1200, n.burnin=800)
#
# emax1 <- mb.emax(network, emax=list(pool="rel", method="random"),
#                     et50=list(pool="rel", method="common"),
#                     positive.scale=TRUE,
#                     n.chain=3, n.iter=1000, n.burnin=600)
#
# emax2 <- mb.emax(network, emax=list(pool="rel", method="common"),
#                     et50=list(pool="const", method="common"),
#                     positive.scale=TRUE,
#                     n.chain=3, n.iter=1000, n.burnin=600)
#
# emax.hill <- mb.emax.hill(network, emax=list(pool="arm", method="common"),
#                              et50=list(pool="arm", method="random"),
#                              hill=list(pool="const", method=2),
#                              positive.scale=TRUE, n.chain=3, n.iter=1000, n.burnin=600)
#
#
# piece <- mb.piecelinear(network, slope.1=list(pool="arm", method="common"),
#                            slope.2=list(pool="rel", method="common"),
#                            knot=list(pool="const", method=1),
#                            positive.scale=TRUE, n.chain=3, n.iter=1000, n.burnin=600)
#
# fract.first <- mb.fract.first(network, slope=list(pool="rel", method="random"),
#                            power=list(pool="const", method="common"),
#                            positive.scale=TRUE, n.chain=3, n.iter=1000, n.burnin=600)
#
#
# # Specify a test model in terms of betas (use user.fun once it is fixed)
# quad <- mb.run(network, fun="quadratic",
#                   beta.1=list(pool="rel", method="common"),
#                   beta.2=list(pool="arm", method="random"),
#                   positive.scale=TRUE,
#                   alpha="study",
#                   n.chain=3, n.iter=1200, n.burnin=800)
#
# # Class data
# #classdata <- osteopain
# # classdata$class[classdata$treatment=="Pl_0"] <- 1
# # classdata$class[classdata$treatment!="Pl_0"] <- 2
# classdata <- goutSUA_CFBcomb
# classnetwork <- mb.network(classdata)
#
# exp.class.fixed <- mb.exponential(classnetwork,
#                                      lambda=list(pool="rel", method="common"),
#                                      positive.scale=TRUE,
#                                n.chain=3, n.iter=1200, n.burnin=800,
#                                class.effect=list("lambda"="common"))
#
# exp.class.random <- mb.exponential(classnetwork,
#                                       lambda=list(pool="rel", method="common"),
#                                       positive.scale=TRUE,
#                                n.chain=3, n.iter=1200, n.burnin=800,
#                                class.effect=list("lambda"="random"))
#
#
# # A model that does not save the required parameters for postestimation
# resdev <- mb.emax(network, emax=list(pool="rel", method="common"),
#                      et50=list(pool="const", method="common"), positive.scale=TRUE,
#                           n.chain=3, n.iter=1200, n.burnin=800,
#                           parameters.to.save=c("resdev")
# )
#
# ################### Testing add_index ################
#
# testthat::test_that("predict.mbnma functions correctly", {
#   model.list <- list(exponential, emax1, emax2, emax.hill, piece, fract.first, quad)
#   treats.list <- list(c(1,5,8,15), c(2,6,9,16:20))
#   ref.resp.list <- list(osteopain[osteopain$treatname=="Placebo_0",],
#                    osteopain[osteopain$treatname=="Celebrex_200",])
#   times.list <- list(c(0:10), c(1,10:20))
#   E0.list <- list(7, "rnorm(n, 7,2)")
#   synth.list <- list("fixed", "random")
#
#   for (i in 1:4) {
#     model <- model.list[[i]]
#     print(paste0("modellist: ", i))
#     for (k in 1:2) {
#       synth <- synth.list[[k]]
#       for (m in 1:2) {
#         E0 <- E0.list[[m]]
#         ref.resp <- ref.resp.list[[m]]
#         treats <- treats.list[[m]]
#         times <- times.list[[m]]
#
#         # Tests using ref.resp
#         testthat::expect_error(predict(model, times=times,
#                                    E0=E0, treats=treats,
#                                    ref.resp=ref.resp, synth=synth
#                                    ), NA)
#
#         pred <- predict(model, times=times,
#                               E0=E0, treats=treats,
#                               ref.resp=ref.resp, synth=synth)
#
#         testthat::expect_equal(length(pred$pred.mat), length(treats))
#         #testthat::expect_equal(names(pred$pred.mat), as.character(treats))
#         testthat::expect_equal(names(pred$pred.mat), model$treatments[treats])
#         testthat::expect_identical(names(pred), c("summary", "pred.mat", "mbnma"))
#         testthat::expect_equal(nrow(pred$pred.mat[[1]]), model$BUGSoutput$n.sims)
#         testthat::expect_equal(nrow(pred$summary[[1]]), length(times))
#         testthat::expect_equal(identical(pred$summary[[1]]$time, times), TRUE)
#       }
#     }
#   }
#
#   # Tests of class models
#   testthat::expect_error(predict(exp.class.random, times=times,
#                              E0=E0, treats=treats,
#                              ref.resp=ref.resp, synth=synth))
#
#   expect_error(predict(exp.class.fixed, times=times,
#                                  E0=E0, treats=treats,
#                                  ref.resp=ref.resp, synth=synth), "Class effects have not")
#
#
#
#   # Tests using ref.resp
#   ref.resp <- list("emax"=-1)
#   testthat::expect_error(predict(emax2, times=times,
#                              E0=7, treats=treats,
#                              ref.resp=ref.resp),
#                NA)
#
#   ref.resp <- list("d.emax"=-1) # incorrect prior name ("d.emax" rather than "emax")
#   testthat::expect_error(predict(emax2, times=times,
#                              E0=7, treats=treats,
#                              ref.resp=ref.resp))
#
#   ref.resp <- list("beta.1"="rnorm(n, -1,1)")
#   testthat::expect_error(predict(quad, times=times,
#                              E0=E0, treats=treats,
#                              ref.resp=ref.resp), NA)
#
#   ref.resp <- list("beta.1"="rnorm(n, -1")
#   testthat::expect_error(predict(quad, times=times,
#                              E0=E0, treats=treats,
#                              ref.resp=ref.resp))
#
#   ref.resp <- list("beta.1"="rnorm(n, -1,1)", "beta.2"="rnorm(n, 1, 0.1)") # beta.2 is not a relative effect in quad
#   testthat::expect_error(predict(quad, times=times,
#                              E0=E0, treats=treats,
#                              ref.resp=ref.resp))
#
#   ref.resp <- list("beta.1"="rnorm(n, -1,1)")
#   testthat::expect_error(predict(quad, times=times,
#                              E0=E0, treats=treats,
#                              ref.resp=ref.resp), NA)
#
#   ref.resp <- list("beta.1"="rnorm(n, -1,1)")
#   testthat::expect_error(predict(emax2, times=times,
#                              E0=E0, treats=treats,
#                              ref.resp=ref.resp))
#
#   ref.resp <- list("emax"="rnorm(n, -1,1)", "et50"="rnorm(n, 1, 0.1)")
#   testthat::expect_error(predict(emax2, times=times,
#                              E0=E0, treats=treats,
#                              ref.resp=ref.resp))
#
#   ref.resp <- list("emax"=-1, "et50"=0.1)
#   testthat::expect_error(predict(emax1, times=times,
#                              E0=E0, treats=treats,
#                              ref.resp=ref.resp),
#                NA)
#
#   ref.resp <- list("slope.2"=0)
#   testthat::expect_error(predict(piece, times=times,
#                                  E0=E0, treats=treats,
#                                  ref.resp=ref.resp),
#                          NA)
#
#   ref.resp <- list("slope"=0)
#   testthat::expect_error(predict(fract.first, times=times,
#                                  E0=E0, treats=treats,
#                                  ref.resp=ref.resp),
#                          NA)
#
#   # Error due to wrong parameters being saved from model
#   testthat::expect_error(predict(resdev, times=times,
#                              E0=E0, treats=treats,
#                              ref.estimate=ref.estimate))
#
#
#   # Expect no error even if ref.resp is NULL
#   testthat::expect_error(predict(emax.hill, times=times,
#                              E0=E0, treats=treats,
#                              ref.resp=NULL), NA)
#
#
#   # Test using different treatment codes (numeric and class)
#   expect_error(predict(emax1, times=c(0:15), treats=c(1,5,10,15)), NA)
#   expect_error(predict(emax1, times=c(0:15), treats=c(1,5,10,15,100)), "numeric treatment codes")
#
#   pred <- predict(emax1, times=c(0:15), treats=c("Pl_0","Du_90","Na_1000"))
#   expect_identical(names(pred$summary), c("Pl_0","Du_90","Na_1000"))
#
#
#   # Test using alogliptin dataset
#   alognet <- mb.network(alog_pcfb)
#   quad <- mb.quadratic(alognet)
#   expect_error(predict(quad), NA)
#
#   quad <- mb.quadratic(alognet, beta.1=list("arm", "common"), beta.2=list("rel", "random"), n.iter=8000)
#   quad <- mb.quadratic(alognet)
#   quad <- mb.quadratic(alognet, beta.2=list("arm", "common"), beta.1=list("rel", "random"), n.iter=8000)
# })



# testthat::test_that("ref.synth functions correctly", {
#   ref.resp <- osteopain[osteopain$treatname=="Placebo_0",]
#
#   testthat::expect_warning(ref.synth(ref.resp, exponential, synth="random",
#                            n.burnin=100, n.iter=200))
#
#   rand <- ref.synth(ref.resp, exponential, synth="random")
#   fix <- ref.synth(ref.resp, exponential, synth="fixed")
#
#   testthat::expect_equal(length(rand), 2)
#   testthat::expect_equal(length(fix), 1)
#
#   rand.emax <- ref.synth(ref.resp, emax1, synth="random")
#   testthat::expect_equal(length(rand.emax), 4)
#
#   # Class models
#   ref.resp <- classdata[classdata$treatname=="Placebo",]
#   testthat::expect_warning(ref.synth(ref.resp, exp.class.fixed, synth="random"), "Synthesis of reference")
# })
#
#
#
#
# testthat::test_that("rankauc functions correctly", {
#   model.list <- list(emax1, emax2, quad)
#   treats.list <- list(c(1,5,8,15),
#                       c(2,6,9,16:20),
#                       c(1:4))
#   int.list <- list(c(0,10), c(5,10), c(1,3))
#   subs.list <- list(10, 5, 40)
#   dec.list <- list(TRUE, FALSE, TRUE)
#   for (i in seq_along(model.list)) {
#     auc <- MBNMAtime:::rankauc(model.list[[i]], decreasing=dec.list[[i]], treats=model.list[[i]]$treatments[treats.list[[i]]],
#              int.range=int.list[[i]], subdivisions=subs.list[[i]], n.iter=100)
#
#     testthat::expect_equal(names(auc), c("summary", "prob.matrix", "rank.matrix", "auc.int"))
#     testthat::expect_equal(nrow(auc[["summary"]]), length(treats.list[[i]]))
#     testthat::expect_equal(nrow(auc[["prob.matrix"]]), ncol(auc[["prob.matrix"]]))
#     testthat::expect_equal(nrow(auc[["prob.matrix"]]), length(treats.list[[i]]))
#     testthat::expect_equal(nrow(auc[["rank.matrix"]]), 100)
#     testthat::expect_equal(colnames(auc[["rank.matrix"]]), model.list[[i]]$treatments[treats.list[[i]]])
#   }
#
#   i <- 1
#   testthat::expect_error(MBNMAtime:::rankauc(model.list[[i]], decreasing=5, treats=treats.list[[i]],
#                     int.range=int.list[[i]], subdivisions=subs.list[[i]], n.iter=100))
#
#   testthat::expect_error(MBNMAtime:::rankauc(model.list[[i]], decreasing=dec.list[[i]], treats=c("Placecbo", "Celebrex"),
#                         int.range=int.list[[i]], subdivisions=subs.list[[i]], n.iter=100))
#
#   testthat::expect_error(MBNMAtime:::rankauc(model.list[[i]], decreasing=dec.list[[i]], treats=treats.list[[i]],
#                         int.range=c(1:10), subdivisions=subs.list[[i]], n.iter=100))
#
#   i <- 2
#   testthat::expect_error(MBNMAtime:::rankauc(model.list[[i]], decreasing=dec.list[[i]], treats=treats.list[[i]],
#                         int.range=c(-5,5), subdivisions=subs.list[[i]], n.iter=100))
#
#   testthat::expect_error(MBNMAtime:::rankauc(model.list[[i]], decreasing=dec.list[[i]], treats=treats.list[[i]],
#                         subdivisions=subs.list[[i]], n.iter=100))
#
#   testthat::expect_error(MBNMAtime:::rankauc(model.list[[i]], decreasing=dec.list[[i]], treats=treats.list[[i]],
#                         int.range=int.list[[i]], subdivisions=-10, n.iter=100))
#
#   # Error due to wrong parameters being saved from model
#   testthat::expect_error(MBNMAtime:::rankauc(resdev, decreasing=dec.list[[i]],
#                         treats=treats.list[[i]],
#                         int.range=int.list[[i]],
#                         subdivisions=subs.list[[i]], n.iter=100))
#
# })
#
#
#
#
#
# testthat::test_that("rank.mbnma functions correctly", {
#   model.list <- list(emax1, emax2, quad)
#   treats.list <- list(c(1,5,8,15), c("Pl_0", "Ce_100", "Ce_200", "Va_20"))
#   i <- 1
#
#   rank <- rank(exponential, params="d.lambda",
#                      direction=-1, treats=treats.list[[i]])
#
#   testthat::expect_equal(names(rank), "d.lambda")
#   testthat::expect_equal(names(rank$d.lambda), c("summary", "prob.matrix", "rank.matrix"))
#   testthat::expect_equal(nrow(rank$d.lambda[["summary"]]), length(treats.list[[i]]))
#   testthat::expect_equal(nrow(rank$d.lambda[["prob.matrix"]]), ncol(rank$d.lambda[["prob.matrix"]]))
#   testthat::expect_equal(nrow(rank$d.lambda[["prob.matrix"]]), length(treats.list[[i]]))
#   testthat::expect_equal(nrow(rank$d.lambda[["rank.matrix"]]), model.list[[i]]$BUGSoutput$n.sims)
#
#   # Check that treatment codes can be character or numeric when estimating AUC
#   expect_error(rank(emax1, params="auc",
#                direction=-1, treats=treats.list[[1]], n.iter=100), NA)
#   expect_error(rank(emax1, params="auc",
#                     direction=-1, treats=treats.list[[2]], n.iter=100), NA)
#
#   if (is.numeric(treats.list[[i]])) {
#     matchtreat <- exponential$treatments[treats.list[[i]]]
#   } else if (is.character(treats.list[[i]])) {
#     matchtreat <- treats.list[[i]]
#   }
#   testthat::expect_equal(colnames(rank$d.lambda[["rank.matrix"]]), matchtreat)
#
#
#
#   rank <- rank(emax1, params=c("d.emax", "d.et50"),
#                      direction=-1, treats=treats.list[[i]])
#
#   testthat::expect_equal(sort(names(rank)), sort(c("d.emax", "d.et50")))
#   testthat::expect_equal(names(rank$d.et50), c("summary", "prob.matrix", "rank.matrix"))
#   testthat::expect_equal(nrow(rank$d.et50[["summary"]]), length(treats.list[[i]]))
#   testthat::expect_equal(nrow(rank$d.et50[["prob.matrix"]]), ncol(rank$d.et50[["prob.matrix"]]))
#   testthat::expect_equal(nrow(rank$d.et50[["prob.matrix"]]), length(treats.list[[i]]))
#   testthat::expect_equal(nrow(rank$d.et50[["rank.matrix"]]), model.list[[i]]$BUGSoutput$n.sims)
#
#   i <- 2
#   rank <- rank(quad, params=c("beta.2", "d.1"),
#                      direction=-1, treats=treats.list[[i]])
#
#   testthat::expect_equal(sort(names(rank)), sort(c("beta.2", "d.1")))
#   testthat::expect_equal(names(rank$d.1), c("summary", "prob.matrix", "rank.matrix"))
#   testthat::expect_equal(nrow(rank$d.1[["summary"]]), length(treats.list[[i]]))
#   testthat::expect_equal(nrow(rank$d.1[["prob.matrix"]]), ncol(rank$d.1[["prob.matrix"]]))
#   testthat::expect_equal(nrow(rank$beta.2[["prob.matrix"]]), length(treats.list[[i]]))
#   testthat::expect_equal(nrow(rank$beta.2[["rank.matrix"]]), model.list[[i]]$BUGSoutput$n.sims)
#
#   if (is.numeric(treats.list[[i]])) {
#     matchtreat <- quad$treatments[treats.list[[i]]]
#   } else if (is.character(treats.list[[i]])) {
#     matchtreat <- treats.list[[i]]
#   }
#   testthat::expect_equal(colnames(rank$beta.2[["rank.matrix"]]), matchtreat)
#   testthat::expect_equal(colnames(rank$d.1[["rank.matrix"]]), matchtreat)
#
#
#   testthat::expect_error(rank(emax1, params=c("beta.1", "beta.2"),
#                           direction=-1, treats=treats.list[[i]]))
#
#   # Class effect models
#   testthat::expect_error(rank(exp.class.fixed,
#                               direction=-1, param="D.lambda", treats=c("1","wer"), level="class"), "classes not included")
#
#   testthat::expect_error(rank(exp.class.fixed,
#                               direction=-1, treats=c("1","2"), param="auc", level="class"), "AUC cannot currently")
#
#   testthat::expect_silent(rank(exp.class.fixed,
#                               direction=-1, treats=c(1,2), param="D.lambda", level="class"))
#
#   expect_silent(rank(exp.class.fixed,
#                      direction=-1, param="D.lambda", level="class"))
#
#   expect_error(rank(emax1, direction=-1, level="class"), "not a class effect model")
# })
