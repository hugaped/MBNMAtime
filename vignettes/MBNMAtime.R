## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE
)

library(MBNMAtime)
library(rmarkdown)
library(knitr)
library(dplyr)

## ---- echo=FALSE---------------------------------------------------------
kable(head(osteopain), digits=2) 

## ---- echo=FALSE---------------------------------------------------------
kable(head(alog_pcfb), digits=2) 

## ---- echo=FALSE---------------------------------------------------------
kable(head(obesityBW_CFB), digits=2) 

## ---- echo=FALSE---------------------------------------------------------
kable(head(goutSUA_CFB), digits=2) 

## ----network.pain--------------------------------------------------------
# Using the pain dataset
network.pain <- mb.network(osteopain, reference = "Pl_0")
print(network.pain)

## ------------------------------------------------------------------------
# Prepare data using the alogliptin dataset
network.alog <- mb.network(alog_pcfb, reference = "placebo")

# Plot network
plot(network.alog)

## ----gout.class.network--------------------------------------------------
plot(network.gout, level = "class", remove.loops = TRUE, label.distance = 5)

## ---- warning=FALSE------------------------------------------------------
plot(network.gout, level = "treatment", v.color = "class", label.distance = 5)

## ------------------------------------------------------------------------
network.gout <- mb.network(goutSUA_CFBcomb, reference="Plac")
plot(network.gout, label.distance = 5)

## ----pain.time-----------------------------------------------------------
# Prepare data using the pain dataset
network.pain <- mb.network(osteopain, reference="Pl_0")

# Draw plot of raw study responses over time
timeplot(network.pain)

## ----obese.time, message=FALSE-------------------------------------------
# Draw plot of raw study responses over time grouped by agent class in the obesity dataset
network.obese <- mb.network(obesityBW_CFB)
timeplot(network.obese, level="class")

## ---- results="hide"-----------------------------------------------------
# Run a linear time-course MBNMA
mbnma <- mb.run(network.pain, fun="linear", 
                   beta.1=list(pool="rel", method="common"))

## ------------------------------------------------------------------------
summary(mbnma)

## ---- eval=FALSE---------------------------------------------------------
#  # An alternative would be to use a linear wrapper for mb.run() which would give the same result
#  mb.linear(network.pain,
#                    slope=list(pool="rel", method="common"))

## ---- results="hide", message=FALSE--------------------------------------
# Run an emax time-course MBNMA pooling absolute effects
mbnma <- mb.run(network.gout, fun="emax", 
                   beta.1=list(pool="arm", method="random"), 
                   beta.2=list(pool="const", method="common"), 
                   intercept=FALSE)

## ------------------------------------------------------------------------
summary(mbnma)

## ---- eval=FALSE---------------------------------------------------------
#  # An alternative would be to use an emax wrapper for mb.run() which would give the same result
#  mb.emax(network.gout,
#             emax=list(pool="arm", method="random"),
#             et50=list(pool="const", method="common"),
#             intercept=FALSE)

## ---- results="hide"-----------------------------------------------------
# Run a piecewise linear time-course MBNMA
mbnma <- mb.run(network.pain, fun="piecelinear", 
                   beta.1=list(pool="rel", method="common"), 
                   beta.2=list(pool="arm", method="common"), 
                   beta.3=list(pool="const", method=1))

## ------------------------------------------------------------------------
summary(mbnma)

## ---- eval=FALSE---------------------------------------------------------
#  # An alternative would be to use a piecewise linear wrapper for mb.run() which would give the same result
#  mb.piecelinear(network.pain,
#                    slope.1=list(pool="rel", method="common"),
#                    slope.2=list(pool="arm", method="common"),
#                    knot=list(pool="const", method=1),
#                    n.iter=5000)

## ---- eval=FALSE---------------------------------------------------------
#  # Run an emax time-course MBNMA that accounts for correlation between time points using a wrapper for mb.run()
#  mbnma <- mb.emax(network.pain,
#                      emax=list(pool="rel", method="common"),
#                      et50=list(pool="const", method="common"),
#                      rho="estimate", covar="CS")

## ---- results="hide", warning=FALSE--------------------------------------
# Run an emax time-course MBNMA with a random class effects on beta.1 (Emax parameters)
# Additional iterations run to ensure MCMC chain convergence
mbnma <- mb.run(network.gout, fun="emax", 
                   beta.1=list(pool="rel", method="random"), 
                   beta.2=list(pool="const", method="common"), 
                   intercept=FALSE, n.iter=20000,
                   class.effect=list(beta.1="random"))

## ------------------------------------------------------------------------
summary(mbnma)

## ---- results="hide", error=TRUE-----------------------------------------
# Run an MBNMA that accounts for correlation between time points using AR1 covariance structure
mbnma <- mb.emax(network.pain, 
                    emax=list(pool="rel", method="common"), 
                    et50=list(pool="const", method="common"), 
                    rho="estimate", covar="AR1")

print(mbnma$model.arg$priors)

## ---- eval=FALSE---------------------------------------------------------
#  # Define prior for rho which permits evaluation
#  new.priors <- list(
#    "rho" = "dunif(0, 1)"
#    )
#  
#  # Run an MBNMA model with new priors
#  mbnma <- mb.emax(network.pain,
#                      emax=list(pool="rel", method="common"),
#                      et50=list(pool="const", method="common"),
#                      rho="estimate", covar="AR1",
#                      priors=new.priors)

## ---- eval=FALSE---------------------------------------------------------
#  # Define relative magnitudes of slope.1 and slope.2
#  rel.size <- c(10,1)
#  
#  mbnma <- mb.piecelinear(network.pain,
#                          slope.1=list(pool="rel", method="random"),
#                          slope.2=list(pool="rel", method="common"),
#                          knot=list(pool="const", method=1),
#                          var.scale = rel.size)

## ---- results="hide", fig.show="hold"------------------------------------
# Run a first-order fractional polynomial time-course MBNMA
mbnma <- mb.fract.first(network.pain, 
                           slope=list(pool="rel", method="random"),
                           power=list(pool="const", method="common")
                           )

# Plot a scatter plot of residual deviance contributions (the default)
devplot(mbnma, n.iter=1000)

## ---- results="hide"-----------------------------------------------------
# Plot boxplots of deviance contributions
devplot(mbnma, dev.type = "dev", plot.type = "box", n.iter=1000)

## ---- results="hide"-----------------------------------------------------
# Plot fitted and observed values with treatment labels
fitplot(mbnma, n.iter=1000)

## ---- results="hide"-----------------------------------------------------
# Run a quadratic time-course MBNMA using the alogliptin dataset
mbnma <- mb.quadratic(network.alog, 
                    beta.1=list(pool="rel", method="random"), 
                    beta.2=list(pool="rel", method="common"),
                    intercept=FALSE)

plot(mbnma)

## ---- eval=FALSE, results="hide"-----------------------------------------
#  # Run a piecewise linear time-course MBNMA with a knot at 1 week
#  mbnma <- mb.piecelinear(network.pain,
#                             slope.1=list(pool="rel", method="common"),
#                             slope.2=list(pool="rel", method="common"),
#                             knot=list(pool="const", method=1))
#  
#  # Rank results based on AUC (calculated 0-10 weeks), more negative slopes considered to be "better"
#  ranks <- rank(mbnma, params=c("auc", "d.slope.2"),
#                      int.range=c(0,10),  direction=-1, n.iter=1000)

## ---- echo=FALSE---------------------------------------------------------
load("vignettedata/ranks.RData")

## ------------------------------------------------------------------------
print(ranks)

## ------------------------------------------------------------------------
# Ranking histograms for AUC
plot(ranks, params = "auc")

## ---- results="hide", message=FALSE--------------------------------------
mbnma <- mb.emax(network.pain, 
                    emax=list(pool="rel", method="common"), 
                    et50=list(pool="const", method="common"))

# Generate a dataset that is made up only of network reference treatment responses over time (in this case placebo)
placebo.data <- network.pain$data.ab[network.pain$data.ab$treatname=="Placebo_0",]

# Predict responses for a selection of treatments using a deterministic E0 and placebo.data to estimate the network reference treatment effect
predict.data <- predict(mbnma, treats=c("Pl_0", "Ce_200", "Du_90", "Et_60", 
                                        "Lu_400", "Na_1000", "Ox_44", "Ro_25",
                                        "Tr_300", "Va_20"),
                        times=c(0:15), E0=10, 
                        ref.resp=placebo.data)

## ------------------------------------------------------------------------
# Summary of posterior median predictions
summary(predict.data)

## ---- results="hide", eval=TRUE, message=FALSE---------------------------
# Define stochastic values for network reference treatment effect on Emax
ref.data <- list("emax"="rnorm(n, -2, 0.15)")

# Predict responses for treatments 1-9 using a stochastic E0 and ref.resp to estimate the network reference treatment effect
predict.resp <- predict(mbnma, treats=c("Pl_0", "Ce_200", "Du_90", "Et_60", 
                                        "Lu_400", "Na_1000", "Ox_44", "Ro_25",
                                        "Tr_300", "Va_20"),
                        times=c(0:15), E0="rnorm(n, 9, 0.05)",
                        ref.resp=ref.data)

## ------------------------------------------------------------------------
plot(predict.resp, overlay.ref=TRUE, disp.obs=TRUE)

## ---- fig.height=3, results="hide"---------------------------------------
# Fit a quadratic time-course MBNMA to the Obesity dataset
mbnma <- mb.quadratic(network.obese, 
                         beta.1 = list(pool="rel", method="common"),
                         beta.2 = list(pool="rel", method="common")
                         )

# Define stochastic values centred at zero for network reference treatment
ref.data <- list(beta.1="rnorm(n, 0, 0.05)", beta.2="rnorm(n, 0, 0.0001)")

# Predict responses over the
predict.resp <- predict(mbnma, times=c(0:50), E0="rnorm(n, 120,1)", treats = c(1,4,15),
                        ref.resp=ref.data)

# Plot predictions
plot(predict.resp, disp.obs = TRUE)

## ---- warning=FALSE------------------------------------------------------
# Loops of evidence within the alogliptin dataset
splits.alog <- mb.nodesplit.comparisons(network.alog)
print(splits.alog)

## ---- eval=FALSE---------------------------------------------------------
#  # Fit a piecewise linear MBNMA with fixed relative effects on slope.1 and slope.2
#  mbnma <- mb.piecelinear(network.pain,
#                             slope.1=list(pool="rel", method="common"),
#                             slope.2=list(pool="rel", method="common"),
#                             knot=list(pool="const", method=0.5),
#                             pd="pd.kl")
#  
#  # Fit a UME model on both slope parameters simultaneously in a piecewise linear MBNMA
#  ume <- mb.piecelinear(network.pain,
#                           slope.1=list(pool="rel", method="common"),
#                           slope.2=list(pool="rel", method="common"),
#                           knot=list(pool="const", method=0.5),
#                           UME=TRUE, pd="pd.kl")
#  
#  # Fit a UME model on slope.1 only in a piecewise linear MBNMA
#  ume.slope.1 <- mb.piecelinear(network.pain,
#                           slope.1=list(pool="rel", method="common"),
#                           slope.2=list(pool="rel", method="common"),
#                           knot=list(pool="const", method=0.5),
#                           UME="slope.1", pd="pd.kl")
#  
#  # Fit a UME model on slope.2 only in a piecewise linear MBNMA
#  ume.slope.2 <- mb.piecelinear(network.pain,
#                           slope.1=list(pool="rel", method="common"),
#                           slope.2=list(pool="rel", method="common"),
#                           knot=list(pool="const", method=0.5),
#                           UME="slope.2", pd="pd.kl")
#  

## ---- echo=FALSE, eval=FALSE---------------------------------------------
#  umelist.1 <- list(paste("Deviance for mbnma:",
#                          round(mbnma$BUGSoutput$median$deviance,2), sep=" "),
#                  paste("Deviance for ume on slope.1 and slope.2:",
#                        round(ume$BUGSoutput$median$deviance,2), sep=" "),
#                  paste("Deviance for ume on slope.1:",
#                        round(ume.slope.1$BUGSoutput$median$deviance,2), sep=" "),
#                  paste("Deviance for ume on slope.2:",
#                        round(ume.slope.2$BUGSoutput$median$deviance,2), sep=" ")
#                  )
#  save(umelist.1, file="umelist.1.Rdata")

## ---- echo=FALSE---------------------------------------------------------
load("vignettedata/umelist.1.RData")
for (i in seq_along(umelist.1)) {
  print(umelist.1[[i]])
}

## ---- eval=FALSE---------------------------------------------------------
#  # Run an Emax MBNMA with random relative effects on emax
#  mbnma <- mb.emax(network.gout,
#                      emax=list(pool="rel", method="random"),
#                      et50=list(pool="const", method="common"),
#                      n.iter=10000, n.thin=10,
#                      intercept=FALSE)
#  
#  # Fit a UME model on Emax parameters
#  ume <- mb.emax(network.gout,
#                    emax=list(pool="rel", method="random"),
#                    et50=list(pool="const", method="common"),
#                    n.iter=10000, n.thin=10,
#                    intercept=FALSE, UME=TRUE)

## ---- echo=FALSE, eval=FALSE---------------------------------------------
#  umelist.2 <- list(paste("Deviance for mbnma:", round(mbnma$BUGSoutput$median$deviance,2), sep=" "),
#                  paste("Deviance for UME:", round(ume$BUGSoutput$median$deviance,2), sep=" ")
#                  )
#  umelist.sd <- list(paste("SD for mbnma:", round(mbnma$BUGSoutput$median$sd.emax,2), sep=" "),
#                     paste("SD for UME:", round(ume$BUGSoutput$median$sd.emax,2), sep=" "))
#  save(umelist.2, umelist.sd, file="umelist.2.RData")

## ---- echo=FALSE---------------------------------------------------------
load("vignettedata/umelist.2.RData")
for (i in seq_along(umelist.2)) {
  print(umelist.2[[i]])
}

## ---- echo=FALSE---------------------------------------------------------
# Compare between-study SD between models
for (i in seq_along(umelist.sd)) {
  print(umelist.sd[[i]])
}

## ---- eval=FALSE---------------------------------------------------------
#  # Nodesplit using an Emax MBNMA
#  nodesplit <- mb.nodesplit(network.pain, fun="emax",
#                               beta.1=list(pool="rel", method="random"),
#                               beta.2=list(pool="const", method="common"),
#                               nodesplit.parameters="all")
#  

## ---- echo=FALSE, eval=FALSE---------------------------------------------
#  for (i in seq_along(nodesplit)) {
#    for (k in seq_along(nodesplit[[i]])) {
#      nodesplit[[i]][[k]]$direct <- NULL
#      nodesplit[[i]][[k]]$indirect <- NULL
#  
#      ggdat <- nodesplit[[i]][[k]]$forest.plot$data
#      nodesplit[[i]][[k]]$forest.plot <- NULL
#      nodesplit[[i]][[k]]$forest.plot <- list(data=ggdat)
#  
#      ggdat <- nodesplit[[i]][[k]]$density.plot$data
#      nodesplit[[i]][[k]]$density.plot <- NULL
#      nodesplit[[i]][[k]]$density.plot <- list(data=ggdat)
#    }
#  }
#  save(nodesplit, file="nodesplit.1.RData")

## ---- echo=FALSE---------------------------------------------------------
load("vignettedata/nodesplit.1.RData")

## ------------------------------------------------------------------------
print(nodesplit)

## ---- fig.height=2.5, fig.show="hold"------------------------------------
# Plot forest plots of direct and indirect results for each node-split comparison
plot(nodesplit, plot.type="forest")

# Plot posterior densities of direct and indirect results for each node-split comparisons
plot(nodesplit, plot.type="density")

## ---- eval=FALSE---------------------------------------------------------
#  # Nodesplit on beta.1 and beta.2 using a piecewise linear MBNM
#  nodesplit <- mb.nodesplit(network.pain, fun="piecelinear",
#                               beta.1=list(pool="rel", method="common"),
#                               beta.2=list(pool="rel", method="common"),
#                               beta.3=list(pool="const", method=0.5),
#                               nodesplit.parameters="all")

## ---- echo=FALSE, eval=FALSE---------------------------------------------
#  for (i in seq_along(nodesplit)) {
#    for (k in seq_along(nodesplit[[i]])) {
#      nodesplit[[i]][[k]]$direct <- NULL
#      nodesplit[[i]][[k]]$indirect <- NULL
#  
#      ggdat <- nodesplit[[i]][[k]]$forest.plot$data
#      nodesplit[[i]][[k]]$forest.plot <- NULL
#      nodesplit[[i]][[k]]$forest.plot <- list(data=ggdat)
#  
#      ggdat <- nodesplit[[i]][[k]]$density.plot$data
#      nodesplit[[i]][[k]]$density.plot <- NULL
#      nodesplit[[i]][[k]]$density.plot <- list(data=ggdat)
#    }
#  }
#  save(nodesplit, file="nodesplit.2.RData")

## ---- echo=FALSE---------------------------------------------------------
load("vignettedata/nodesplit.2.RData")

## ---- fig.height=2.5-----------------------------------------------------
print(nodesplit)

plot(nodesplit, plot.type="forest")

