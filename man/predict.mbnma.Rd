% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbnma-class.R
\name{predict.mbnma}
\alias{predict.mbnma}
\title{Predict effects over time in a given population based on MBNMA time-course
models}
\usage{
\method{predict}{mbnma}(
  object,
  times = seq(0, max(object$model.arg$jagsdata$time, na.rm = TRUE), length.out = 30),
  E0 = 0,
  treats = NULL,
  level = "treatment",
  ref.resp = NULL,
  synth = "common",
  lim = "cred",
  ...
)
}
\arguments{
\item{object}{An S3 object of \code{class("mbnma")} generated by running
a time-course MBNMA model}

\item{times}{A sequence of positive numbers indicating which time points to
predict mean responses for (or at which to conduct a node-split if used with \code{mb.nodesplit()})}

\item{E0}{An object to indicate the value(s) to use for the response at time = 0
in the prediction. This can take a number of different formats depending
on how it will be used/calculated. The default is 0 but this may lead
to non-sensical predictions if Ratio of Means are modeled.
\itemize{
\item \code{numeric()} A single numeric value representing the deterministic response at time = 0
\item \code{formula()} A formula representing a stochastic distribution for the response
at time = 0. This is specified as a random number generator
(RNG) given as a string, and can take any RNG distribution for which a function exists
in R. For example: \code{~rnorm(n, 7, 0.5)}.
}}

\item{treats}{A character vector of treatment/class names or a numeric vector of treatment/class codes (as coded
in \code{mbnma}) that indicates which treatments/classes to calculate predictions for. If left as \code{NULL} then
predictions will be calculated for all treatments/classes. Whether the vector should correspond to treatments or
classes depends on the value of \code{level}.}

\item{level}{Can take either \code{"treatment"} to make predictions for treatments, or \code{"class"} to make predictions for classes (in
which case \code{object} must be a class effect model).}

\item{ref.resp}{An object to indicate the value(s) to use for the reference treatment response in MBNMA models
in which the reference treatment response is not estimated within the model (i.e. those that model any time-
course parameters using \code{pool="rel"}). This can take a number of different formats depending
on how it will be used/calculated. There are two approaches for this:
\enumerate{
\item The reference response can be estimated from a dataset of studies investigating the reference
treatment using meta-analysis. This dataset could be a set of observational
studies that are specific to the population on which to make
predictions, or it could be a subset of the study arms within the MBNMA dataset
that investigate the reference treatment. The data should be provided to \code{ref.resp} as a
\code{data.frame()} containing the data in long format (one row per observation). See \code{\link[=ref.synth]{ref.synth()}}
\item Values for the reference treatment response can be assigned to different time-course parameters
within the model that have been modelled using consistency relative effects (\code{pool="rel"}).
These are given as a list, in which each named element corresponds to a time-course
parameter modelled in \code{mbnma}, specified on the corresponding scale (i.e. specified on the log scale if modelled
on the log scale using Ratios of Means). Their values can be either of the following:
}
\itemize{
\item \code{numeric()} A numeric value representing the deterministic value of the time-course parameter in
question in individuals given the reference treatment. \code{0} is used as the default, which assumes no
effect of time on the reference treatment (i.e. mean differences / relative effects versus the
reference treatment are modeled).
\item \code{formula()} A formula representing a stochastic distribution for the value of the time-course
parameter in question. This is specified as a random number generator (RNG) given as a formula,
and can take any RNG distribution for which a function exists in R. For example: \code{~rnorm(n, -3, 0.2)}.
}}

\item{synth}{A character object that can take the value \code{"common"} or \code{"random"} that
specifies the the type of pooling to use for synthesis of \code{ref.resp}. Using \code{"random"} rather
than \code{"common"} for \code{synth} will result in wider 95\\\% CrI for predictions.}

\item{lim}{Specifies calculation of either 95\% credible intervals (\code{lim="cred"}) or 95\% prediction intervals (\code{lim="pred"}).}

\item{...}{Arguments to be sent to R2jags for synthesis of the network
reference treatment effect (using \code{\link[=ref.synth]{ref.synth()}})}
}
\value{
An S3 object of class \code{mb.predict} that contains the following
elements:
\itemize{
\item \code{summary} A named list of data frames. Each data frame contains
a summary of predicted responses at follow-up times specified in \code{times}
for each treatment specified in \code{treats}
\item \code{pred.mat} A named list of
matrices. Each matrix contains the MCMC results of predicted responses at
follow-up times specified in \code{times} for each treatment specified in
\code{treats}
}
}
\description{
Used to predict effects over time for different treatments or to predict
the results of a new study. For MBNMA models that include consistency
relative effects on time-course parameters, this is calculated by combining
relative treatment effects with a given reference treatment response
(specific to the population of interest).
}
\details{
By default the network reference treatment baseline (\code{E0}) and time-course
parameter values are set to zero so that \code{predict()} estimates mean differences
(/relative treatment effects) over time versus the network reference treatment.

\code{ref.resp} only needs to be specified if \code{mbnma} has
been estimated using consistency relative effects (\code{pool="rel"}) for
any time-course parameters, as these inform the absolute values of the
network reference treatment parameters which can then be added to the
relative effects to calculate specific predictions.
}
\examples{
\donttest{
# Create an mb.network object from a dataset
network <- mb.network(osteopain)

# Run an MBNMA model with an Emax time-course
emax <- mb.run(network,
  fun=temax(pool.emax="rel", method.emax="common",
    pool.et50="abs", method.et50="common"))

# Predict responses using a stochastic baseline (E0) and a distribution for the
#network reference treatment
preds <- predict(emax, times=c(0:10),
  E0=~rnorm(n, 7, 0.5),
  ref.resp=list(emax=~rnorm(n, -0.5, 0.05)))
summary(preds)

# Predict responses using the original dataset to estimate the network reference
#treatment response
paindata.ref <- osteopain[osteopain$treatname=="Placebo_0",]
preds <- predict(emax, times=c(5:15),
  E0=10,
  ref.resp=paindata.ref)
summary(preds)

# Repeat the above prediction but using a random effects meta-analysis of the
#network reference treatment response
preds <- predict(emax, times=c(5:15),
  E0=10,
  ref.resp=paindata.ref,
  synth="random")
summary(preds)
}

}
