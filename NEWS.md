# MBNMAtime 0.2.3

## Additions/changes
- Truncated priors are used as the default priors for time-course parameters that can only take positive values (e.g. `et50`) for all functions
- `texp()` has been removed - `titp()` is a more stable parameterisation of this function

## Bug fixes

### Major
- Error fixed preventing natural splines from functioning, causing removal from CRAN

### Minor
- Return values added to documentation for functions in which they were missing

# MBNMAtime 0.2.2

## Additions/changes
- Can now specify numeric values for time-course parameters in the `method` argument. Can be useful for discrete values that cannot be estimated (e.g. fractional polynomial powers, Hill parameter).
- Fractional polynomial powers in `tfpoly()` can only take numeric values from set defined in Jansen 2015.
- Integrated Two-Component Prediction (ITP) function (`titp()`) added
- `get.relative()` can be used to combine two MBNMA models to allow different time-course functions to be fitted to a different set of treatments (see examples in the vignette)
- New priors that restrict posterior to positive values where necessary can be easily incorporated. 
- `binplot()` can be used to plot the results of NMAs conducted at multiple time bins. This can be particularly useful to explore which time-course functions might be appropriate, and to check the validity of MBNMA predictions.
- `mb.nodesplit()` can be performed at specific time-points, in addition to by time-course parameter
- `corparam` set to `FALSE` as default

## Bug fixes

### Minor
- Error with `overlay.nma` argument in `plot.mb.predict()` fixed


# MBNMAtime 0.2.1

## Additions/changes
- `get.relative()` function can be used to calculate relative effects/mean differences between treatments/classes
- `cumrank()` added for cumulative ranking plots. Also calculates SUCRA values for each treatment and time-course parameter
at specified follow-up times (even those at which treatments have not been compared within any study)
- Studies reporting change from baseline or absolute means can now be specified in `mb.network()`, or
will be automatically inferred from the data (studies with no time=0 are assumed to report change
from baseline)
- Model intercept (response at time=0) is now conditional on change from baseline *for each study*
- `texp()` now implements 2-parameter exponential function (though the simpler 1-parameter model remains the default)

## Bug fixes

### Major
- Error with `predict()` not properly incorporating absolute time-course parameters fixed

### Minor
- Error with `model.file` input length fixed for `mb.run()`


# MBNMAtime 0.2.0

## Additions/changes
- Added variance adjustment (`covar="varadj"`) for correlation between time-points - this is now the default in `mb.run()`
- Added log linear time-course function (`tloglin()`)
- Added spline functions (piecewise linear splines, B-splines, restricted cubic splines, natural splines)
- Added `overlay.nma` option to `predict()` to allow plotting of "lumped" NMA results over MBNMA predictions
- Modelling can now incorporate Standardised Mean Differences (`link="smd"`) or Ratios of Means (`link="log"`) to allow modelling of studies with different scales
- `lower_better` argument used instead of `decreasing` for rankings
- Time-course functions given to `mb.run()` are now given as `class("timefun")` and time-course parameters are specified within these functions
- Predictions from `predict()` can now be ranked
- Forest plots now also plot posterior densities using `ggdist::stat_halfeye()`
- Neater outputs when using `print()` or `summary()`
- Wishart prior used to model correlations between within-study baseline effects on different time-course parameters, in addition
to relative effects on different time-course parameters.


## Bug fixes

### Major
- Corrected calculation for Bayesian p-value in `mb.nodesplit()`


# MBNMAtime 0.1.3

## Additions/changes
- Added citation file
- `plot.mb.network()` now uses a `layout` argument that takes an igraph layout function instead of `layout_in_circle` (which was a logical argument). This allows any igraph layout to be plotted rather than just a circle (e.g. `igraph::as_star()`)
- Objects returned from `plot.mb.network` now have specific igraph attributes assigned to them, which can be easily changed by the user.
- `user.fun` now takes a formula as an argument (for example `~ (beta.1 * dose) + (beta.2 * dose^2)`) rather than a string.
- `mb.network` objects are now stored within lists of most other mb class objects for easy reference of data format

## Bug fixes

### Major
- Exponential function models were not working previously but the dose-response function has been rewritten so that it runs the model correctly.
- Ensured comparisons are cycled through correctly in mb.nodesplit
- Ensured `timeplot` raw responses can be plotted by either arm (`plotby="arm"`) or relative (`plotby="rel"`) effects.

# MBNMAtime 0.1.2

## First release of package

Welcome to MBNMAtime. Ready for release into the world. I hope it can be of service to you! For dose-response MBNMA, also check out the sister package, MBNMAdose.
