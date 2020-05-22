# MBNMAtime 0.1.4

## Additions/changes
- Corrected calculation for Bayesian p-value in `mb.nodesplit()`

## Bug fixes

### Major


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
