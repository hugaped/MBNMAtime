# MBNMAtime 0.1.3

## Additions/changes
- Added citation file
- `plot.mb.network()` now uses a `layout` argument that takes an igraph layout function instead of `layout_in_circle` (which was a logical argument). This allows any igraph layout to be plotted rather than just a circle (e.g. `igraph::as_star()`)
- Objects returned from `plot.mb.network` now have specific igraph attributes assigned to them, which can be easily changed by the user.

## Bug fixes

### Major
- Exponential function models were not working previously but the dose-response function has been rewritten so that it runs the model correctly.
- Ensured comparisons are cycled through correctly in mb.nodesplit


# MBNMAtime 0.1.2

## First release of package

Welcome to MBNMAtime. Ready for release into the world. I hope it can be of service to you! For dose-response MBNMA, also check out the sister package, MBNMAdose.
