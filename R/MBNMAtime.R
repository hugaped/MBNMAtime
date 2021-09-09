#' MBNMAtime for Model-Based Network Meta-Analysis of longitudinal (time-course) data
#'
#' @description
#' MBNMAtime provides a collection of useful commands that allow users to run time-course
#' Model-Based Network Meta-Analysis (MBNMA).
#'
#' @section Introduction:
#' MBNMAtime allows meta-analysis of studies with multiple follow-up measurements that can
#' account for time-course for a single or multiple treatment comparisons.
#'
#' Including all available follow-up measurements within a study makes use of all the available
#' evidence in a way that maintains connectivity between treatments, and it does so in a way
#' that explains time-course, thus explaining heterogeneity and inconsistency that may be
#' present in a standard Network Meta-Analysis (NMA). All models and analyses are implemented
#' in a Bayesian framework, following an extension of the standard NMA methodology presented by
#' \insertCite{lu2004}{MBNMAtime} and are run in JAGS \insertCite{jags}{MBNMAtime}. Correlation between
#' time-points can be accounted for in the mdoelling framework. For full details of time-course MBNMA
#' methodology see \insertCite{pedder2019;textual}{MBNMAtime}.
#'
#' @section Workflow:
#' Functions within `MBNMAtime` follow a clear pattern of use:
#'
#' 1. Load your data into the correct format using \code{\link{mb.network}}
#' 2. Specify a suitable time-course function and analyse your data using \code{\link{mb.run}}
#' 3. Test for consistency using \code{\link{mb.nodesplit}} or by fitting Unrelated Mean Effects models
#' 4. Examine model results using forest plots and treatment rankings
#' 5. Use your model to predict responses using \code{\link{predict.mbnma}}
#'
#' At each of these stages there are a number of informative plots that can be generated to help make sense of your data and the models that you are fitting.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \donttest{
#' # Generate an "mb.network" object that stores data in the correct format
#' network <- mb.network(osteopain)
#'
#' # Generate a network plot
#' plot(network, label.distance=3)
#'
#' # Analyse data using mb.run()
#' result <- mb.run(network, fun=tloglin())
#'
#' # Time-course parameters can be explicitly specified
#' # Correlation between time-points can be accounted for
#' result <- mb.run(network,
#'   fun=temax(pool.emax="rel", method.emax="common",
#'     pool.et50="rel", method.et50="common"),
#'   rho="dunif(0,1)")
#'
#' # Explore model fit statistics - plot residual deviances
#' devplot(result, n.iter=500)
#'
#' # Generate a forest plot for model results
#' plot(result)
#'
#' # Predict responses
#' pred <- predict(result, time=c(0:10), E0=8,
#'   ref.resp=subset(osteopain, treatment=="Pl_0"))
#'
#' # Plot predicted response
#' plot(pred, disp.obs=TRUE)
#'
#' # Rank by Area Under the time-course Curve
#' ranks <- rank(result, param="auc", lower_better=TRUE, n.iter=500)
#'
#' # Plot histogram of rankings
#' plot(ranks)
#' }
#'
#' @keywords internal
"_PACKAGE"
