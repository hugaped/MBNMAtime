#' MBNMAtime for Model-Based Network Meta-Analysis of longitudinal (time-course) data
#'
#' @description
#' MBNMAtime provides a collection of useful commands that allow users to run time-course
#' Model-Based Network Meta-Analyses (MBNMA) or Model-Based Meta-Analyses (MBMA).
#'
#' @section Introduction:
#' MBNMAtime allows meta-analysis of studies with multiple follow-up measurements that can
#' account for time-course for a single or multiple treatment comparisons.
#'
#' Including all available follow-up measurements within a study makes use of all the available
#' evidence in a way that maintains connectivity between treatments, and it does so in a way
#' that explains time-course, thus explaining heterogeneity and inconsistency that may be
#' present in a standard Network Meta-Analysis (NMA). All models and analyses are implemented
#' in a Baysian framework, following an extension of the standrd NMA methodology presented by
#' \insertCite{lu2004}{MBNMAtime} and are run in JAGS \insertCite{jags}{MBNMAtime}. For full details of time-course MBNMA
#' methodology see \insertCite{pedder2019;textual}{MBNMAtime}.
#'
#' @section Workflow:
#' Functions within `MBNMAtime` follow a clear pattern of use:
#'
#' 1. Load your data into the correct format using \code{\link{mb.network}}
#' 2. Analyse your data using \code{\link{mb.run}}, or any of the available wrapper time-course functions
#' 3. Test for consistency using functions like \code{\link{mb.nodesplit}}
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
#' network <- mb.network(osteopain, ref="Pl_0")
#'
#' # Generate a network plot
#' plot(network, label.distance=3)
#'
#' # Analyse data using mb.run()
#' result <- mb.run(network, fun="emax",
#'   beta.1=list(pool="rel", method="common"),
#'   beta.2=list(pool="arm", method="common"),
#'   positive.scale=TRUE)
#'
#' # ...or achieve the same result by using a wrapper function for mb.run()
#' result <- mb.emax(network,
#'   emax=list(pool="rel", method="common"),
#'   et50=list(pool="arm", method="common"),
#'   positive.scale=TRUE)
#'
#' # Explore model fit statistics - plot residual deviances
#' devplot(result, n.iter=500)
#'
#' # Generate a forest plot for model results
#' plot(result)
#'
#' # Predict responses
#' pred <- predict(result, time=c(0:10), baseline=10,
#'   ref.data=list("emax"="rnorm(nsims,-2,0.5)"),
#'   treats=c("Pl_0", "Ce_400", "Et_5", "Ox_44", "Tr_300"))
#'
#' # Plot predicted response
#' plot(pred, disp.obs=TRUE)
#'
#' # Rank by Area Under the time-course Curve
#' ranks <- rank(result, param="auc", direction=-1, n.iter=500)
#'
#' # Plot histogram of rankings
#' plot(ranks)
#' }
#'
#' @keywords internal
"_PACKAGE"
