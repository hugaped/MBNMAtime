##############################################
#### Functions for class("mb.network") ####
##############################################





#' @describeIn mb.network Generate a network plot
#'
#' @param x An object of class `mb.network`.
#' @param layout An igraph layout specification. This is a function specifying an igraph
#'   layout that determines the arrangement of the vertices (nodes). The default
#'   `igraph::as_circle()` arranged vertices in a circle. Two other useful layouts for
#'   network plots are: `igraph::as_star()`, `igraph::with_fr()`. Others can be found
#'   in \code{\link[igraph]{layout_}}
#' @param edge.scale A number to scale the thickness of connecting lines
#'   (edges). Line thickness is proportional to the number of studies for a
#'   given comparison. Set to `0` to make thickness equal for all comparisons.
#' @param v.color Can take either `"connect"` (the default) to indicate that nodes should
#'   only be coloured if they are connected to the network reference treatment (indicates
#'   network connectivity) or `"class"` to colour nodes by class (this requires that the
#'   variable `class` be included in the dataset).
#' @param v.scale A number with which to scale the size of the nodes. If the variable `N`
#'   (to indicate the numbers of participants at each observation) is included in the
#'   dataset then the size of the nodes will be proportional to the number of participants
#'   within a treatment/class in the network *at the earliest time point reported in each study*.
#' @param label.distance A number scaling the distance of labels from the nodes
#'   to improve readability. The labels will be directly on top of the nodes if
#'   the default of `0` is used. Option only applicable if `layout_in_circle` is
#'   set to `TRUE`.
#' @param level A string indicating whether nodes/facets should represent `treatment`
#'   or `class` in the plot. Can be used to examine the expected impact of modelling
#'   class/agent effects.
#' @param remove.loops A boolean value indicating whether to include loops that
#'   indicate comparisons within a node.
#' @param legend A boolean value indicating whether or not to plot a legend with class names if `v.color="class"`
#' @param legend.x Can be either a string or a numerical x-coordinate indicating where the legend should be
#'   plotted (see \code{\link[graphics]{legend}}).
#' @param legend.y A numerical y-coordinate indicating where the legend should be plotted - only required if `legend.x` is also
#'   a numeric co-ordinate.
#' @param ... Options for plotting in `igraph`.
#'
#' @return Returns an object of class `"igraph"`, which can be modified by other
#' functions within the `igraph` package.
#'
#'
#' @details The S3 method `plot()` on an `mb.network` object generates a
#'   network plot that shows how different treatments are connected within the
#'   network via study comparisons. This can be used to identify how direct and
#'   indirect evidence are informing different treatment comparisons. Depends on
#'   \code{\link[igraph]{igraph}}.
#'
#' @examples
#'
#' \donttest{
#'
#' # Create an mb.network object from the data
#' network <- mb.network(osteopain)
#'
#' # Arrange network plot in a star with the reference treatment in the centre
#' plot(network, layout=igraph::as_star())
#'
#' # Generate a network plot at the class level that removes loops indicating comparisons
#' #within a node
#' goutnet <- mb.network(goutSUA_CFB)
#' plot(goutnet, level="class", remove.loops=TRUE)
#'
#' # Generate a network plot at the treatment level that colours nodes by class
#' plot(goutnet, v.color="class", remove.loops=TRUE)
#'
#' # Plot network in which node size is proportional to number of participants
#' alognet <- mb.network(alog_pcfb)
#' plot(alognet, v.scale=2)
#'
#' }
#'
#' @export
plot.mb.network <- function(x, edge.scale=1, label.distance=0,
                            level="treatment", remove.loops=FALSE, v.color="connect",
                            v.scale=NULL, layout=igraph::in_circle(), legend=TRUE,
                            legend.x="bottomleft", legend.y=NULL,
                            ...)
{
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mb.network", add=argcheck)
  checkmate::assertClass(layout, "igraph_layout_spec", add=argcheck)
  checkmate::assertNumeric(edge.scale, finite=TRUE, len=1, add=argcheck)
  checkmate::assertNumeric(label.distance, finite=TRUE, len=1, add=argcheck)
  checkmate::assertNumeric(v.scale, lower = 0, finite=TRUE, null.ok=TRUE, len=1, add=argcheck)
  checkmate::assertChoice(level, choices = c("treatment", "class"), add=argcheck)
  checkmate::assertChoice(v.color, choices = c("connect", "class"), add=argcheck)
  checkmate::assertLogical(remove.loops, len=1, add=argcheck)
  checkmate::assertLogical(legend, len=1, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Generate comparisons (using get.latest.time and mb.contrast?
  df <- get.latest.time(x)$data.ab


  # Check if level="class" that classes are present in dataset
  if (level=="class") {
    if (!("classes" %in% names(x))) {
      stop("`level` has been set to `class` but there are no class codes given in the dataset")
    }
    nodes <- x[["classes"]]

    # Replace treatment column with class codes
    df$treatment <- df$class

  } else if (level=="treatment") {
    nodes <- x[["treatments"]]
  }

  # Calculate participant numbers (if v.scale not NULL)
  if (!is.null(v.scale)) {
    if (!("n" %in% names(df))) {
      stop("`n` not included as a column in dataset. Vertices/nodes will all be scaled to be the same size.")
    }

    data.early <- x$data.ab[x$data.ab$fupcount==1,]
    size.vec <- vector()
    for (i in seq_along(nodes)) {
      size.vec <- append(size.vec, sum(data.early$n[data.early[[level]]==i]))
    }
    # Scale size.vec by the max node.size
    size.vec <- size.vec/ (max(size.vec)/20)

    node.size <-
      stats::setNames(size.vec, nodes)
    node.size <- as.matrix(node.size*v.scale)
  } else {
    node.size <- NULL
  }


  comparisons <- mb.comparisons(df)
  #treatments <- x[["treatments"]]

  # Code to make graph.create as an MBNMA command if needed
  g <- igraph::graph.empty()
  g <- g + igraph::vertex(nodes)
  ed <- t(matrix(c(comparisons[["t1"]], comparisons[["t2"]]), ncol = 2))
  ed <- factor(as.vector(ed), labels=nodes)
  edges <- igraph::edges(ed, weight = comparisons[["nr"]], arrow.mode=0)
  #edges <- igraph::edges(as.vector(ed), weight = comparisons[["nr"]], arrow.mode=0)
  g <- g + edges
  #g <- igraph::add.edges(g, edges[[1]], weight = (comparisons[["nr"]]*10), arrow.mode=0)

  # Check network is connected and produce warning message if not
  disconnects <- check.network(g)
  if (v.color=="connect") {
    igraph::V(g)$color <- "SkyBlue2"
    igraph::V(g)$color[which(names(igraph::V(g)) %in% disconnects)] <- "white"
  } else if (v.color=="class") {
    if (!("classes" %in% names(x))) {
      stop("`level` has been set to `class` but there are no class codes given in the dataset")
    }

    # Get large vector of distinct colours using Rcolorbrewer
    cols <- genmaxcols()

    if (level=="treatment") {
      igraph::V(g)$color <- cols[x$classkey$class]
    } else if (level=="class") {
      igraph::V(g)$color <- cols[unique(x$classkey$class)]
    }
  }

  # Add attributes
  igraph::V(g)$label.dist <- label.distance
  if (!is.null(node.size)) {igraph::V(g)$size <- node.size}
  igraph::E(g)$width <- edge.scale * comparisons[["nr"]]
  igraph::E(g)$color <- "grey20"

  if (remove.loops==TRUE) {
    g <- igraph::simplify(g, remove.multiple = FALSE, remove.loops = TRUE)
  }

  # Change label locations if layout_in_circle
  laycheck <- as.character(layout)[2]
  if (any(
    grepl("layout_in_circle", laycheck) |
    grepl("layout_as_star", laycheck))) {
    lab.locs <- radian.rescale(x=seq(1:length(nodes)), direction=-1, start=0)
    igraph::V(g)$label.degree <- lab.locs
  }

  # Plot netgraph
  g$layout <- igraph::layout_(g, layout)
  igraph::plot.igraph(g, ...)

  # Add legend
  if (v.color=="class" & legend==TRUE) {
    graphics::legend(x=legend.x, y=legend.y, legend=x$classes, pt.bg=unique(igraph::V(g)$color), pch=21, pt.cex=1.5, cex=0.8)
  }

  return(invisible(g))
}




#' Print mb.network information to the console
#'
#' @param x An object of class `mb.network`.
#' @param ... further arguments passed to or from other methods
#'
#' @return Prints the contents of `x` to the console.
#'
#' @export
print.mb.network <- function(x,...) {
  nn <- names(x)
  ll <- length(x)
  if (length(nn) != ll)
    nn <- paste("Component", seq.int(ll))
  for (i in seq_len(ll)) {
    cat(nn[i], ":\n")
    if (is.data.frame((x[[i]]))) {
      print(x[[i]], max=ncol(x[[i]])*6, ...)
    } else {
      print(x[[i]], ...)
    }
    cat("\n")
  }
  invisible(x)
}




#' Print summary mb.network information to the console
#'
#' @param object An object of class `mb.network`.
#' @param ... further arguments passed to or from other methods
#'
#' @return Prints summary details of `x` to the console.
#'
#' @export
summary.mb.network <- function(object,...) {

  # Print basic network statistics to the console
  cat(crayon::underline(crayon::bold("Description:", object$description, "\n")))
  cat("Number of studies:", length(unique(object$data.ab$studyID)), "\n")
  cat("Number of treatments:", length(object$treatments), "\n")

  if ("classes" %in% names(object)) {
    cat("Number of classes:", length(object$classes), "\n")
  }

  # Single line per study
  study.df <- object$data.ab %>%
    dplyr::group_by(studyID) %>%
    dplyr::slice_tail()

  # FUPs per study
  cat("Median (min, max) follow-up measurements per study: ", stats::median(study.df$fups),
      " (", min(study.df$fups), ", ", max(study.df$fups), ")\n", sep="")

  # Study duration
  cat("Median (min, max) study duration: ", stats::median(study.df$time),
      " (", min(study.df$time), ", ", max(study.df$time), ")\n", sep="")

  # Single line per study & treatment
  treat.df <- object$data.ab %>%
    dplyr::group_by(studyID, treatment) %>%
    dplyr::slice_tail() %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarise(med=stats::median(fups), min=min(fups), max=max(fups),
                     medt=stats::median(time), mint=min(time), maxt=max(time)) %>%
    dplyr::arrange(treatment)

  treat.df <- data.frame("Treatment"=as.character(factor(treat.df$treatment, labels=object$treatments)),
                         "FUPstudy"=paste0(treat.df$med, " (", treat.df$min, ", ", treat.df$max, ")"),
                         "TIMEstudy"=paste0(treat.df$medt, " (", treat.df$mint, ", ", treat.df$maxt, ")")
  )

  # Print data frame
  print(knitr::kable(treat.df, col.names = c("Treatment", "Follow-up measurements per study", "Study duration")))
  cat("\n")

  # Check network is connected at treatment-level
  g <- suppressWarnings(plot.invisible(object, level="treatment"))
  connects <- is.finite(igraph::shortest.paths(igraph::as.undirected(g),
                                               to=1))
  if (any(connects==FALSE)) {
    cat("Treatment-level network is", crayon::bold(crayon::red("DISCONNECTED"), "\n"))
  } else {
    cat("Treatment-level network is", crayon::bold(crayon::green("CONNECTED"), "\n"))
  }

  invisible(object)
}


