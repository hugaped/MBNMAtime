# Graphical functions
# Author: Hugo Pedder
# Date created: 2018-09-10



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
#' @param ... Options for plotting in `igraph`.
#'
#' @details The S3 method `plot()` on an `mb.network` object generates a
#'   network plot that shows how different treatments are connected within the
#'   network via study comparisons. This can be used to identify how direct and
#'   indirect evidence are informing different treatment comparisons. Depends on
#'   \code{\link[igraph]{igraph}}.
#'
#' @examples
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
#' @export
plot.mb.network <- function(x, edge.scale=1, label.distance=0,
                           level="treatment", remove.loops=FALSE, v.color="connect",
                           v.scale=NULL, layout=igraph::in_circle(),
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
  checkmate::reportAssertions(argcheck)

  # Generate comparisons (using get.latest.time and mb.contrast?
  data <- get.latest.time(x)


  # Check if level="class" that classes are present in dataset
  if (level=="class") {
    if (!("classes" %in% names(x))) {
      stop("`level` has been set to `class` but there are no class codes given in the dataset")
    }
    nodes <- x[["classes"]]

    # Replace treatment column with class codes
    data$treatment <- data$class

  } else if (level=="treatment") {
    nodes <- x[["treatments"]]
  }

  # Calculate participant numbers (if v.scale not NULL)
  if (!is.null(v.scale)) {
    if (!("N" %in% names(data))) {
      stop("`N` not included as a column in dataset. Vertices/nodes will all be scaled to be the same size.")
    }

    data.early <- x$data.ab[x$data.ab$fupcount==1,]
    size.vec <- vector()
    for (i in seq_along(nodes)) {
      size.vec <- append(size.vec, sum(data.early$N[data.early[[level]]==i]))
    }
    # Scale size.vec by the max node.size
    size.vec <- size.vec/ (max(size.vec)/20)

    node.size <-
      stats::setNames(size.vec, nodes)
    node.size <- as.matrix(node.size*v.scale)
  } else {
    node.size <- NULL
    }


  comparisons <- mb.comparisons(data)
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
  layout <- igraph::layout_(g, layout)
  igraph::plot.igraph(g,
                      layout = layout,
                      ...
  )

  # # Plot netgraph
  # if (layout_in_circle==TRUE) {
  #   lab.locs <- radian.rescale(x=seq(1:length(nodes)), direction=-1, start=0)
  #   igraph::V(g)$label.degree <- lab.locs
  #   igraph::plot.igraph(g,
  #        layout = igraph::layout_in_circle(g),
  #        ...
  #   )
  # } else {
  #   igraph::plot.igraph(g, ...)
  # }

  return(invisible(g))
}





#' Get large vector of distinct colours using Rcolorbrewer
genmaxcols <- function() {

  cols1 <- RColorBrewer::brewer.pal(9, "Set1")
  cols2 <- RColorBrewer::brewer.pal(9, "Pastel1")
  cols <- c(rbind(cols1, cols2))

  cols1 <- RColorBrewer::brewer.pal(8, "Set2")
  cols2 <- RColorBrewer::brewer.pal(8, "Pastel2")
  cols <- c(cols, c(rbind(cols1, cols2)))

  cols <- c(cols, RColorBrewer::brewer.pal(12, "Set3"))

  return(cols)
}









#' Calculate position of label with respect to vertex location within a circle
#'
#' Useful for graphs drawn using `igraph` to reposition labels relative to vertices when vertices
#' are laid out in a circle (as is common in network plots). `igraph` interprets position within
#' `vertex.label.degree` as radians, so it is necessary to convert locations into radian values. This
#' is the main role of this function.
#'
#' @param x A numeric vector of positions around a circle, typically sequentially numbered.
#' @param start A number giving the offset from 12 o'clock in radians for the label locations.
#' @param direction Either `1` for clockwise numbering (based on the order of `x`) or `-1` for
#' anti-clockwise.
#'
#' @examples
#' MBNMAtime:::radian.rescale(c(1:10), start=0, direction=1)
#'
#' @references
#' https://gist.github.com/kjhealy/834774/a4e677401fd6e4c319135dabeaf9894393f9392c
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}




#' Plots predicted responses from a time-course MBNMA model
#'
#' @param x An object of class `"mb.predict"` generated by
#'   `predict("mbnma")`
#' @param disp.obs A boolean object to indicate whether to show shaded sections
#'   of the plot for where there is observed data (`TRUE`) or not (`FALSE`)
#' @param overlay.ref A boolean object indicating whether to overlay a line
#'   showing the median network reference treatment response over time on the
#'   plot (`TRUE`) or not (`FALSE`). The network reference treatment (treatment
#'   1) must be included in `predict` for this to display the network reference
#'   treatment properly.
#' @param col A character indicating the colour to use for shading if `disp.obs`
#'   is set to `TRUE`. Can be either `"blue"`, `"green"`, or `"red"`
#' @param max.col.scale Rarely requires adjustment. The maximum count of
#'   observations (therefore the darkest shaded color) only used if `disp.obs` is
#'   used. This allows consistency of shading between multiple plotted graphs.
#'   It should always be at least as high as the maximum count of observations
#'   plotted
#' @param ... Arguments for `ggplot`
#' @inheritParams plot.mb.rank
#'
#' @details For the S3 method `plot()`, if `disp.obs` is set to `TRUE` it is
#'   advisable to ensure predictions in `predict` are estimated using an even
#'   sequence of time points to avoid misrepresentation of shaded densities.
#'   Shaded counts of observations will be relative to the treatment plotted in
#'   each panel rather than to the network reference treatment if `disp.obs` is
#'   set to `TRUE`.
#'
#' @examples
#' \donttest{
#' # Create an mb.network object from a dataset
#' network <- mb.network(alog_pcfb)
#'
#' # Run an MBNMA model with an Emax time-course
#' emax <- mb.emax(network,
#'   emax=list(pool="rel", method="common"),
#'   et50=list(pool="rel", method="common"))
#'
#' # Predict responses using the original dataset to estimate the network reference
#' #treatment response
#' df.ref <- alog_pcfb[alog_pcfb$treatment=="placebo",]
#' predict <- predict(emax, times=c(0:15), baseline=10, ref.estimate=df.ref)
#'
#' # Plot the predicted responses with observations displayed on plot as green shading
#' plot(predict, disp.obs=TRUE, overlay.ref=FALSE, col="green")
#'
#' # Plot the predicted responses with the median network reference treatment response overlayed
#' #on the plot
#' plot(predict, disp.obs=FALSE, overlay.ref=TRUE)
#' }
#'
#' @export
plot.mb.predict <- function(x, disp.obs=FALSE, overlay.ref=TRUE,
                           col="blue", max.col.scale=NULL, treat.labs=NULL, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mb.predict", add=argcheck)
  checkmate::assertLogical(disp.obs, len=1, add=argcheck)
  checkmate::assertLogical(overlay.ref, len=1, add=argcheck)
  checkmate::reportAssertions(argcheck)

  pred <- x[["summary"]]

  data <- pred[[1]]
  data[["treat"]] <- rep(0, nrow(data))
  data <- data[0,]
  for (i in seq_along(pred)) {
    cut <- pred[[i]]
    #cut[["treat"]] <- rep(as.numeric(names(pred)[i]), nrow(cut))
    cut[["treat"]] <- rep(names(pred)[i], nrow(cut))
    data <- rbind(data, cut)
  }

  # Keep only relevant columns
  data <- data[, which(names(data) %in%
                         c("time", "2.5%", "50%", "97.5%", "treat"))]


  # Add treatment labels
  if (!is.null(treat.labs)) {
    data$treat <- factor(data$treat, labels=treat.labs)
  } else if (is.null(treat.labs)) {
    treat.labs <- names(pred$summary)
  }

  # Required for overlaying ref treatment effect
  if (overlay.ref==TRUE) {
    ref.treat <- x$network$treatments[1]

    if (!(ref.treat %in% names(pred))) {
      stop(paste0("Reference treatment (", ref.treat, ") must be included in `x` in order for it to be plotted"))
    }

    #data[["ref.median"]] <- rep(pred[["1"]][[/"50%"]], length(pred))
    data[["ref.median"]] <- rep(pred[[ref.treat]][["50%"]], length(pred))
    #data <- data[data$treat!=1,]
    data <- data[data$treat!=ref.treat,]
    treat.labs <- treat.labs[treat.labs!=ref.treat]
    x[["summary"]][[ref.treat]] <- NULL
  }


  # Plot ggplot axes
  g <- ggplot2::ggplot(data, ggplot2::aes(x=time, y=`50%`, ymin=`2.5%`, ymax=`97.5%`), ...)

  # Add shaded regions for observations in original dataset
  if (disp.obs==TRUE) {
    #g <- disp.obs(g, network, x, col=col, max.col.scale=max.col.scale)
    g <- disp.obs(g, predict=x, col=col, max.col.scale=max.col.scale)
  }

  # Overlay reference treatment effect
  if (overlay.ref==TRUE) {
    g <- g + ggplot2::geom_line(ggplot2::aes(y=ref.median, colour="Reference Mean"), size=1)
    message(paste0("Reference treatment in plots is ", ref.treat))
  }

  # Add overlayed lines and legends
  g <- g + ggplot2::geom_line(ggplot2::aes(linetype="Predicted Mean")) +
    ggplot2::geom_line(ggplot2::aes(y=`2.5%`, linetype="95% CrI")) +
    ggplot2::geom_line(ggplot2::aes(y=`97.5%`, linetype="95% CrI"))

  g <- g + ggplot2::facet_wrap(~factor(treat)) +
    ggplot2::labs(y="Predicted response", x="Time")

  g <- g + ggplot2::scale_linetype_manual(name="",
                                 values=c("Predicted Mean"="solid",
                                          "95% CrI"="dashed"))

  g <- g + ggplot2::scale_color_manual(name="",
                              values=c("Reference Mean"="red"))

  g <- g + theme_mbnma()

  return(g)
}




#' Overlays observations as shaded regions on a time-course plot
#'
#' @param g An object of `c("gg", "ggplot")` that is a plot of predicted responses from an MBNMA model.
#' @inheritParams plot.mb.predict
#'
#' @noRd
disp.obs <- function(g, predict, col="blue", max.col.scale=NULL) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(g, c("gg", "ggplot"), add=argcheck)
  #checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertClass(predict, "mb.predict", add=argcheck)
  checkmate::assertInt(max.col.scale, lower=1, null.ok = TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  #treats <- as.numeric(names(predict[["summary"]]))
  #treats <- which(predict$mbnma$network$treatments %in% names(predict$summary))
  treats <- names(predict[["summary"]])
  #data.ab <- jagstonetwork(predict$mbnma)
  data.ab <- predict$network$data.ab

  #raw.data <- network[["data.ab"]]
  raw.data <- dplyr::arrange(data.ab, studyID, fupcount, arm)
  raw.data$studyID <- as.character(raw.data$studyID)

  predict.data <- predict[["summary"]][[1]]
  predict.data[["treat"]] <- NA
  predict.data[["count"]] <- NA
  predict.data <- predict.data[0,]
  for (i in seq_along(predict[["summary"]])) {
    temp <- predict[["summary"]][[i]]
    # temp[["treat"]] <- rep(as.numeric(names(predict[["summary"]][i])),
    #                        nrow(temp))
    # temp[["treat"]] <- rep(which(predict$mbnma$network$treatments %in% names(predict$summary)[i]),
    #                        nrow(temp))
    temp[["treat"]] <- rep(names(predict$summary)[i],
                           nrow(temp))
    temp[["count"]] <- rep(NA, nrow(temp))

    dplyr::arrange(temp, time)
    if (temp[["time"]][1]>0) {
      temp[["count"]][1] <- nrow(raw.data[raw.data$treatment==which(predict$network$treatments %in% temp[["treat"]][1]) &
                                            raw.data[["time"]]<=temp[["time"]][1] &
                                            raw.data[["time"]]>0
                                          ,])
    } else {temp[["count"]][1] <- 0  }

    for (k in 2:nrow(temp)) {
      temp[["count"]][k] <- nrow(raw.data[raw.data$treatment==which(predict$network$treatments %in% temp[["treat"]][k]) &
                                            raw.data[["time"]]<=temp[["time"]][k] &
                                            raw.data[["time"]]>temp[["time"]][k-1]
                                          ,])
    }

    predict.data <- rbind(predict.data, temp)
  }

  # Check max.col.scale
  n.cut <- max(predict.data$count)
  if (!is.null(max.col.scale)) {
    if (!is.numeric(max.col.scale)) {
      stop("max.col.scale must be a number greater than the maximum number of observed counts in the plotted treatments")
    }
    if (n.cut > max.col.scale) {
      stop("max.col.scale must be a number greater than the maximum number of observed counts in the plotted treatments")
    }
    n.cut <- max.col.scale
  }

  # Generate colours
  #cols <- col.scale(n.cut=max(predict.data$count), col=col)
  cols <- alpha.scale(n.cut=n.cut, col=col)

  for (treat in seq_along(treats)) {
    subset <- predict.data[predict.data$treat==treats[treat],]

    # Start assuming lowest time = 0
    for (m in 2:nrow(subset)) {
      g <- g + ggplot2::geom_ribbon(data=subset[subset$time<=subset$time[m] &
                                         subset$time>=subset$time[m-1]
                                       ,],
                           ggplot2::aes(x=time, ymin=`2.5%`, ymax=`97.5%`),
                           fill=cols[subset$count[m]+1])
    }

  }

  return(g)
}





#' Generates a scale of colours
#'
#' FUNCTION IS DEPRACATED
#' @noRd
col.scale <- function(n.cut, col="blue",
                      rgb.start=NULL, rgb.end=NULL) {
  # Set colour intensities
  if (col=="blue") {
    max.rgb <- c(0,0,200)
    min.rgb <- c(200,200,250)
  } else if (col=="red") {
    max.rgb <- c(120,0,0)
    min.rgb <- c(250,200,200)
  } else if (col=="green") {
    max.rgb <- c(0,100,0)
    min.rgb <- c(200,250,200)
  }

  if (!is.null(rgb.start) & !is.null(rgb.end)) {
    max.rgb <- rgb.start
    min.rgb <- rgb.end
  }

  #min.rgb <- rgb.end
  cut <- (min.rgb-max.rgb)/(n.cut-1)

  #rgb <- c(255,255,255)
  col.mat <- c(255,255,255)
  rgb <- min.rgb
  col.mat <- cbind(col.mat, rgb)
  for (i in 1:n.cut) {
    rgb <- rgb - cut
    col.mat <- cbind(col.mat, abs(rgb))
  }
  #col.mat <- rbind(col.mat, c(0, rep(255, (n.cut+1))))
  col.mat <- rbind(col.mat, c(0, rep(255, (ncol(col.mat)-1))))

  hexcol <- vector()
  for (i in 1:(n.cut)) {
    hexcol <- append(hexcol, rgb(col.mat[1,i], col.mat[2,i],
                                 col.mat[3,i], col.mat[4,i], maxColorValue=255
                                 ))
  }
  return(hexcol)
}





#' Generates colours with varying degrees of transparency
#'
#' Identical to MBNMAdose
#'
#' @param ncut A number indicating the number of different counts in the dataset
#' @param col Colour to use for shading
#'
#' @noRd
alpha.scale <- function(n.cut, col="blue") {
  # Run checks
  checkmate::assertIntegerish(n.cut, lower=1, len=1)

  # Set colour intensities
  if (is.character(col)) {
    if (col=="blue") {
      hue <- c(0,0,200)
    } else if (col=="red") {
      hue <- c(120,0,0)
    } else if (col=="green") {
      hue <- c(0,100,0)
    } else {
      stop("Permitted colours are `blue`, `red`, `green`, or an RGB code")
    }
  } else if (is.numeric(col)) {
    if (length(col)!=3) {
      stop("Specified RGB code must have length 3")
    }
    if (any(col>255) | any(col<0)) {
      stop("RGB values must be between 0 and 255")
    }
    hue <- col
  }

  #min.rgb <- rgb.end
  cut <- (255-0)/(n.cut)
  alpha <- 0
  alpha.vec <- alpha

  for (i in 1:n.cut) {
    alpha <- alpha+cut
    alpha.vec <- append(alpha.vec, alpha)
  }

  hexcol <- vector()
  for (i in 1:(n.cut+1)) {
    hexcol <- append(hexcol, grDevices::rgb(hue[1], hue[2],
                                 hue[3], alpha.vec[i], maxColorValue=255
    ))
  }

  return(hexcol)
}








#' Plot histograms of rankings from MBNMA models
#' @param x An object of class `"mb.rank"` generated by `rank.mbnma()`
#' @param treat.labs A vector of treatment labels in the same order as treatment codes.
#' Easiest to use treatment labels stored by `mb.network()`
#' @param ... Arguments to be sent to `ggplot::geom_bar()`
#' @inheritParams rank.mbnma
#'
#' @return A series of histograms that show rankings for each treatment/agent/prediction, with a
#' separate panel for each parameter
#' The object returned is a list containing a separate element for each parameter in `params`
#' which is an object of class `c("gg", "ggplot")`.
#'
#' @examples
#' \donttest{
#' # Create an mb.network object from a dataset
#' network <- mb.network(osteopain)
#'
#' # Run an MBNMA model with an Emax time-course
#' emax <- mb.emax(network,
#'   emax=list(pool="rel", method="common"),
#'   et50=list(pool="arm", method="common"),
#'   positive.scale=TRUE)
#'
#' # Calculate treatment rankings
#' ranks <- rank(emax,
#'   param=c("auc", "d.emax", "beta.et50"),
#'   int.range=c(0,15),
#'   treats=c(1:10), n.iter=500)
#'
#' # Plot histograms for ranking by AUC
#' plot(ranks, param="auc")
#'
#' # Plot histograms for ranking by d.emax
#' plot(ranks, param="d.emax")
#' }
#'
#' @export
plot.mb.rank <- function(x, params=NULL, treat.labs=NULL, ...) {
  # ... are commands to be sent to geom_histogram

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mb.rank", add=argcheck)
  checkmate::assertCharacter(params, null.ok=TRUE, add=argcheck)
  checkmate::assertCharacter(treat.labs, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  output <- list()

  if (is.null(params)) {
    params <- names(x)
  }

  if (is.null(treat.labs)) {
    treat.labs <- as.character(colnames(x[[1]]$rank.matrix))
  } else if (!is.null(treat.labs)) {
    if (length(treat.labs)!=ncol(x[[1]]$rank.matrix)) {
      stop("`treat.labs` must be the same length as the number of treatments that have been ranked in `x`")
    }
  }

  for (param in seq_along(params)) {

    rank.mat <- x[[params[param]]]$rank.matrix
    #treats <- colnames(rank.mat)
    treats <- c(1:ncol(rank.mat))

    ranks.param <- vector()
    treat <- vector()
    for (i in seq_along(treats)) {
      treat <- append(treat, rep(treats[i], nrow(rank.mat)))
      ranks.param <- append(ranks.param, rank.mat[,i])
    }
    data <- data.frame("ranks"=ranks.param, "treat"=treat)

    # if (!is.null(treat.labs)) {
    #   data$treat <- factor(data$treat, labels=treat.labs)
    # } else {
    #   data$treat <- factor(as.numeric(as.character(data$treat)))
    # }
    data$treat <- factor(data$treat, labels=treat.labs)

    g <- ggplot2::ggplot(data, ggplot2::aes(x=data$ranks)) +
      ggplot2::geom_bar(...) +
      ggplot2::xlab("Rank (1 = best)") +
      ggplot2::ylab("MCMC iterations") +
      ggplot2::facet_wrap(~treat) +
      ggplot2::ggtitle(params[param]) +
      theme_mbnma()

    graphics::plot(g)

    output[[params[param]]] <- g
  }

  return(invisible(output))
}









#' Forest plot for results from time-course MBNMA models
#'
#' Generates a forest plot for time-course parameters of interest from results from time-course MBNMA models.
#'
#' @inheritParams devplot
#' @param x An S3 object of class `"mbnma"` generated by running
#' a time-course MBNMA model
#' @param params A character vector of time-course parameters to plot.
#' Parameters must be given the same name as monitored nodes in `mbnma` and must vary by treatment or class. Can be set to
#' `NULL` to include all available time-course parameters estimated by `mbnma`.
#' @param treat.labs A character vector of treatment labels. If left as `NULL` (the default) then
#' labels will be used as defined in the data.
#' @param class.labs A character vector of class labels if `mbnma` was modelled using class effects
#' If left as `NULL` (the default) then labels will be used as defined in the data.
#' @param ... Arguments to be sent to `ggplot`
#'
#' @return A forest plot of class `c("gg", "ggplot")` that has separate panels for different time-course parameters
#'
#' @examples
#'\donttest{
#' # Create an mb.network object from a dataset
#' network <- mb.network(alog_pcfb)
#'
#' # Run an MBNMA model with an Emax time-course
#' emax <- mb.emax(network,
#'   emax=list(pool="rel", method="common"),
#'   et50=list(pool="rel", method="common"))
#'
#' # Generate forest plot
#' plot(emax)
#'
#' # Plot results for only one time-course parameter
#' plot(emax, params="d.emax")
#' }
#' @export
plot.mbnma <- function(x, params=NULL, treat.labs=NULL, class.labs=NULL, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mbnma", add=argcheck)
  checkmate::assertCharacter(params, null.ok=TRUE, add=argcheck)
  checkmate::assertCharacter(treat.labs, null.ok=TRUE, add=argcheck)
  checkmate::assertCharacter(class.labs, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check that specified params are monitored in model
  if (!all(params %in% x[["parameters.to.save"]])) {
    stop(paste0("Variable 'params': Must contain elements of set {", paste(x[["parameters.to.save"]], collapse = ", "), "}"))
  }

  # Check that specified params vary over treatment
  for (i in seq_along(params)) {
    if (length(x[["BUGSoutput"]][["median"]][[params[i]]])<=1) {
      stop(paste0(params[i], " does not vary over treatment/class and cannot be plotted"))
    }
  }


  # Add all available params if is.null(params)
  if (is.null(params)) {
    params <- vector()

    # Add d
    params <- append(params,
                     x[["parameters.to.save"]][grep("^d\\.", x[["parameters.to.save"]])]
    )

    # Add D
    params <- append(params,
                     x[["parameters.to.save"]][grep("^D\\.", x[["parameters.to.save"]])]
    )

    # Add BETA
    params <- append(params,
                     x[["parameters.to.save"]][grep("^BETA\\.", x[["parameters.to.save"]])]
    )

    # Add beta
    for (i in seq_along(x[["BUGSoutput"]][["root.short"]])) {
      if (grepl("^beta\\.", x[["BUGSoutput"]][["root.short"]][i]) &
          length(x[["BUGSoutput"]][["long.short"]][[i]])>1) {
        params <- append(params, x[["BUGSoutput"]][["root.short"]][i])
      }
    }
    if (length(params)==0) {
      stop("No time-course consistency parameters can be identified from the model")
    }
  }


  # Compile parameter data into one data frame
  mb.sum <- as.data.frame(x[["BUGSoutput"]][["summary"]])
  plotdata <- mb.sum[0,]
  for (i in seq_along(params)) {
    paramdata <- mb.sum[grepl(paste0("^", params[i]),rownames(mb.sum)),]
    paramdata[["timeparam"]] <- rep(params[i], nrow(paramdata))
    plotdata <- rbind(plotdata, paramdata)
  }
  plotdata[["param"]] <- as.numeric(gsub("(.+\\[)([0-9]+)(\\])", "\\2", rownames(plotdata)))

  # Change param labels for agents
  treatdat <- plotdata[grepl("^d\\.", rownames(plotdata)) | grepl("^beta\\.", rownames(plotdata)),]
  if (!is.null(treat.labs)) {
    treatcodes <- as.numeric(gsub("(^.+\\[)([0-9]+)(\\])", "\\2", rownames(treatdat)))
    if (length(treat.labs)!=max(treatcodes)) {
      stop("`treat.labs` length does not equal number of treatments that have been modelled for this time-course parameter")
    } else {
      t.labs <- treat.labs[sort(unique(treatcodes))]
    }
  } else if ("treatments" %in% names(x)) {
    t.labs <- x$network[["treatments"]]
  } else {
    t.labs <- sort(unique(treatdat$param))
  }

  # Change param labels for classes
  classdat <- plotdata[grepl("^D\\.", rownames(plotdata)) | grepl("^BETA\\.", rownames(plotdata)),]
  c.labs <- vector()
  if (nrow(classdat)!=0) {
    if (!is.null(class.labs)) {
      classcodes <- as.numeric(gsub("(^.+\\[)([0-9]+)(\\])", "\\2", rownames(classdat)))
      c.labs <- class.labs[classcodes]
    } else if ("classes" %in% names(x)) {
      c.labs <- x[["classes"]][x[["classes"]]!="Placebo"]
    } else {
      c.labs <- sort(unique(classdat$param))
    }
  }

  # Increase param number for classes
  ntreat <- ifelse(nrow(treatdat)>0, max(treatdat$param), 0)
  plotdata$param[grepl("^D\\.", rownames(plotdata)) | grepl("^BETA\\.", rownames(plotdata))] <-
    plotdata$param[grepl("^D\\.", rownames(plotdata)) | grepl("^BETA\\.", rownames(plotdata))] + ntreat

  # Attach labels
  if (nrow(treatdat)>0) {
    all.labs <- c(t.labs, c.labs)
  } else {all.labs <- c.labs}
  plotdata$param <- factor(plotdata$param, labels=all.labs)

  if (any(is.na(levels(plotdata$param)))) {
    stop("`treat.labs` or `class.labs` have not been specified correctly")
  }

  g <- ggplot2::ggplot(plotdata, ggplot2::aes(y=`50%`, x=factor(plotdata$param))) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=`2.5%`, ymax=`97.5%`)) +
    ggplot2::coord_flip()

  g <- g + ggplot2::facet_wrap(~timeparam, scales="free")

  # Axis labels
  g <- g + ggplot2::xlab("Treatment / Class") +
    ggplot2::ylab("Effect size") +
    theme_mbnma()

  graphics::plot(g, ...)
  return(invisible(g))
}





#' Check if all nodes in the network are connected (identical to MBNMAdose)
#' @noRd
check.network <- function(g, reference=1) {
  connects <- is.finite(igraph::shortest.paths(igraph::as.undirected(g),
                                               to=reference))
  treats <- rownames(connects)[connects==FALSE]

  if (length(treats>0)) {
    warning(paste0("The following treatments/agents are not connected to the network reference:\n",
                   paste(treats, collapse = "\n")))
  }
  return(treats)
}






#' Get data for nodesplit plots from mb.nodesplit
#' @noRd
nodesplit.plotdata <- function(nodesplit, type) {
  type <- paste0(type, ".plot")
  plotdata <- data.frame()
  time.param <- vector()
  comparison <- vector()
  for (param in seq_along(nodesplit[[1]])) {
    for (split in seq_along(nodesplit)) {
      temp <- nodesplit[[split]][[param]][[type]]$data
      temp$time.param <- rep(names(nodesplit[[1]])[param], nrow(temp))
      temp$comparison <- rep(names(nodesplit)[split], nrow(temp))
      plotdata <- rbind(plotdata, temp)
    }
  }
  return(plotdata)
}





#' Plot raw responses over time by treatment or class
#'
#' @param network An object of class `mb.network`.
#' @param plotby A character object that can take either `"arm"` to indicate that raw responses
#' should be plotted separately for each study arm, or `"rel"` to indicate that the relative
#' effects within each study should be plotted. In this way the time-course of both the absolute
#' effects and the relative effects can be examined.
#' @param ... Arguments to be sent to `ggplot`
#' @inheritParams plot.mb.network
#'
#' @return The function returns an object of `class(c("gg", "ggplot")`. Characteristics
#' of the object can therefore be amended as with other plots generated by `ggplot()`.
#' A message will indicate if data are assumed to be change from baseline (i.e. if there
#' are no responses in the data at time=0). In this case responses will be set to 0 at
#' time=0.
#'
#' @details Plots can be faceted by either treatment (`level="treatment"`) or class
#' (`level="class"`) to investigate similarity of treatment responses within classes/agents.
#' Points represent observed responses and lines connect between observations within the
#' same study and arm.
#'
#' @examples
#' # Make network
#' network <- mb.network(goutSUA_CFB)
#'
#' # Use timeplot to plot responses grouped by treatment
#' timeplot(network)
#'
#' # Use timeplot ot plot resposes grouped by class
#' timeplot(network, level="class")
#'
#' @export
timeplot <- function(network, level="treatment", plotby="arm", ...) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mb.network", add=argcheck)
  checkmate::assertChoice(level, choices = c("treatment", "class"), add=argcheck)
  checkmate::assertChoice(plotby, choices = c("arm", "rel"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (level=="class" & !("classes" %in% names(network))) {
    stop("`level` cannot be set to class if there is no class variable included in the dataset/network")
  }

  plotdata <- network$data.ab

  if (plotby=="arm") {
    # Identify if data is CFB or not
    base.dat <- plotdata[plotdata$fupcount==1,]
    if (any(base.dat$time!=0)) {
      base.dat$y <- rep(0, nrow(base.dat))
      base.dat$se <- rep(0, nrow(base.dat))
      base.dat$time <- rep(0, nrow(base.dat))

      plotdata <- rbind(base.dat, plotdata)
      plotdata <- dplyr::arrange(plotdata, studyID, time, treatment)

      # Do not run this function with pylr loaded!!
      plotdata <- plotdata %>%
        dplyr::group_by(studyID, time) %>%
        dplyr::mutate(arm = sequence(dplyr::n()))

      message(cat("Absence of observations with time=0 in network - data assumed to be change from baseline:",
                  "plotting response=0 at time=0", sep="\n"))
    }

    g <- ggplot2::ggplot(plotdata,
                         ggplot2::aes(x=plotdata$time, y=plotdata$y,
                                      group=(paste(plotdata$studyID, plotdata$arm, sep="_")))) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size=1)

    if (level=="treatment") {
      g <- g + ggplot2::facet_wrap(~factor(.data$treatment, labels=network$treatments))
      #g <- g + ggplot2::facet_wrap(ggplot2::vars(treatment))
    } else if (level=="class") {
      g <- g + ggplot2::facet_wrap(~factor(.data$class, labels=network$classes))
    }

  } else if (plotby=="rel") {

    if (level=="treatment") {

      diffs <- plotdata %>%
        dplyr::mutate(Rx.Name = factor(network$treatments[.data$treatment], levels=network$treatments))

    } else if (level=="class") {

      diffs <- plotdata %>%
        dplyr::mutate(Rx.Name = factor(network$classes[.data$class], levels=network$classes))
    }

    diffs <- diffs %>%
      dplyr::inner_join(diffs, by=c("studyID", "time")) %>%
      dplyr::filter(.data$treatment.x < .data$treatment.y) %>%
      dplyr::mutate(pairDiff = .data$y.x - .data$y.y)

    diffs <- diffs %>%
      dplyr::bind_rows(diffs %>%
                         dplyr::group_by(.data$studyID, .data$arm.x, .data$arm.y) %>%
                         dplyr::slice(1) %>%
                         dplyr::mutate(time=0, pairDiff=0))

    g <- ggplot2::ggplot(data=diffs, ggplot2::aes(x=.data$time, y=.data$pairDiff, group=.data$studyID)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_grid(rows=ggplot2::vars(.data$Rx.Name.y), cols=ggplot2::vars(.data$Rx.Name.x))

  }


  g <- g + ggplot2::xlab("Time") + ggplot2::ylab("Response") +
    theme_mbnma()

  graphics::plot(g, ...)
  return(invisible(g))

}






#' Plot deviance contributions from an MBNMA model
#'
#' @param mbnma An S3 object of class `"mbnma"` generated by running
#' a time-course MBNMA model
#' @param dev.type Deviances to plot - can be either residual deviances (`"resdev"`, the
#' default) or deviances (`"dev"`)
#' @param plot.type Deviances can be plotted either as scatter points (`"scatter"` - using
#' `geom_point()`, the default) or as boxplots (`"box"`)
#' @param xaxis A character object that indicates whether deviance contributions should be plotted
#' by time (`"time"`) or by follow-up count (`"fup"`)
#' @param facet A boolean object that indicates whether or not to facet by treatment
#' @param n.iter The number of iterations to update the model whilst monitoring additional parameters (if necessary).
#' Must be a positive integer. Default is the value used in `mbnma`.
#' @param n.thin The thinning rate. Must be a positive integer. Default is the value used in `mbnma`.
#' @param ... Arguments to be sent to `ggplot2::geom_point()` or `ggplot2::geom_boxplot`
#' @inheritParams predict.mbnma
#'
#' @return Generates a plot of deviance contributions and returns a list containing the
#' plot (as an object of class `c("gg", "ggplot")`), and a data.frame of posterior mean
#' deviance/residual deviance contributions for each observation.
#'
#' @details
#' Deviances should only be plotted for models that have converged successfully. If deviance
#' contributions have not been monitored in `mbnma$parameters.to.save` then additional
#' iterations will have to be run to get results for these.
#'
#' Deviance contributions cannot be calculated for models with a multivariate likelihood (i.e.
#' those that account for correlation between observations) because the covariance matrix in these
#' models is treated as unknown (if `rho="estimate"`) and deviance contributions will be correlated.
#'
#' @examples
#' \donttest{
#' # Make network
#' network <- mb.network(alog_pcfb)
#'
#' # Run MBNMA
#' mbnma <- mb.quadratic(network)
#'
#' # Plot residual deviance contributions in a scatterplot
#' devplot(mbnma)
#'
#' # Plot deviance contributions in boxplots at each follow-up measurement
#' # Monitor for 500 additional iterations
#' devplot(mbnma, dev.type="dev", plot.type="box", xaxis="fup", n.iter=500)
#' }
#' @export
devplot <- function(mbnma, dev.type="resdev", plot.type="scatter",
                    xaxis="time", facet=TRUE,
                    n.iter=mbnma$BUGSoutput$n.iter, n.thin=mbnma$BUGSoutput$n.thin,
                    ...) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertChoice(dev.type, choices = c("dev", "resdev"), add=argcheck)
  checkmate::assertChoice(plot.type, choices = c("scatter", "box"), add=argcheck)
  checkmate::assertLogical(facet, add=argcheck)
  checkmate::assertChoice(xaxis, choices = c("time", "fup"), add=argcheck)
  checkmate::assertInt(n.iter, lower=1, null.ok = TRUE, add=argcheck)
  checkmate::assertInt(n.thin, lower=1, upper=n.iter, null.ok = TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check if !is.null(rho) (in which case cannot monitor deviances)
  # if (!is.null(mbnma$model.arg$rho)) {
  if (!is.null(mbnma$model.arg$rho)) {
    stop("Correlation between time points has been modelled using a multivariate likelihood.
         Deviances cannot be calculated for this model.")
  }

  if (!(dev.type %in% mbnma$parameters.to.save)) {
    msg <- paste0("`", dev.type, "` not monitored in mbnma$parameters.to.save.\n",
                  "additional iterations will be run in order to obtain results for `", dev.type, "`")
    message(msg)

    dev.df <- mb.update(mbnma, param=dev.type, n.iter=n.iter, n.thin=n.thin)

  } else {

    dev.df <- get.theta.dev(mbnma, param=dev.type)
  }

  # Changes fup to time if xaxis=time
  if (xaxis=="time") {
    dev.df$fup <- apply(dev.df, MARGIN=1, FUN=function(x) mbnma$model$data()$time[x[1], x[3]])
    xlab <- "Time"
  } else if (xaxis=="fup") {
    xlab <- "Follow-up count"
  }

  if (facet==TRUE) {
    dev.df$facet <- apply(dev.df, MARGIN=1, FUN=function(x) mbnma$model$data()$treat[x[1], x[2]])
  }

  # Plots the residual deviances over time grouped by study and arm
  if (plot.type=="scatter") {
    g <- ggplot2::ggplot(dev.df, ggplot2::aes(x=fup, y=mean), group=(paste(study, arm, sep="_"))) +
      ggplot2::geom_point(...)
  } else if (plot.type=="box") {
    g <- ggplot2::ggplot(dev.df, ggplot2::aes(x=factor(fup), y=mean)) +
      ggplot2::geom_boxplot(...)
  }

  # Add facets
  if (facet==TRUE) {
    g <- g + ggplot2::facet_wrap(~factor(facet, labels=mbnma$network$treatments), scales="free_x")
  }

  # Add axis labels
  g <- g + ggplot2::xlab(xlab) +
    ggplot2::ylab("Posterior mean") +
    theme_mbnma()

  graphics::plot(g)
  return(invisible(list("graph"=g, "dev.data"=dev.df)))
}



#' Plot fitted values from MBNMA model
#'
#' @param treat.labs A character vector of treatment labels with which to name graph panels.
#' Can use `mb.network()[["treatments"]]` with original dataset if in doubt.
#' @param disp.obs A boolean object to indicate whether raw data responses should be
#' plotted as points on the graph
#' @param ... Arguments to be sent to `ggplot2::geom_point()` or `ggplot2::geom_line()`
#' @inheritParams R2jags::jags
#' @inheritParams devplot
#'
#' @return Generates a plot of fitted values from the MBNMA model and returns a list containing
#' the plot (as an object of class `c("gg", "ggplot")`), and a data.frame of posterior mean
#' fitted values for each observation.
#'
#' @details
#' Fitted values should only be plotted for models that have converged successfully.
#' If fitted values (`theta`) have not been monitored in `mbnma$parameters.to.save`
#' then additional iterations will have to be run to get results for these.
#'
#' @examples
#' \donttest{
#' # Make network
#' network <- mb.network(osteopain)
#'
#' # Run MBNMA
#' mbnma <- mb.emax(network,
#'   emax=list(pool="rel", method="common"),
#'   et50=list(pool="const", method="common"))
#'
#' # Plot fitted values from the model with treatment labels
#' # Monitor fitted values for 500 additional iterations
#' fitplot(mbnma, treat.labs=network$treatments, n.iter=500)
#' }
#'
#' @export
fitplot <- function(mbnma, treat.labs=NULL, disp.obs=TRUE,
                    n.iter=mbnma$BUGSoutput$n.iter, n.thin=mbnma$BUGSoutput$n.thin, ...) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertLogical(disp.obs, add=argcheck)
  checkmate::assertInt(n.iter, lower=1, null.ok = TRUE, add=argcheck)
  checkmate::assertInt(n.thin, lower=1, upper=n.iter, null.ok = TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (!("theta" %in% mbnma$parameters.to.save)) {
    msg <- paste0("`theta` not monitored in mbnma$parameters.to.save.\n",
                  "additional iterations will be run in order to obtain results")
    message(msg)

    theta.df <- mb.update(mbnma, param="theta", n.iter=n.iter, n.thin=n.thin)
  } else {
    theta.df <- get.theta.dev(mbnma, param="theta")
  }

  # Get treatment codes
  theta.df <- get.treattimes(mbnma, theta.df, type="treat")

  # Get times
  theta.df <- get.treattimes(mbnma, theta.df, type="time")

  # Get y
  theta.df <- get.treattimes(mbnma, theta.df, type="y")

  # Check if intercept is modelled
  base.dat <- theta.df[theta.df$fup==1,]
  if (mbnma$model.arg$intercept==FALSE | all(base.dat$time!=0)) {
    if (any(base.dat$time!=0)) {
      message(cat("Absence of observations with time=0 in network - data assumed to be change from baseline:",
                  "plotting response=0 at time=0", sep="\n"))
    }

    base.dat$y <- rep(0, nrow(base.dat))
    base.dat$time <- rep(0, nrow(base.dat))
    base.dat$fup <- rep(0, nrow(base.dat))
    base.dat$mean <- rep(0, nrow(base.dat))

    theta.df <- rbind(theta.df, base.dat)
  }

  # Generate plot
  g <- ggplot2::ggplot(theta.df,
                       ggplot2::aes(x=time, y=mean, group=paste(study, arm, sep="_"))) +
    ggplot2::geom_line()

  # Overlay observed responses
  if (disp.obs==TRUE) {
    g <- g + ggplot2::geom_point(ggplot2::aes(y=y), size=1)
  }

  # Add facet labels
  if (!is.null(treat.labs)) {
    if (length(treat.labs)!=max(theta.df$treat)) {
      stop("treat.labs must be the same length as the number of treatment codes in the model")
    }
    g <- g + ggplot2::facet_wrap(~factor(treat, labels=treat.labs))
  } else {
    g <- g + ggplot2::facet_wrap(~factor(treat, labels=mbnma$network$treatments))
  }

  # Add axis labels
  g <- g + ggplot2::xlab("Time") +
    ggplot2::ylab("Response") +
    theme_mbnma()

  graphics::plot(g)

  #return(invisible(list("graph"=g, "theta.data"=theta.df)))
}



# Get matching treatments or times for a data.frame from an MBNMA model
get.treattimes <- function(mbnma, add.df, type="treat") {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertChoice(type, choices = c("time", "treat", "y"), add=argcheck)
  checkmate::assertDataFrame(add.df, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Get codes
  index.df <- reshape2::melt(mbnma$model$data()[[type]])
  index.df <- index.df[stats::complete.cases(index.df),]

  # Match codes to data frame
  if (type=="treat") {
    var <- "arm"
  } else if (type=="time") {
    var <- "fup"
  }

  if (type %in% c("treat", "time")) {
    studyarm.index <- paste(index.df[,1], index.df[,2], sep="_")
    studyarm.df <- paste(add.df$study, add.df[[var]], sep="_")
    add.df[[type]] <- index.df[,3][match(studyarm.df, studyarm.index)]
  } else if (type=="y") {
    studyarm.index <- paste(index.df[,1], index.df[,2], index.df[,3], sep="_")
    studyarm.df <- paste(add.df$study, add.df$arm, add.df$fup, sep="_")
    add.df[[type]] <- index.df[,4][match(studyarm.df, studyarm.index)]
  }

  return(add.df)
}






#' Extracts fitted values or deviance contributions into a data.frame with indices
#' @noRd
get.theta.dev <- function(mbnma, param="theta") {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertChoice(param, choices = c("dev", "resdev", "theta"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  dev.df <- as.data.frame(
    mbnma$BUGSoutput$summary[grep(paste0("^", param, "\\["),
                                  rownames(mbnma$BUGSoutput$summary)),]
  )

  if (mbnma$type=="time") {
    # Takes the study, arm and follow-up measurement identifier for each residual deviance point
    id <- matrix(unlist(strsplit(
      gsub(
        "([a-z]+\\[)([0-9]+)(,)([0-9]+)(,)([0-9]+)(\\])",
        "\\2.\\4.\\6", rownames(dev.df)),
      split="\\.")),
      byrow=TRUE,
      ncol=3
    )

    dev.df$fup <- as.numeric(id[,3])

  } else if (mbnma$type=="dose") {
    # Takes the study, arm and follow-up measurement identifier for each residual deviance point
    id <- matrix(unlist(strsplit(
      gsub(
        "([a-z]+\\[)([0-9]+)(,)([0-9]+)(\\])",
        "\\2.\\4", rownames(dev.df)),
      split="\\.")),
      byrow=TRUE,
      ncol=2
    )
  }

  dev.df$study <- as.numeric(id[,1])
  dev.df$arm <- as.numeric(id[,2])

  return(dev.df)
}







#' @describeIn mb.nodesplit Plot outputs from nodesplit models
#'
#' @param x An object of `class("mb.nodesplit")`
#' @param plot.type A character string that can take the value of `"forest"` to plot
#' only forest plots, `"density"` to plot only density plots, or left as `NULL` (the
#' default) to plot both types of plot.
#' @param params A character vector corresponding to a time-course parameter(s) for which to plot results.
#' If left as `NULL` (the default), nodes-split results for all time-course parameters will be plotted.
#'
#' @details The S3 method `plot()` on an `mb.nodesplit` object generates either
#' forest plots of posterior medians and 95\\% credible intervals, or density plots
#' of posterior densities for direct and indirect evidence.
#'
#' @return Plots the desired graph(s) and returns an object (or list of objects if
#' `plot.type=NULL`) of `class(c("gg", "ggplot"))`, which can be edited using `ggplot` commands.
#'
#' @export
plot.mb.nodesplit <- function(x, plot.type=NULL, params=NULL, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mb.nodesplit", add=argcheck)
  checkmate::assertChoice(plot.type, choices = c("density", "forest"), null.ok=TRUE, add=argcheck)
  checkmate::assertCharacter(params, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check if params are within nodesplit
  if (!is.null(params)) {
    for (i in seq_along(params)) {
      if (!(params[i] %in% names(x[[1]]))) {
        stop(paste0(params[i], " is not a time-course parameter in `x` for which node-split results are available"))
      }
    }
  } else if (is.null(params)) {
    params <- names(x[[1]])
  }


  if (is.null(plot.type)) {
    plot.type <- c("density", "forest")
  }

  forestdata <- x[[1]]$forest.plot$data[0,]
  densitydata <- x[[1]]$density.plot$data[0,]
  forestfacet <- vector()
  forestparam <- vector()
  densityfacet <- vector()
  densityparam <- vector()
  for (i in seq_along(x)) {
    for (k in seq_along(params)) {
      comp <- paste(x[[i]][[k]]$comparison, collapse=" vs ")
      temp <- x[[i]][[k]]$forest.plot$data
      forestfacet <- append(forestfacet, rep(comp, nrow(temp)))
      forestparam <- append(forestparam, rep(params[k], nrow(temp)))
      forestdata <- rbind(forestdata, temp)

      temp <- x[[i]][[k]]$density.plot$data
      densityfacet <- append(densityfacet, rep(comp, nrow(temp)))
      densityparam <- append(densityparam, rep(params[k], nrow(temp)))
      densitydata <- rbind(densitydata, temp)
    }
  }
  forestdata$comp <- forestfacet
  densitydata$comp <- densityfacet
  forestdata$param <- forestparam
  densitydata$param <- densityparam

  plotlist <- list()
  for (k in seq_along(params)) {
    if ("forest" %in% plot.type) {
      plotdata <- forestdata[forestdata$param==params[k],]
      forest <-
        ggplot2::ggplot(data=plotdata, ggplot2::aes(x=plotdata$source, y=plotdata$med,
                                                    ymin=plotdata$l95, ymax=plotdata$u95)) +
        ggplot2::geom_pointrange() +
        ggplot2::coord_flip() +  # flip coordinates (puts labels on y axis)
        ggplot2::xlab("") + ggplot2::ylab("Treatment effect (95% CrI)") +
        ggplot2::theme(axis.text = ggplot2::element_text(size=10),
                       axis.title = ggplot2::element_text(size=12),
                       title=ggplot2::element_text(size=18)) +
        ggplot2::theme(plot.margin=ggplot2::unit(c(1,1,1,1),"cm")) +
        ggplot2::facet_wrap(~factor(plotdata$comp)) +
        ggplot2::ggtitle(paste0("Forest plot of node-split for ", params[k])) +
        theme_mbnma()

      graphics::plot(forest, ...)
      plotlist[[length(plotlist)+1]] <- forest
    }
    if ("density" %in% plot.type) {
      plotdata <- densitydata[densitydata$param==params[k],]

      density <- ggplot2::ggplot(plotdata, ggplot2::aes(x=plotdata$value,
                                                        linetype=plotdata$Estimate, fill=plotdata$Estimate)) +
        ggplot2::geom_density(alpha=0.2) +
        ggplot2::xlab("Treatment effect (95% CrI)") +
        ggplot2::ylab("Posterior density") +
        ggplot2::theme(strip.text.x = ggplot2::element_text(size=12)) +
        ggplot2::theme(axis.text = ggplot2::element_text(size=12),
                       axis.title = ggplot2::element_text(size=14)) +
        ggplot2::facet_wrap(~factor(plotdata$comp)) +
        ggplot2::ggtitle(paste0("Posterior densities of node-split for ", params[k])) +
        ggplot2::guides(fill=ggplot2::guide_legend((title="Evidence Source")),
                        linetype=ggplot2::guide_legend((title="Evidence Source"))) +
        theme_mbnma()

      graphics::plot(density, ...)
      plotlist[[length(plotlist)+1]] <- density
    }
  }

  return(invisible(plotlist))
}
