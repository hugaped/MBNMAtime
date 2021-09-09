# Graphical functions
# Author: Hugo Pedder
# Date created: 2018-09-10




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




#' Predict "lumped" NMA predictions for plotting
#'
#' @param pred An object of `class("mb.predict")`
#' @param incl.range A numeric vector of length 2 representing a range between which time-points should be synthesised
#' @inheritParams plot.mb.predict
#' @inheritParams mb.run
#'
#' @noRd
overlay.nma <- function(pred, incl.range, method="common", link="identity", ...) {

  # Declare global variable
  cor <- NULL
  dif <- NULL

  network <- pred$network
  nmanet <- network$data.ab[network$data.ab$time>=incl.range[1] & network$data.ab$time<=incl.range[2],]

  # Take the follow-up closes to mean(incl.range) if multiple fups are within the range
  nmanet <- nmanet %>% dplyr::group_by(studyID) %>%
    dplyr::mutate(dif=time-mean(incl.range)) %>%
    dplyr::arrange(dif) %>%
    dplyr::group_by(studyID, arm) %>%
    dplyr::slice_head()

  if (!1 %in% nmanet$treatment) {
    stop("Network reference treatment (treatment=1) must be evaluated at a time-point within 'incl.range' to allow overlay.nma")
  }

  # Check network connectedness
  comparisons <- mb.comparisons(nmanet)
  nodes <- unique(sort(nmanet$treatment))
  nodes <- network$treatments[nodes]
  g <- igraph::graph.empty()
  g <- g + igraph::vertex(nodes)
  ed <- t(matrix(c(comparisons[["t1"]], comparisons[["t2"]]), ncol = 2))
  ed <- factor(as.vector(ed), labels=nodes)
  edges <- igraph::edges(ed, weight = comparisons[["nr"]], arrow.mode=0)
  g <- g + edges

  discon <- check.network(g)

  # Drop treatments disconnected to network reference
  if (length(discon)>0) {
    message(paste("The following treatments are disconnected from the network reference for studies in 'incl.range' and will be excluded from overlay.nma:",
                  paste(discon, collapse = "\n"), sep="\n"))

    drops <- which(network$treatments %in% drops)
    nmanet <- nmanet[!nmanet$treatment %in% drops,]

    nodes <- unique(sort(nmanet$treatment))
    nodes <- network$treatments[nodes]
  }

  # Recode treatments from 1-X
  nmanet$treatment <- as.numeric(factor(nmanet$treatment))

  # Run model (incl write priors)
  nma <- nma.run(data.ab=nmanet, method=method, link=link, ...)

  if (any(nma$BUGSoutput$summary[,"Rhat"]>1.02)) {
    warn <- rownames(nma$BUGSoutput$summary)[which(nma$BUGSoutput$summary[,"Rhat"]>1.02)]
    warning(crayon::bold(crayon::red(paste0("Rhat values for the following parameters are >1.02 - suggests problems with NMA model convergence:\n",
                  paste(warn, collapse="\n")))))
  }

  # Predict NMA results at specific time point
  timedif <- abs(pred$times - mean(incl.range))
  if (!any(timedif < (mean(incl.range) - incl.range[1]))) {
    stop("No time-point in predicted values is within 'incl.range'.\nMust specify at least one of 'times' in 'predict()' to be within 'incl.range'")
  }
  timeindex <- order(timedif)[1]
  predref <- pred$pred.mat[[1]][[timeindex]]

  if (method=="common") {
    predtrt <- sample(predref, size=nma$BUGSoutput$n.sims, replace=TRUE) + nma$BUGSoutput$sims.list$d[,-1]
  } else if (method=="random") {

    if (is.vector(nma$BUGSoutput$sims.list$d[,-1])) {
      cols <- 1
    } else {
      cols <- ncol(nma$BUGSoutput$sims.list$d[,-1])
    }

    mat <- array(dim=c(nma$BUGSoutput$n.sims,
                       cols,
                       2)
                 )
    mat[,,1] <- nma$BUGSoutput$sims.list$d[,-1]
    mat[,,2] <- nma$BUGSoutput$sims.list$sd
    predtrt <- apply(mat, MARGIN=c(1,2), FUN=function(x) stats::rnorm(1, x[1], x[2]))
    predtrt <- sample(predref, size=nma$BUGSoutput$n.sims, replace=TRUE) + predtrt
  }
  if (is.vector(predtrt)) {
    predtrt <- matrix(predtrt, ncol=1)
  }
  colnames(predtrt) <- nodes[-1]

  pred.df <- as.data.frame(t(apply(predtrt, MARGIN=2, FUN=function(x){stats::quantile(x, probs = c(0.025,0.5,0.975))})))
  pred.df$time <- pred$times[timeindex]
  pred.df$treat <- rownames(pred.df)

  # Select only treatments included in pred
  pred.df <- pred.df[pred.df$treat %in% names(pred$summary),]

  if (method=="random") {
    sd <- nma$BUGSoutput$summary[rownames(nma$BUGSoutput$summary)=="sd",]
    sd <- round(sd, 3)
    sd <- paste(sd[5], " (", sd[3], ", ", sd[7], ")", sep="")
  } else {
    sd <- NULL
  }

  return(list(pred.df=pred.df,
              totresdev=round(nma$BUGSoutput$median$totresdev,1), ndat=nrow(nmanet),
              dic=round(nma$BUGSoutput$median$totresdev+nma$BUGSoutput$pD,1),
              sd=sd))

  # DO BITS FOR RANDOM AND ADD ANNOTATIONS TO MODEL WITH RESIDUAL DEVIANCE AND N DATA POINTS AND SD (95%CRI)
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
#' @param plotby A character object that can take either `"arm"` to indicate that raw responses
#' should be plotted separately for each study arm, or `"rel"` to indicate that the relative
#' effects within each study should be plotted. In this way the time-course of both the absolute
#' effects and the relative effects can be examined.
#' @param ... Arguments to be sent to `ggplot()`
#' @inheritParams plot.mb.network
#' @inheritParams mb.run
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
#' \donttest{
#' # Make network
#' goutnet <- mb.network(goutSUA_CFB)
#'
#' # Use timeplot to plot responses grouped by treatment
#' timeplot(goutnet)
#'
#' # Use timeplot ot plot resposes grouped by class
#' timeplot(goutnet, level="class")
#'
#' # Plot matrix of relative effects
#' timeplot(goutnet, level="class", plotby="rel")
#'
#' # Plot using Standardised Mean Differences
#' copdnet <- mb.network(copd)
#' timeplot(copdnet, plotby="rel", link="smd")
#'
#' }
#'
#' @export
timeplot <- function(network, level="treatment", plotby="arm", link="identity", ...) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mb.network", add=argcheck)
  checkmate::assertChoice(level, choices = c("treatment", "class"), add=argcheck)
  checkmate::assertChoice(plotby, choices = c("arm", "rel"), add=argcheck)
  checkmate::assertChoice(link, choices = c("identity", "smd", "rom"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Define global variables
  poolvar <- NULL
  Rx.Name.x <- NULL
  Rx.Name.y <- NULL

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
    g <- g + ggplot2::ylab("Response")

  } else if (plotby=="rel") {

    diffs <- plotdata %>%
      dplyr::mutate(Rx.Name = factor(network$treatments[.data$treatment], levels=network$treatments))

    if (link=="identity") {
      diffs <- diffs %>%
        dplyr::inner_join(diffs, by=c("studyID", "time")) %>%
        dplyr::filter(.data$treatment.x < .data$treatment.y) %>%
        dplyr::mutate(pairDiff = .data$y.y - .data$y.x)
    } else if (link=="rom") {
      diffs <- diffs %>%
        dplyr::inner_join(diffs, by=c("studyID", "time")) %>%
        dplyr::filter(.data$treatment.x < .data$treatment.y) %>%
        dplyr::mutate(pairDiff = log(.data$y.y/.data$y.x))
    } else if (link=="smd") {
      if (!"n" %in% names(diffs)) {
        stop("'n' must be included in dataset if using link='smd'")
      }

      diffs <- diffs %>%
        dplyr::inner_join(diffs, by=c("studyID", "time")) %>%
        dplyr::filter(.data$treatment.x < .data$treatment.y) %>%
        dplyr::mutate(
          var.y = (.data$se.y * (.data$n.y)^0.5)^2,
          var.x = (.data$se.x * (.data$n.x)^0.5)^2,
          poolvar = ((.data$var.y * (.data$n.y-1)) + (.data$var.x * (.data$n.x-1))) / (.data$n.y + .data$n.x-2),
          pairDiff = (.data$y.y - .data$y.x) / (poolvar^0.5)
          )
    }


    diffs <- diffs %>%
      dplyr::bind_rows(diffs %>%
                         dplyr::group_by(.data$studyID, .data$arm.x, .data$arm.y) %>%
                         dplyr::slice(1) %>%
                         dplyr::mutate(time=0, pairDiff=0))

    if (level=="class") {

      diffs <- diffs %>%
        dplyr::mutate(Rx.Name.x = network$classkey$class[match(Rx.Name.x, network$classkey$treatment)],
                      Rx.Name.y = network$classkey$class[match(Rx.Name.y, network$classkey$treatment)]
                      )
    }

    # g <- ggplot2::ggplot(data=diffs, ggplot2::aes(x=.data$time, y=.data$pairDiff, group=.data$studyID)) +
    #   ggplot2::geom_line() +
    #   ggplot2::geom_point() +
    #   ggplot2::facet_grid(rows=ggplot2::vars(.data$Rx.Name.y), cols=ggplot2::vars(.data$Rx.Name.x))

    g <- ggplot2::ggplot(data=diffs, ggplot2::aes(x=.data$time, y=.data$pairDiff,
                                                  # group=.data$studyID),
                                                  group=paste(.data$studyID, .data$arm.x, .data$arm.y, sep="_")),
                         ...) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_grid(rows=ggplot2::vars(.data$Rx.Name.y), cols=ggplot2::vars(.data$Rx.Name.x))

    if (link=="identity") {
      g <- g + ggplot2::ylab("Response")
    } else if (link=="rom") {
      g <- g + ggplot2::ylab("Log-Ratio of Means")
    } else if (link=="smd") {
      g <- g + ggplot2::ylab("Standardised Mean Difference")
    }
  }

  g <- g + ggplot2::xlab("Time") +
    theme_mbnma()

  graphics::plot(g)
  return(invisible(g))

}






#' Plot deviance contributions from an MBNMA model
#'
#' @param mbnma An S3 object of class `"mbnma"` generated by running
#' a time-course MBNMA model
#' @param dev.type Deviances to plot - can be either residual deviances (`"resdev"`) or deviances (`"dev"`, the default)
#' @param plot.type Deviances can be plotted either as scatter points (`"scatter"` - using
#' `geom_point()`) or as boxplots (`"box"`, the default)
#' @param xaxis A character object that indicates whether deviance contributions should be plotted
#' by time (`"time"`) or by follow-up count (`"fup"`)
#' @param facet A boolean object that indicates whether or not to facet by treatment
#' @param n.iter The number of iterations to update the model whilst monitoring additional parameters (if necessary).
#' Must be a positive integer. Default is the value used in `mbnma`.
#' @param n.thin The thinning rate. Must be a positive integer. Default is the value used in `mbnma`.
#' @param ... Arguments to be sent to `ggplot2::ggplot()`
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
#' alognet <- mb.network(alog_pcfb)
#'
#' # Run MBNMA
#' mbnma <- mb.run(alognet, fun=tpoly(degree=2), intercept=FALSE)
#'
#' # Plot residual deviance contributions in a scatterplot
#' devplot(mbnma)
#'
#' # Plot deviance contributions in boxplots at each follow-up measurement
#' # Monitor for 500 additional iterations
#' devplot(mbnma, dev.type="dev", plot.type="box", xaxis="fup", n.iter=500)
#' }
#' @export
devplot <- function(mbnma, dev.type="dev", plot.type="box",
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
  if (mbnma$model.arg$covar %in% c("AR1", "CS")) {
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
    dev.df$fup <- apply(dev.df, MARGIN=1, FUN=function(x) mbnma$model.arg$jagsdata$time[x[1], x[3]])
    xlab <- "Time"
  } else if (xaxis=="fup") {
    xlab <- "Follow-up count"
  }

  if (facet==TRUE) {
    dev.df$facet <- apply(dev.df, MARGIN=1, FUN=function(x) mbnma$model.arg$jagsdata$treat[x[1], x[2]])
  }

  # Plots the residual deviances over time grouped by study and arm
  if (plot.type=="scatter") {
    g <- ggplot2::ggplot(dev.df, ggplot2::aes(x=fup, y=mean), group=(paste(study, arm, sep="_")), ...) +
      ggplot2::geom_point()
  } else if (plot.type=="box") {
    g <- ggplot2::ggplot(dev.df, ggplot2::aes(x=factor(fup), y=mean), ...) +
      ggplot2::geom_boxplot()
  }

  # Add facets
  if (facet==TRUE) {
    g <- g + ggplot2::facet_wrap(~factor(facet, labels=mbnma$network$treatments), scales="free_x")
  }

  # Add axis labels
  g <- g + ggplot2::xlab(xlab) +
    ggplot2::ylab("Posterior mean") +
    ggplot2::geom_hline(yintercept = 0, lty="dashed") +
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
#' @param ... Arguments to be sent to `ggplot2::ggplot()`
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
#' painnet <- mb.network(osteopain)
#'
#' # Run MBNMA
#' mbnma <- mb.run(painnet,
#'   fun=temax(pool.emax="rel", method.emax="common",
#'     pool.et50="abs", method.et50="random"))
#'
#' # Plot fitted values from the model
#' # Monitor fitted values for 500 additional iterations
#' fitplot(mbnma, n.iter=500)
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
                       ggplot2::aes(x=time, y=mean, group=paste(study, arm, sep="_")), ...) +
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
  #index.df <- reshape2::melt(mbnma$model$data()[[type]])
  index.df <- reshape2::melt(mbnma$model.arg$jagsdata[[type]])
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









#' Plot illustrative time-course functions (CURRENTLY NOT FUNCTIONAL FOR ALL TIME_COURSE FUNCTIONS)
#'
#' Can be used to plot potential time-course functions to identify if they may
#' be a suitable fit for the data.
#'
#' @param x An object of `class("timefun")`
#' @inheritParams tuser
#'
#' @return An illustrative plot of the time-course function with the parameters specified
#' in `beta.1`, `beta.2`, `beta.3` and `beta.4`
#'
#' @noRd
plot.timefun <- function(x=tpoly(degree=1), beta.1=0, beta.2=0,
                         beta.3=0, beta.4=0) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "timefun", add=argcheck)
  checkmate::assertNumeric(beta.1, add=argcheck)
  checkmate::assertNumeric(beta.2, add=argcheck)
  checkmate::assertNumeric(beta.3, add=argcheck)
  checkmate::assertNumeric(beta.4, add=argcheck)
  checkmate::reportAssertions(argcheck)

  time <- seq(0,100)

  funstr <- x$jags
  #funstr <- gsub("beta", "tempbeta", funstr)
  funstr <- gsub("\\[i\\,k\\]", "", funstr)
  funstr <- gsub("\\[i\\,m\\]", "", funstr)

  if (any(c("rcs", "ns", "bs") %in% x$name)) {
    spline <- genspline(time, spline=x$name, knots=x$knots)
    funstr <- gsub("spline\\[i\\,m", "spline[", funstr)
  }

  if (any(c("user", "fpoly") %in% x$name)) {
    stop("Plotting for 'tfpoly()' and 'tuser()' not currently supported")
  }

  nparam <- x$nparam

  df <- data.frame("y"=NA, "time"=NA)
  y <- eval(parse(text=funstr))
  df <- rbind(df, stats::setNames(cbind(y,time), names(df)))
  df <- df[-1,]

  vals <- vector()
  for (i in 1:nparam) {
    vals <- append(vals, get(paste0("beta.", i)))
  }

  g <- ggplot2::ggplot(df, ggplot2::aes(x=time, y=y)) +
    ggplot2::geom_line(size=1) +
    ggplot2::labs(subtitle=paste(paste0("beta.", 1:nparam, " = ", vals), collapse="    ")) +
    ggplot2::xlab("Time") +
    theme_mbnma()

  graphics::plot(g)
  return(invisible(g))
}






plot.invisible <- function(...){
  ff <- tempfile()
  grDevices::png(filename=ff)
  res <- graphics::plot(...)
  grDevices::dev.off()
  unlink(ff)
  return(res)
}




#' MBNMA ggplot2 theme style
#' @noRd
theme_mbnma <- function(...) {
  ggplot2::theme_bw(...) +
    ggplot2::theme(
      # change stuff here
      panel.background  = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill="white", colour=NA),
      legend.background = ggplot2::element_rect(fill="transparent", colour=NA),
      legend.key = ggplot2::element_rect(fill="transparent", colour=NA),

      # From multinma
      #panel.border = ggplot2::element_rect(colour = "grey70", fill = NA),
      panel.grid.major = ggplot2::element_line(colour = "grey95"),
      panel.grid.minor = ggplot2::element_line(colour = "grey95"),
      strip.background = ggplot2::element_rect(colour = "black",
                                               fill = "lightsteelblue1"),
      strip.text = ggplot2::element_text(colour = "black")
    )
}
