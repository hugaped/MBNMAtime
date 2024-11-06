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
#' @return A numeric vector of rescaled values
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
#' FUNCTION IS DEPRECATED
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
#' @param timebins A numeric vector defining the boundaries of the time bins within which to perform
#' a standard NMA. Length must be >=2. See details.
#' @inheritParams plot.mb.predict
#' @inheritParams mb.run
#' @inheritParams mb.predict
#'
#' @details `timebins` indicate regions of the data (defined as "time bins") over which it may be reasonable to "lump" different
#' follow-up times from different studies together and assume a standard NMA model.
#'
#' @noRd
overlay.nma <- function(pred, timebins, method="common", link="identity", lim="cred",
                        plottype="pred", ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  if (plottype=="pred") {
    checkmate::assertClass(pred, c("mb.predict"), add=argcheck)
  }
  checkmate::assertNumeric(timebins, lower=0, null.ok=FALSE, add=argcheck)
  checkmate::assertChoice(lim, choices=c("cred", "pred"), add=argcheck)
  checkmate::assertChoice(plottype, choices=c("rel", "pred"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Declare global variable
  cor <- NULL
  dif <- NULL

  network <- pred$network

  outlist <- list()

  # Loop across timebins
  for (bin in 2:length(timebins)) {

    skip <- FALSE

    incl.range <- c(timebins[bin-1], timebins[bin])
    str.incl.range <- paste(c("time=", incl.range[1], " and time=", incl.range[2]), collapse="")

    if (plottype=="pred") {
      # Predict NMA results at specific time point
      timedif <- abs(pred$times - mean(incl.range))
      if (!any(timedif < (mean(incl.range) - incl.range[1]))) {
        warning(paste(c("No time-point in predicted values is between ", str.incl.range,
                        ". Cannot overlay.nma for this time bin."), collapse=""))

        skip <- TRUE
      }
      timeindex <- order(timedif)[1]
      predref <- pred$pred.mat[[1]][[timeindex]]

    } else if (plottype=="rel") {
      timeindex <- 1
      predref <- 0
    }

    if (skip==FALSE) {
      nmanet <- network$data.ab[network$data.ab$time>incl.range[1] & network$data.ab$time<=incl.range[2],]

      # Take the follow-up closes to mean(incl.range) if multiple fups are within the range
      nmanet <- nmanet %>% dplyr::group_by(studyID) %>%
        #dplyr::mutate(dif=time-mean(incl.range)) %>%
        dplyr::mutate(dif=time-ifelse(plottype=="pred", pred$times[timeindex], mean(incl.range))) %>%
        dplyr::arrange(dif) %>%
        dplyr::group_by(studyID, arm) %>%
        dplyr::slice_head()

      if (!1 %in% nmanet$treatment) {
        message(paste(c("overlay.nma not possible between ", str.incl.range,
                        ". Data for network reference (treatment=1) not available in this time bin"),
                      collapse = ""))

        skip <- TRUE
      }
    }

    if (skip==FALSE) {

      # Store max range of data over which NMA is performed
      if (plottype=="pred") {
        dat.range <- range(c(nmanet$time, pred$times[timeindex]))
      } else if (plottype=="rel") {
        dat.range <- range(c(nmanet$time))
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

      discon <- suppressWarnings(check.network(g))

      # Drop treatments disconnected to network reference
      if (length(discon)>0) {
        # message(paste("The following treatments are disconnected from the network reference for studies in 'incl.range' and will be excluded from overlay.nma:",
        #               paste(discon, collapse = "\n"), sep="\n"))

        drops <- which(network$treatments %in% discon)
        nmanet <- nmanet[!nmanet$treatment %in% drops,]

        nodes <- unique(sort(nmanet$treatment))
        nodes <- network$treatments[nodes]
      }

      # Recode treatments from 1-X
      nmanet$treatment <- as.numeric(factor(nmanet$treatment))

      # Print message
      message(paste0("Running overlay.nma for ", str.incl.range))

      # Run model (incl write priors)
      nma <- nma.run(data.ab=nmanet, method=method, link=link, pD=TRUE,
                     ...)

      if (method=="common" | "cred" %in% lim) {
        predtrt <- sample(predref, size=nma$BUGSoutput$n.sims, replace=TRUE) + nma$BUGSoutput$sims.list$d[,-1]
      } else if (method=="random" & "pred" %in% lim) {

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
      if (plottype=="pred") {
        pred.df <- pred.df[pred.df$treat %in% names(pred$summary),]
      }

      # Add range of data to pred.df
      pred.df$tmin <- dat.range[1]
      pred.df$tmax <- dat.range[2]

      # Create width tolerance for plotting
      if (dat.range[1]==dat.range[2]) {
        scale <- range(timebins)
        width <- (scale[2] - scale[1])/200

        pred.df$tmin <- pred.df$tmin - width
        pred.df$tmax <- pred.df$tmax + width
      }

      if (method=="random") {
        sd <- nma$BUGSoutput$summary[rownames(nma$BUGSoutput$summary)=="sd",]
        sd <- round(sd, 3)
        sd <- paste(sd[5], " (", sd[3], ", ", sd[7], ")", sep="")
      } else {
        sd <- NULL
      }

      listnam <- paste0("bin", bin-1)
      outlist[[listnam]] <- list(pred.df=pred.df,
                                 nma.summary=nma$BUGSoutput$summary,
                                 totresdev=round(nma$BUGSoutput$median$totresdev,1), ndat=nrow(nmanet),
                                 dic=round(nma$BUGSoutput$median$totresdev+nma$BUGSoutput$pD,1),
                                 sd=sd,
                                 incl.range=dat.range)
    }
  }

  return(outlist)
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
#' should be plotted separately for each study arm, or `"rel"` to indicate that the within-study
#' relative effects/treatment differences should be plotted. In this way the time-course of both the absolute
#' effects and the relative effects can be examined.
#' @param ... Arguments to be sent to `ggplot()`
#' @inheritParams plot.mb.network
#' @inheritParams mb.run
#'
#' @return The function returns an object of `class(c("gg", "ggplot")`. Characteristics
#' of the object can therefore be amended as with other plots generated by `ggplot()`.
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
timeplot <- function(network, level="treatment",
                     plotby="arm", link="identity",
                     ...) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mb.network", add=argcheck)
  checkmate::assertChoice(level, choices = c("treatment", "class"), add=argcheck)
  checkmate::assertChoice(plotby, choices = c("arm", "rel"), add=argcheck)
  checkmate::assertChoice(link, choices = c("identity", "smd", "log"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Define global variables
  poolvar <- NULL
  Rx.Name.x <- NULL
  Rx.Name.y <- NULL

  if (level=="class" & !("classes" %in% names(network))) {
    stop("`level` cannot be set to class if there is no class variable included in the dataset/network")
  }

  plotdata <- network$data.ab
  cfb <- network$cfb

  if (plotby=="arm") {

    studies <- unique(plotdata$studyID)
    for (i in seq_along(studies)) {
      base.dat <- subset(plotdata, fupcount==1 & studyID==studies[i])
      if (cfb[i]==TRUE & any(base.dat$time!=0)) {
        base.dat$y <- rep(0, nrow(base.dat))
        base.dat$se <- rep(0, nrow(base.dat))
        base.dat$time <- rep(0, nrow(base.dat))

        plotdata <- rbind(base.dat, plotdata)
      }
    }
    plotdata <- dplyr::arrange(plotdata, studyID, time, treatment)

    # Do not run this function with pylr loaded!!
    plotdata <- plotdata %>%
      dplyr::group_by(studyID, time) %>%
      dplyr::mutate(arm = sequence(dplyr::n()))

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
        dplyr::inner_join(diffs, by=c("studyID", "time"), relationship="many-to-many") %>%
        dplyr::filter(.data$treatment.x < .data$treatment.y) %>%
        dplyr::mutate(pairDiff = .data$y.y - .data$y.x)
    } else if (link=="log") {
      diffs <- diffs %>%
        dplyr::inner_join(diffs, by=c("studyID", "time"), relationship="many-to-many") %>%
        dplyr::filter(.data$treatment.x < .data$treatment.y) %>%
        dplyr::mutate(pairDiff = log(.data$y.y/.data$y.x))
    } else if (link=="smd") {
      if (!"n" %in% names(diffs)) {
        stop("'n' must be included in dataset if using link='smd'")
      }

      diffs <- diffs %>%
        dplyr::inner_join(diffs, by=c("studyID", "time"), relationship="many-to-many") %>%
        dplyr::filter(.data$treatment.x < .data$treatment.y) %>%
        dplyr::mutate(
          var.y = (.data$se.y * (.data$n.y)^0.5)^2,
          var.x = (.data$se.x * (.data$n.x)^0.5)^2,
          poolvar = ((.data$var.y * (.data$n.y-1)) + (.data$var.x * (.data$n.x-1))) / (.data$n.y + .data$n.x-2),
          pairDiff = (.data$y.y - .data$y.x) / (poolvar^0.5)
          )
    }

    # Handle cfb data
    studies <- unique(diffs$studyID)
    diffs$bl <- FALSE
    for (i in seq_along(cfb)) {
      if (cfb[i]==FALSE) {
        diffs$bl[diffs$studyID==studies[i] & diffs$fupcount.x==1 & diffs$fupcount.y==1] <- TRUE
      }
    }

    diffs <- diffs %>% dplyr::bind_rows(diffs %>% dplyr::group_by(.data$studyID,
                                                                  .data$arm.x, .data$arm.y) %>%
                                          dplyr::slice(1) %>%
                                          dplyr::mutate(time=dplyr::case_when(bl==0 ~ 0),
                                                        pairDiff=dplyr::case_when(bl==0 ~ 0),
                                                        bl=dplyr::case_when(bl==0 ~ 1)))

    diffs <- diffs[!is.na(diffs$time),]

    if (level=="class") {

      diffs <- diffs %>%
        dplyr::mutate(Rx.Name.x = network$classkey$class[match(Rx.Name.x, network$classkey$treatment)],
                      Rx.Name.y = network$classkey$class[match(Rx.Name.y, network$classkey$treatment)]
                      )
    }

    g <- ggplot2::ggplot(data=diffs, ggplot2::aes(x=.data$time, y=.data$pairDiff,
                                                  # group=.data$studyID),
                                                  group=paste(.data$studyID, .data$arm.x, .data$arm.y, sep="_")),
                         ...) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_grid(rows=ggplot2::vars(.data$Rx.Name.y), cols=ggplot2::vars(.data$Rx.Name.x))

    if (link=="identity") {
      g <- g + ggplot2::ylab("Response")
    } else if (link=="log") {
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




#' Plot relative effects from NMAs performed at multiple time-bins
#'
#' @param legend `TRUE`/`FALSE` to indicate whether a legend should be plotted.
#' @inheritParams mb.run
#' @inheritParams plot.mb.predict
#' @inheritParams predict.mbnma
#'
#' @return Plots treatment effects from NMAs performed within discrete time bins. The
#' object returned is a list containing the plot and a sublist of NMA results and
#' predictions from each time bin specified in `overlay.nma`.
#'
#' @details
#' Performs several standard NMAs at different time "bins", time periods within
#' which treatment effects are assumed to be constant over time. Separate NMAs
#' are then performed within each time bin on data points from studies that fall
#' within the time bin (only a single follow-up time is taken from each study
#' to avoid double counting).
#'
#' Note that the wider the time bin boundaries specified by the user, the
#' larger the potential range of included follow-up times and this can
#' introduce heterogeneity or inconsistency.
#'
#' Results are plotted versus the network reference and are plotted on the
#' specified link scale. Each time bin window is marked on the plot by
#' vertical dashed lines. The NMA estimates within each time bin are plotted
#' as a horizontal solid black line (the posterior median) with a shaded region
#' indicating the 95% credible interval (prediction intervals can instead
#' be plotted). The width of these shaded regions is equal to the range of study
#' time-points included in the NMA performed within that timebin, which
#' may therefore be more narrow than the time bin specified in the `binplot()`
#' command due to the follow-up times at which data is available in included
#' studies.
#'
#' @inheritSection plot.mb.predict Overlaying NMA results
#'
#' @examples
#' \donttest{
#' # Create an mb.network object from a dataset
#' alognet <- mb.network(alog_pcfb)
#'
#' # Plot relative effects from NMAs calculated for a single time-bins
#' # Do not plot time-bin boundaries
#' binplot(alognet, overlay.nma=c(0,5), plot.bins=FALSE)
#'
#' # Plot relative effects from NMAs at multiple time-bins
#' # With random treatment effects
#' binplot(alognet, overlay.nma=c(5,10,15,20),
#'   method="random")
#' }
#'
#' @export
binplot <- function(network, overlay.nma=c(0, stats::quantile(network$data.ab$time)),
                        method="common", link="identity", lim="cred", plot.bins=TRUE,
                    legend=TRUE, ...) {

  checkmate::assertNumeric(overlay.nma, min.len = 2)

  # Ensure timebins are unique
  overlay.nma <- unique(overlay.nma)

  # Create pred to match overlay.nma
  pred <- list()
  pred$times <- seq(0,max(overlay.nma), length.out=100)
  pred$pred.mat <- list(matrix(0, ncol=100, nrow=1))
  pred$network <- network

  class(pred) <- "mb.predict"

  nma <- overlay.nma(pred, timebins=overlay.nma, method=method, link=link, lim=lim,
                     plottype="rel", ...)

  if (length(nma)==0) {
    stop("No NMA can be performed between time points specified in overlay.nma")
  }

  # Bind data frames
  plot.df <- nma[[1]]$pred.df[0,]
  for (i in seq_along(nma)) {
    plot.df <- rbind(plot.df, nma[[i]]$pred.df)
  }

  g <- ggplot2::ggplot(data=plot.df, ggplot2::aes(x=tmin, xend=tmax, ymin=`2.5%`, ymax=`97.5%`, y=`50%`)) +
    ggplot2::geom_rect(ggplot2::aes(ymin=`2.5%`, ymax=`97.5%`, xmin=tmin, xmax=tmax,
                                           fill="95% Interval"),
                              alpha=1) +
    ggplot2::geom_segment(ggplot2::aes(y=`50%`, yend=`50%`, x=tmin, xend=tmax, color="Posterior median"),
                          linewidth=0.8) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=0), linetype="dotted")

  colorvals <- c("Posterior median"="gray0")

  capt <- paste0("Results versus ", network$treatments[1],
                 "\nWidth of shaded region denotes range of follow-up times included in data for each NMA",
                 "\nHeight of shaded region denotes the 95% interval of the posterior distribution")

  if (plot.bins==TRUE) {
    capt <- paste0(capt, "\nVertical dashed lines indicate time bin boundaries")

    g <- g + ggplot2::geom_vline(xintercept=overlay.nma, linetype="dashed", alpha=0.5) +
      ggplot2::scale_x_continuous(breaks=unique(c(0, overlay.nma)))
  }

  g <- g +
    ggplot2::labs(caption=capt) +
    ggplot2::scale_fill_manual(name="", values=c("95% Interval"="lightblue"))

  g <- g + ggplot2::facet_wrap(~factor(treat)) +
    ggplot2::labs(y="Treatment effect (on link scale)", x="Time") +
    ggplot2::scale_color_manual(name="", values=colorvals) +
    theme_mbnma() +
    ggplot2::theme(legend.position="none")

  if (legend==TRUE) {
    # legend.png <- png::readPNG("man/figures/binplot_legend.PNG")
    legend.png <- png::readPNG(system.file("figures/binplot_legend.PNG", package = "MBNMAtime"))
    g2 <- grid::rasterGrob(legend.png, interpolate = TRUE)
    g <- gridExtra::arrangeGrob(g, g2, ncol=2, widths=c(10, 1))
  }

  graphics::plot(g)

  out <- list("graph"=g)
  out[["overlay.nma"]] <- nma

  return(invisible(out))

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
                    n.iter=round(mbnma$BUGSoutput$n.iter/4), n.thin=mbnma$BUGSoutput$n.thin,
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
                    n.iter=round(mbnma$BUGSoutput$n.iter/4), n.thin=mbnma$BUGSoutput$n.thin, ...) {
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

  if (any(c("ns", "bs", "ls") %in% x$name)) {
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
    ggplot2::geom_line(linewidth=1) +
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







#' Plot cumulative ranking curves from MBNMA models
#'
#' @param x An object of class `"mb.rank"` generated by `rank.mbnma()`
#' @param sucra A logical object to indicate whether Surface Under Cumulative Ranking Curve (SUCRA)
#' values should be calculated and returned as a data frame. Areas calculated
#' using trapezoid approach.
#' @param ... Arguments to be sent to `ggplot::geom_line()`
#'
#' @return Line plots showing the cumulative ranking probabilities for each agent/class for
#' the ranked dose response paramtere in `x`. The object returned is a list which contains the plot
#' (an object of `class(c("gg", "ggplot")`) and a data frame of SUCRA values
#' if `sucra = TRUE`.
#'
#' @examples
#' \donttest{
#' # Using the alogliptin data
#' network <- mb.network(alog_pcfb)
#'
#' # Estimate rankings  from an Emax dose-response MBNMA
#' emax <- mb.run(network, fun=temax())
#' ranks <- rank(emax, param=c("emax"))
#'
#' # Plot cumulative rankings for both dose-response parameters simultaneously
#' # Note that SUCRA values are also returned
#' cumrank(ranks)
#' }
#' @export
cumrank <- function(x, sucra=TRUE, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mb.rank", add=argcheck)
  checkmate::assertLogical(sucra, null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  output <- list()

  params <- x$param

  df <- data.frame()

  cum.mat <- x$cum.matrix
  treats <- colnames(cum.mat)

  melt <- reshape2::melt(cum.mat)
  melt$param <- params

  df <- rbind(df, melt)

  df$Parameter <- factor(df$param)

  g <- ggplot2::ggplot(df, ggplot2::aes(x=Var1, y=value, linetype=Parameter, colour=Parameter), ...) +
    ggplot2::geom_line(linewidth=1)

  g <- g + ggplot2::facet_wrap(~factor(Var2)) +
    ggplot2::xlab("Rank (1 = best)") +
    ggplot2::ylab("Cumulative probability") +
    ggplot2::labs(linetype="Parameter", colour="Parameter") +
    theme_mbnma()

  graphics::plot(g)

  output <- list("cumplot"=g)

  # Calculate AUC
  if (sucra==TRUE) {
    df.auc <- df %>%
      dplyr::group_by(df$Var2, df$param) %>%
      dplyr::do(data.frame(sucra=calcauc(.))) %>%
      dplyr::ungroup()

    # Normalise SUCRA values to 0,1
    df.auc$sucra <- df.auc$sucra/nrow(df.auc)

    names(df.auc) <- c("treatment", "parameter", "sucra")

    output[["sucra"]] <- df.auc

    print(df.auc)
  }

  return(invisible(output))
}
