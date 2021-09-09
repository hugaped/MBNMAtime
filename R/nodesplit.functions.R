# Functions for node-split models
# Author: Hugo Pedder
# Date created: 2018-09-10


#' Identify comparisons in loops that fulfil criteria for node-splitting
#'
#' Identify comparisons informed by both direct and indirect evidence from
#' independent sources, which therefore fulfil the criteria for testing for
#' inconsistency via node-splitting. Follows the method of van Valkenhoef \insertCite{vanvalkenhoef2016;textual}{MBNMAtime}.
#'
#' @param data A data frame containing variables `studyID` and `treatment` (as
#'   numeric codes) that indicate which treatments are used in which studies.
#'
#' @details Similar to \code{\link[gemtc]{mtc.nodesplit}} but uses a fixed
#'   reference treatment and therefore suggests fewer loops in which to test for
#'   inconsistency. Heterogeneity can also be parameterised as inconsistency and
#'   so testing for inconsistency in additional loops whilst changing the
#'   reference treatment would also be identifying heterogeneity. Depends on \code{\link[igraph]{igraph}}.
#'
#' @return A data frame of comparisons that are informed by direct and indirect
#'   evidence from independent sources. Each row of the data frame is a
#'   different treatment comparison. Numerical codes in `t1` and `t2` correspond
#'   to treatment codes.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' data <- data.frame(studyID=c(1,1,2,2,3,3,4,4,5,5,5),
#'   treatment=c(1,2,1,3,2,3,3,4,1,2,4)
#'   )
#'
#' # Identify comparisons informed by direct and indirect evidence
#' inconsistency.loops(data)
#' @export
inconsistency.loops <- function(data)
{
  # Assert checks
  checkmate::assertDataFrame(data)

  treatments <- factor(unique(data$treatment))

  data <- data %>%
    dplyr::group_by(studyID) %>%
    dplyr::mutate(design=list(as.numeric(treatment)))


  comparisons <- ref.comparisons(data)

  splits1 <- vector()
  splits2 <- vector()
  paths <- vector()
  loops <- vector()

  for (i in 1:nrow(comparisons)) {

    drops <- comparisons[-i,]

    # Alternative graph create (non-gemtc)
    g <- igraph::graph.empty()
    g <- g + igraph::vertex(levels(treatments))
    #g <- g + igraph::edges.create(new.comparisons, arrow.mode=0)

    ed <- t(matrix(c(drops[["t1"]], drops[["t2"]]), ncol = 2))
    edges <- igraph::edges(as.vector(ed), arrow.mode=0)
    g <- g + edges

    # Check whether there is still an indirect connection once direct evidence studies are removed
    if (as.logical(is.finite(igraph::shortest.paths(igraph::as.undirected(g),
                                                    comparisons[i,1], comparisons[i,2]))) == TRUE) {

      # Identify the path made by the indirect evidence
      path <- as.numeric(igraph::shortest_paths(igraph::as.undirected(g),
                                                comparisons[i,1], comparisons[i,2],
                                                weights=NA
      )[["vpath"]][[1]])

      loop <- sort(path)

      splits1 <- append(splits1, comparisons[["t1"]][i])
      splits2 <- append(splits2, comparisons[["t2"]][i])
      paths <- append(paths, paste(path, collapse="->"))
      loops <- append(loops, paste(loop, collapse="->"))

    }
  }

  splits <- data.frame("t1"=splits1, "t2"=splits2, "path"=paths, "loops"=loops, stringsAsFactors = TRUE)

  # Ensures only one comparison given per inconsistent loop
  splits <- splits[seq(dim(splits)[1],1),]
  splits <- splits[duplicated(splits[["loops"]])==FALSE, 1:3]

  if (nrow(splits)==0 | (nrow(splits)==1 & any(is.na(splits$t1)))) {
    stop("No closed loops of treatments arising from independent sources of evidence are present in the data - testing for consistency is not possible in this network")
  }

  return(splits)
}


#' Identify unique comparisons relative to study reference treatment within a
#' network
#'
#' Identify unique contrasts relative to each study reference within a network.
#' Repetitions of the same treatment comparison are grouped together.
#'
#' @param data A data frame containing variables `studyID` and `treatment` (as
#'   numeric codes) that indicate which treatments are used in which studies.
#'
#' @return A data frame of unique comparisons in which each row represents a
#'   different comparison. `t1` and `t2` indicate the treatment codes that make
#'   up the comparison. `nr` indicates the number of times the given comparison
#'   is made within the network.
#'
#'   If there is only a single observation for each study within the dataset
#'   (i.e. as for standard network meta-analysis) `nr` will represent the number
#'   of studies that compare treatments `t1` and `t2`.
#'
#'   If there are multiple observations for each study within the dataset (as in
#'   MBNMAtime) `nr` will represent the number of time points in the
#'   dataset in which treatments `t1` and `t2` are compared.
#'
#' @examples
#' data <- data.frame(studyID=c(1,1,2,2,3,3,4,4,5,5,5),
#'   treatment=c(1,2,1,3,2,3,3,4,1,2,4)
#'   )
#'
#' # Identify comparisons informed by direct and indirect evidence
#' MBNMAtime:::ref.comparisons(data)
ref.comparisons <- function(data)
{
  # Assert checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data, add=argcheck)
  checkmate::assertNames(names(data), must.include=c("studyID", "treatment"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  checkmate::assert(
    checkmate::checkFactor(data[["treatment"]]),
    checkmate::checkNumeric(data[["treatment"]])
  )

  # if (all(names(data) %in% c("studyID", "treatment") != TRUE)) {
  #   stop("data must contain variables 'studyID' and 'treatment'")
  # }
  #
  # if (!(is.factor(data[["treatment"]]) | is.numeric(data[["treatment"]]))) {
  #   stop("`treatment` must be either factor or numeric")
  # }

  data <- dplyr::arrange(data, studyID, treatment)

  t1 <- vector()
  t2 <- vector()
  for (i in seq_along(unique(data[["studyID"]]))) {
    subset <- subset(data, studyID==unique(data[["studyID"]])[i])
    for (k in 2:nrow(subset)) {
      t1 <- append(t1, subset[["treatment"]][1])
      t2 <- append(t2, subset[["treatment"]][k])
      if (is.na(subset[["treatment"]][k])) {
        stop()
      }
    }
  }

  comparisons <- data.frame(t1 = t1, t2 = t2)
  comparisons <- comparisons %>% dplyr::group_by(t1, t2) %>%
    dplyr::mutate(nr = dplyr::n())
  comparisons <- unique(comparisons)
  comparisons <- dplyr::arrange(comparisons, t1, t2)

  return(comparisons)
}


#' Identify unique comparisons within a network (identical to MBNMAdose)
#'
#' Identify unique contrasts within a network that make up all the head-to-head comparisons. Repetitions
#' of the same treatment comparison are grouped together.
#'
#' @param data A data frame containing variables `studyID` and `treatment` (as numeric codes) that
#' indicate which treatments are used in which studies.
#'
#' @return A data frame of unique comparisons in which each row represents a different comparison.
#' `t1` and `t2` indicate the treatment codes that make up the comparison. `nr` indicates the number
#' of times the given comparison is made within the network.
#'
#' If there is only a single observation for each study within the dataset (i.e. as for standard
#' network meta-analysis) `nr` will represent the number of studies that compare treatments `t1` and
#' `t2`.
#'
#' If there are multiple observations for each study within the dataset (as in time-course MBNMA)
#' `nr` will represent the number of time points in the dataset in which treatments `t1` and `t2` are
#' compared.
#'
#' @examples
#' data <- data.frame(studyID=c(1,1,2,2,3,3,4,4,5,5,5),
#'   treatment=c(1,2,1,3,2,3,3,4,1,2,4)
#'   )
#'
#' # Identify comparisons informed by direct and indirect evidence
#' mb.comparisons(data)
#' @export
mb.comparisons <- function(data)
{
  # Assert checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data, add=argcheck)
  checkmate::assertNames(names(data), must.include = c("studyID", "treatment"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  data <- dplyr::arrange(data, studyID, treatment)

  t1 <- vector()
  t2 <- vector()

  for (i in seq_along(data[["studyID"]])) {

    k <- i+1

    while (k<=nrow(data) &
           data[["studyID"]][k] == data[["studyID"]][i] &
           !is.null(data[["studyID"]][k])) {

      # Ensures ordering of t1 to t2 is lowest to highest
      t <- sort(c(data[["treatment"]][i], data[["treatment"]][k]))

      t1 <- append(t1, t[1])
      t2 <- append(t2, t[2])

      k <- k+1
    }

  }

  comparisons <- data.frame("t1"=t1, "t2"=t2)

  comparisons <- comparisons %>%
    dplyr::group_by(t1, t2) %>%
    dplyr::mutate(nr=dplyr::n())

  comparisons <- unique(comparisons)
  comparisons <- dplyr::arrange(comparisons, t1, t2)

  return(comparisons)

}


#' Identify comparisons in time-course MBNMA datasets that fulfil criteria for node-splitting
#'
#' Identify comparisons informed by both direct and indirect evidence from independent sources in MBNMA
#' datasets with repeated measurements in each study. These comparisons are therefore those which
#' fulfil the criteria for testing for inconsistency via node-splitting, following the method of van
#' Valkenhoef \insertCite{vanvalkenhoef2016;textual}{MBNMAtime}.
#'
#' @inheritParams mb.run
#'
#' @inherit inconsistency.loops references details return
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' # Create mb.network object
#' network <- mb.network(osteopain)
#'
#' # Identify comparisons informed by direct and indirect evidence
#' mb.nodesplit.comparisons(network)
#' @export
mb.nodesplit.comparisons <- function(network)
{

  # Assert checks
  checkmate::assertClass(network, classes="mb.network", null.ok=FALSE)

  # First identify loops of evidence
  data <- get.latest.time(network)
  #data$treatment <- factor(data[["treatment"]], labels=network[["treatments"]])
  #data$treatment <- factor(as.numeric(data[["treatment"]]), labels=network[["treatments"]])
  #levels(data$treatment) <- treatments

  splits <- inconsistency.loops(data)

  return(splits)

}




#' Perform node-splitting on a MBNMA time-course network
#'
#' Within a MBNMA time-course network, split contributions into direct and indirect evidence and test
#' for consistency between them. Closed loops of treatments in which it is possible to test for
#' consistency are those in which direct and indirect evidence are available from independent sources
#' van Valkenhoef \insertCite{vanvalkenhoef2016;textual}{MBNMAtime}.
#'
#'
#' @param comparisons A data frame specifying the comparisons to be split (one row per comparison).
#' The frame has two columns indicating each treatment for each comparison: `t1` and `t2`.
#' @param nodesplit.parameters A character vector of named time-course parameters on which to
#' node-split (e.g. c("beta.1", "beta.2")). Can use "all" to split on all time-course parameters.
#' @param ... Arguments to be sent to `mb.run()`
#' @inheritParams mb.run
#'
#' @inherit mb.run details
#' @inherit inconsistency.loops references
#'
#' @return A an object of `class("mb.nodesplit")` that is a list containing elements
#' `d.X.Y` (treatment 1 = `X`, treatment 2 = `Y`). Each element (corresponding to each
#' comparison) contains additional numbered elements corresponding to each parameter in the
#' time-course function on which node splitting was performed. These elements then contain:
#' * `overlap matrix` MCMC results for the difference between direct and indirect evidence
#' * `p.values` Bayesian p-value for the test of consistency between direct and indirect evidence
#' * `quantiles`
#' * `forest.plot`
#' * `density.plot`
#' * `direct` MCMC results for the direct evidence
#' * `indirect` MCMC results for the indirect evidence
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Create mb.network object
#' painnet <- mb.network(osteopain)
#'
#' # Identify comparisons informed by direct and indirect evidence
#' splits <- mb.nodesplit.comparisons(painnet)
#'
#' # Fit a log-linear time-course MBNMA (takes a while to run)
#' result <- mb.nodesplit(painnet, comparisons=splits, nodesplit.parameters="all",
#'   fun=tloglin(pool.rate="rel", method.rate="common"),
#'   rho="dunif(0,1)", covar="varadj"
#'   )
#'
#' # Fit an emax time-course MBNMA with a node-split on emax parameters only
#' result <- mb.nodesplit(painnet, comparisons=splits, nodesplit.parameters="emax",
#'   fun=temax(pool.emax="rel", method.emax="common",
#'     pool.et50="rel", method.et50="common"))
#'
#' # Inspect results
#' print(result)
#' summary(result)
#'
#' # Plot results
#' plot(result)
#' }
#' @export
mb.nodesplit <- function(network, comparisons=mb.nodesplit.comparisons(network),
                            nodesplit.parameters="all", fun=tpoly(degree = 1),
                            ...
)
{
  # Required packages: overlapping

  # network = an object of class mb.network
  # comparisons = data frame specifying the comparisons to be split. The frame must contain columns: 't1' and 't2'
  # nodesplit.parameters = a vector of parameters on which to nodesplit (e.g. c("beta.1", "beta.2")). Can use "all" to split on all parameters.

  # Assert checks - other assertions will be raised in nested functions
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mb.network", null.ok=FALSE, add=argcheck)
  checkmate::assertDataFrame(comparisons, any.missing=FALSE, null.ok=FALSE, ncols=3, add=argcheck)
  checkmate::assertCharacter(nodesplit.parameters, null.ok=FALSE, unique=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check betas are specified correctly and prepare format for subsequent functions
  # for (i in 1:4) {
  #   betaname <- paste0("beta.", i)
  #   if (!is.null(get(betaname))) {
  #     assign(paste0(betaname, ".str"), compound.beta(get(betaname)))
  #   } else if (is.null(get(betaname))) {
  #     assign(paste0(betaname, ".str"), NULL)
  #   }
  # }

  # Get treatment labels
  trt.labs <- network$treatments

  # Check nodesplit.parameters has true values
  if (is.null(nodesplit.parameters)) {
    stop("No parameters have been specified on which to node-split")
  }

  if (!all(nodesplit.parameters %in% c("all", fun$params[which(fun$apool=="rel")]))) {
    stop("Parameter specified for nodesplit.parameters must be a parameter modelled using relative effects specified\nwithin the model")
  }


  if (nodesplit.parameters=="all") {
    UME <- TRUE
    nodesplit.parameters <- fun$params[which(fun$apool=="rel")]
  } else {
    UME <- nodesplit.parameters
  }

  if (length(nodesplit.parameters)==0) {
    stop("Parameter specified for nodesplit.parameters must be a parameter modelled using relative effects specified\nwithin the model")
  }


  ########### CHECKS OF DATASET FOR VALIDITY OF NODE-SPLITTING (possibly use Val Valkenhoef automation) ############
  print("running checks")

  # Ensure lowest treatment code in t1
  comparisons[["t1"]] <- as.numeric(comparisons[["t1"]])
  comparisons[["t2"]] <- as.numeric(comparisons[["t2"]])
  for (i in seq_along(comparisons[["t1"]])) {
    comparisons[i,1:2] <- sort(as.matrix(comparisons[,1:2])[i,])
  }


  ############# Run NMA model #############

  result.nma <- do.call(mb.run, args=list(network=network, fun=fun,
                                         ...))


  ######### Loop over all node splits in network ########

  nodesplit.results <- list()
  network.temp <- network

  for (row in seq_along(comparisons[["t1"]])) {

    comp <- as.numeric(c(comparisons[["t1"]][row], comparisons[["t2"]][row]))

    print(paste0("Comparison ", row,"/",nrow(comparisons)))
    print(paste0("Calculating nodesplit for: ", network$treatments[comp[1]], " vs ", network$treatments[comp[2]]))
    print(paste0("Treatment code: ", comp[1], " vs ", comp[2]))



    #######################################
    ######### For NMA model  ###########
    #######################################

    print("Running NMA model")

    nma.dif <- list()
    for (param in seq_along(nodesplit.parameters)) {
      if (grepl("beta", nodesplit.parameters[param])) {
        index <- which(nodesplit.parameters[param] %in% fun$params)
        node1 <- paste0("d.", index, "[", comp[1], "]")
        node2 <- paste0("d.", index, "[", comp[2], "]")
      } else{
        node1 <- paste0(nodesplit.parameters[param], "[", comp[1], "]")
        node2 <- paste0(nodesplit.parameters[param], "[", comp[2], "]")
      }

      nma1 <- result.nma$BUGSoutput$sims.matrix[,colnames(result.nma$BUGSoutput$sims.matrix)==node1]
      nma2 <- result.nma$BUGSoutput$sims.matrix[,colnames(result.nma$BUGSoutput$sims.matrix)==node2]

      nma.dif[[nodesplit.parameters[param]]] <- nma2 - nma1
    }


    #######################################
    ######### For direct model  ###########
    #######################################

    # Change network reference treatment to estimate direct effects
    if (!exists("results.dir") |
        network.temp$treatment[1] != network$treatment[comp[1]]) {

      if (network.temp$treatment[1] != network$treatment[comp[1]]) {
        data.temp <- network.temp$data.ab
        data.temp$treatment <- factor(data.temp$treatment, labels=network.temp$treatments)
        network.temp <- mb.network(data.temp, reference=network$treatment[comp[1]])
      }

      # Run UME model to estimate direct effects
      result.dir <- do.call(mb.run, args=list(network=network.temp, fun=fun, UME=UME,
                                              ...))
    }

    # Store required model parameter values
    comp.temp <- which(network$treatments[comp[2]] == network.temp$treatments)
    dir.dif <- list()
    for (param in seq_along(nodesplit.parameters)) {
      if (grepl("beta", nodesplit.parameters[param])) {
        index <- which(nodesplit.parameters[param] %in% fun$params)
        node <- paste0("d.", index, "[1,", comp.temp, "]")
      } else{
        node <- paste0(nodesplit.parameters[param], "[1,", comp.temp, "]")
      }

      dir.dif[[nodesplit.parameters[param]]] <-
        result.dir$BUGSoutput$sims.matrix[,colnames(result.dir$BUGSoutput$sims.matrix)==node]
    }

    print("Direct complete")



    #########################################
    ####### For MBNMA indirect model ########
    #########################################

    df <- network[["data.ab"]]

    # Remove comparisons to split on
    df <- drop.comp(df=df, comp=comp)

    if (!(1 %in% df$treatment)) {
      string <- paste("Reference treatment removed for node-split. Treatments have been reordered with the next lowest coded treatment as the reference:\ntreatment ",
                     min(df$treatment, na.rm=TRUE))
      warning(string)
    }

    ind.net <- mb.network(df)

    result.ind <- do.call(mb.run, args=list(network=ind.net, fun=fun, ...))

    ind.dif <- list()
    for (param in seq_along(nodesplit.parameters)) {
      if (grepl("beta", nodesplit.parameters[param])) {
        index <- which(nodesplit.parameters[param] %in% fun$params)
        node1 <- paste0("d.", index, "[", comp[1], "]")
        node2 <- paste0("d.", index, "[", comp[2], "]")
      } else{
        node1 <- paste0(nodesplit.parameters[param], "[", comp[1], "]")
        node2 <- paste0(nodesplit.parameters[param], "[", comp[2], "]")
      }

      ind1 <- result.ind$BUGSoutput$sims.matrix[,colnames(result.ind$BUGSoutput$sims.matrix)==node1]
      ind2 <- result.ind$BUGSoutput$sims.matrix[,colnames(result.ind$BUGSoutput$sims.matrix)==node2]

      ind.dif[[nodesplit.parameters[param]]] <- ind2 - ind1
    }

    print("Indirect complete")

    #########################################
    ######### Calculate Differences ##########
    #########################################

    nodesplit.comparison <- list()

    for (i in seq_along(ind.dif)) {

      overlap.mat <- list("direct"=dir.dif[[i]], "indirect"=ind.dif[[i]])
      overlap <- overlapping::overlap(overlap.mat, plot=FALSE)
      #p.values <- overlap$OV
      diff <- sum((dir.dif[[i]]-ind.dif[[i]])>0) / length(dir.dif[[i]])
      p.values <- min(diff, 1-diff)

      plot.df <- data.frame(source=c(rep("NMA", length(nma.dif[[i]])),
                                     rep("Direct", length(dir.dif[[i]])),
                                     rep("Indirect", length(ind.dif[[i]]))
                                     ),
                            value=c(nma.dif[[i]],
                                    dir.dif[[i]],
                                    ind.dif[[i]]
                                    )) %>%
        dplyr::mutate(Source=factor(source))


      # Quantiles
      quantile_dif <- stats::quantile(ind.dif[[i]] - dir.dif[[i]], c(0.025, 0.5, 0.975))
      quantile_dir <- stats::quantile(dir.dif[[i]], c(0.025, 0.5, 0.975))
      quantile_ind <- stats::quantile(ind.dif[[i]], c(0.025, 0.5, 0.975))
      quantile_nma <- stats::quantile(nma.dif[[i]], c(0.025, 0.5, 0.975))
      quantiles <- list("difference" = quantile_dif, "direct"=quantile_dir, "indirect"=quantile_ind, "nma"=quantile_nma)


      # Add plots for overlap and forest in return
      nodesplit <- list("comparison"= c(trt.labs[comp[2]], trt.labs[comp[1]]),
                        "parameter"=paste(names(ind.dif)[i],
                                          paste("[",comp[1], ",", comp[2], "]", sep=""),
                                          sep=""),
                        "overlap matrix"=overlap.mat,
                        "p.values"=p.values, "quantiles"=quantiles,
                        #"forest.plot"=gg, "density.plot"=dens,
                        "mcmc"=plot.df)

      nodesplit.comparison[[names(ind.dif)[i]]] <- nodesplit

    }

    nodesplit.results[[paste(trt.labs[comp[2]], trt.labs[comp[1]], sep=" vs ")]] <- nodesplit.comparison
  }

  class(nodesplit.results) <- "nodesplit"

  return(nodesplit.results)
}





#' Drops arms with comp treatments to generate dataset for indirect MBNMA
#' @noRd
drop.comp <- function(df, comp) {
  #x <- 1
  studies <- unique(df$studyID)

  for (i in seq_along(studies)) {
    # Separate single study
    subset <- subset(df, df$studyID==studies[i])

    if (all(comp %in% subset$treatment)) {
      # Remove study from dataset
      df <- subset(df, df$studyID!=studies[i])

      if (subset$narm[1]>2) {
        subset <- subset(subset, subset$treatment!=sample(comp,1))
        #subset <- subset(subset, subset$treatment!=comp[abs(x)+1])
        #x <- x-1

        # Reinsert ammended study
        df <- rbind(df, subset)
      }
    }
  }
  return(df)
}





# FUNCTION IS DEPRACATED!!!
check.design <- function(data, path) {
  # Ensure that no TWO comparisons in the loop share the same set of supporting studies
  # (not sure any of these steps are actually necessary)

  comb <- utils::combn(path, 2)
  designs <- list()
  for (comp in 1:ncol(comb)) {

    # Identify unique designs that include each comparison
    temp <- data[["design"]][apply(data,1, function(x) all(comb[,comp] %in% x[["design"]]))]

    if (length(temp)>0) {
      designs[[length(designs) + 1]] <-
        sort(unlist(unique(lapply(temp, function(x) paste(x, collapse="_")))))
    }
  }

  # Check that the sets of designs for each comparison in the loop are not identical
  # This is SLOW and could probably be changed to improve run time
  fail.design <- FALSE
  for (z in seq_along(designs)) {
    count <- sum(sapply(designs, FUN = identical, designs[[z]]))
    if (count>=2) {
      fail.design <- TRUE
    }
  }
  return(fail.design)
}


