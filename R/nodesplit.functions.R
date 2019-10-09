# Functions for nodesplit models
# Author: Hugo Pedder
# Date created: 2018-09-10


#' Identify comparisons in loops that fulfill criteria for node-splitting
#'
#' Identify comparisons informed by both direct and indirect evidence from
#' independent sources, which therefore fulfill the criteria for testing for
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

  splits <- data.frame("t1"=splits1, "t2"=splits2, "path"=paths, "loops"=loops)

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


#' Identify comparisons in time-course MBNMA datasets that fulfill criteria for node-splitting
#'
#' Identify comparisons informed by both direct and indirect evidence from independent sources in MBNMA
#' datasets with repeated measurements in each study. Thes comparisons are therefore those which
#' fulfill the criteria for testing for inconsistency via node-splitting, following the method of van
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
#' nodesplit (e.g. c("beta.1", "beta.2")). Can use "all" to split on all time-course parameters.
#' @param ... Arguments to be sent to `mb.run()`
#' @inheritParams mb.emax.hill
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
#' # Create mb.network object
#' network <- mb.network(osteopain)
#'
#' # Identify comparisons informed by direct and indirect evidence
#' splits <- mb.nodesplit.comparisons(network)
#'
#' # Fit an exponential time-course MBNMA
#' result <- mb.nodesplit(network, comparisons=splits, nodesplit.parameters="all",
#'   fun="exponential",
#'   beta.1=list(pool="rel", method="common"))
#'
#' # Fit an emax time-course MBNMA with a nodesplit on emax parameters only
#' result <- mb.nodesplit(network, comparisons=splits, nodesplit.parameters="beta.1",
#'   fun="emax",
#'   beta.1=list(pool="rel", method="random"),
#'   beta.2=list(pool="rel", method="common")
#'   )
#'
#' # Inspect results
#' print(result)
#' summary(result)
#'
#' # Plot results
#' plot(result)
#' @export
mb.nodesplit <- function(network, comparisons=mb.nodesplit.comparisons(network),
                            nodesplit.parameters="all", fun="linear", user.fun=NULL,
                            beta.1=list(pool="rel", method="common"), beta.2=NULL, beta.3=NULL, beta.4=NULL,
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
  for (i in 1:4) {
    betaname <- paste0("beta.", i)
    if (!is.null(get(betaname))) {
      assign(paste0(betaname, ".str"), compound.beta(get(betaname)))
    } else if (is.null(get(betaname))) {
      assign(paste0(betaname, ".str"), NULL)
    }
  }

  # Get treatment labels
  trt.labs <- network$treatments


  # Checks - NEEDS TO BE CHANGED TO REFLECT MODEL CODE...COULD HAPPEN AFTER mb.write and then
  #used to look inside the written model code to check numbers of parameters etc.
  # "d.2" is not a parameter in current model with single prior - it is "beta.2"
  #if (is.null(parameters.to.save)) {
  #  if (fun=="linear" | fun=="exponential") {
  #    parameters.to.save=c("d.1")
  #  } else if (fun=="emax") {
  #    parameters.to.save=c("d.1", "d.2")
  #  }
  #}

  # Check nodesplit.parameters has true values
  if (is.null(nodesplit.parameters)) {
    stop("No parameters have been specified on which to nodesplit")
  }

  if (!all(nodesplit.parameters %in% c("all", "beta.1", "beta.2", "beta.3", "beta.4"))) {
    stop("Parameter specified for nodesplit.parameters does not exist in model. They must be of the type 'beta.1', 'beta.2', etc.")
  }

  # Check that nodesplit.parameters and class effects are not on same parameter
  #for (i in seq_along(names(class.effect))) {
  #  if (names(class.effect)[i] %in% nodesplit.parameters) {
  #    stop(paste0("Node splitting cannot be applied to parameters with class effects. ", names(class.effect)[i], " has been given a class effect"))
  #  }
  #}

  if (nodesplit.parameters=="all") {
    UME <- TRUE
    nodesplit.parameters <- vector()
    treat.params <- c("beta.1", "beta.2", "beta.3", "beta.4")

    for (i in seq_along(treat.params)) {
      if (!is.null(get(treat.params[i]))) {
        if (get(treat.params[i])$pool == "rel") {
          nodesplit.parameters <- append(nodesplit.parameters, treat.params[i])
        }
      }
    }

  } else {
    UME <- nodesplit.parameters
  }


  # Separate ... arguments for mb.write and mb.run


  ########### CHECKS OF DATASET FOR VALIDITY OF NODE-SPLITTING (possibly use Val Valkenhoef automation) ############
  print("running checks")

  # Ensure lowest treatment code in t1
  comparisons[["t1"]] <- as.numeric(comparisons[["t1"]])
  comparisons[["t2"]] <- as.numeric(comparisons[["t2"]])
  for (i in seq_along(comparisons[["t1"]])) {
    comparisons[i,1:2] <- sort(comparisons[i,1:2])
  }


  ######## Write JAGS scripts ########

  #These can be removed later as now incorporated into mb.run
  # Only needs to capture key time-course parameters in model for gen.paramaters.to.save
  model.ind <- mb.write(fun=fun, user.fun=user.fun, beta.1=beta.1.str, beta.2=beta.2.str,
                           beta.3=beta.3.str, beta.4=beta.4.str, UME=FALSE
                           )

  #model.dir <- mb.write(fun=fun, user.fun=user.fun, alpha=alpha, beta.1=beta.1, beta.2=beta.2,
  #                         beta.3=beta.3, beta.4=beta.4,
  #                         positive.scale=positive.scale, intercept=intercept, rho=rho, covar=covar,
  #                         class.effect=class.effect, UME=UME)

  # parameters.to.save <- gen.parameters.to.save(parameters.to.save=parameters.to.save,
  #                                              model.params = c(1,2,3,4),
  #                                              model = model.ind
  # )

  parameters.to.save <-
    gen.parameters.to.save(model.params=c(1,2,3,4), model=model.ind)


  ############# Run NMA model #############

  result.nma <- mb.run(network, parameters.to.save=parameters.to.save, fun=fun,
                       beta.1=beta.1, beta.2=beta.2, beta.3=beta.3, beta.4=beta.4,
                       ...
                       )


  ######### Loop over all node splits in network ########

  nodesplit.results <- list()
  network.temp <- network

  for (row in seq_along(comparisons[["t1"]])) {

    comp <- as.numeric(c(comparisons[["t1"]][row], comparisons[["t2"]][row]))


    #######################################
    ######### For NMA model  ###########
    #######################################

    print("Running NMA model")

    nma.dif <- list()
    nodesplit.parameters <- sort(nodesplit.parameters)

    for (param in 1:4) {
      if (paste0("beta.", param) %in% nodesplit.parameters) {

        node1 <- paste0("d.", param, "[", comp[1], "]")
        node2 <- paste0("d.", param, "[", comp[2], "]")

        nma1 <- result.nma$BUGSoutput$sims.matrix[,colnames(result.nma$BUGSoutput$sims.matrix)==node1]
        nma2 <- result.nma$BUGSoutput$sims.matrix[,colnames(result.nma$BUGSoutput$sims.matrix)==node2]

        nma.dif[[paste0("beta.", param)]] <- nma2 - nma1
      }
    }


    #######################################
    ######### For direct model  ###########
    #######################################

    # dir.dif <- list()
    # nodesplit.parameters <- sort(nodesplit.parameters)
    #
    # for (param in 1:4) {
    #   if (paste0("beta.", param) %in% nodesplit.parameters) {
    #     node <- paste0("d.", param, "[", comp[1], ",", comp[2], "]")
    #     #dir.dif[[length(dir.dif)+1]] <-
    #     dir.dif[[paste0("beta.", param)]] <-
    #       result.dir$BUGSoutput$sims.matrix[,colnames(result.dir$BUGSoutput$sims.matrix)==node]
    #   }
    # }

    # Change network reference treatment to estimate direct effects
    if (!exists("results.dir") |
        network.temp$treatment[1] != network$treatment[comp[1]]) {

      if (network.temp$treatment[1] != network$treatment[comp[1]]) {
        data.temp <- network.temp$data.ab
        data.temp$treatment <- factor(data.temp$treatment, labels=network.temp$treatments)
        network.temp <- mb.network(data.temp, reference=network$treatment[comp[1]])
      }

      # Run UME model to estimate direct effects
      result.dir <- mb.run(network.temp, fun=fun, parameters.to.save=parameters.to.save,
                              beta.1=beta.1, beta.2=beta.2, beta.3=beta.3, beta.4=beta.4,
                              UME=UME,
                              ...
      )
    }

    # Store required model parameter values
    comp.temp <- which(network$treatments[comp[2]] == network.temp$treatments)
    dir.dif <- list()
    nodesplit.parameters <- sort(nodesplit.parameters)

    for (param in 1:4) {
      if (paste0("beta.", param) %in% nodesplit.parameters) {
        node <- paste0("d.", param, "[1,", comp.temp, "]")
        #dir.dif[[length(dir.dif)+1]] <-
        dir.dif[[paste0("beta.", param)]] <-
          result.dir$BUGSoutput$sims.matrix[,colnames(result.dir$BUGSoutput$sims.matrix)==node]
      }
    }

    print("Direct complete")
    print(paste0("Treatment name ", network$treatments[comp[1]], " vs ", network$treatments[comp[2]]))
    print(paste0("Treatment code ", comp[1], " vs ", comp[2]))



    #########################################
    ####### For MBNMA indirect model ########
    #########################################

    data <- network[["data.ab"]]

    # Remove comparisons to split on
    data <- drop.comp(data=data, comp=comp)

    if (!(1 %in% data$treatment)) {
      string <- paste("Reference treatment removed for node-split. Treatments have been reordered with the next lowest coded treatment as the reference:\ntreatment ",
                     min(data$treatment, na.rm=TRUE))
      warning(string)
    }

    data <- mb.network(data)

    result.ind <- mb.run(data, parameters.to.save=parameters.to.save, fun=fun,
                            beta.1=beta.1, beta.2=beta.2, beta.3=beta.3, beta.4=beta.4,
                            ...
    )

    ind.dif <- list()
    nodesplit.parameters <- sort(nodesplit.parameters)

    for (param in 1:4) {
      if (paste0("beta.", param) %in% nodesplit.parameters) {

        node1 <- paste0("d.", param, "[", comp[1], "]")
        node2 <- paste0("d.", param, "[", comp[2], "]")

        ind1 <- result.ind$BUGSoutput$sims.matrix[,colnames(result.ind$BUGSoutput$sims.matrix)==node1]
        ind2 <- result.ind$BUGSoutput$sims.matrix[,colnames(result.ind$BUGSoutput$sims.matrix)==node2]

        #ind.dif[[length(ind.dif)+1]] <- ind2 - ind1
        ind.dif[[paste0("beta.", param)]] <- ind2 - ind1
      }
    }


    print("Indirect complete")
    print(paste0("Treatment name ", network$treatments[comp[1]], " vs ", network$treatments[comp[2]]))
    print(paste0("Treatment code ", comp[1], " vs ", comp[2]))


    #########################################
    ######### Calculate Differences ##########
    #########################################

    nodesplit.comparison <- list()

    for (i in seq_along(ind.dif)) {

      overlap.mat <- list("direct"=dir.dif[[i]], "indirect"=ind.dif[[i]])
      overlap <- overlapping::overlap(overlap.mat, plot=FALSE)
      p.values <- overlap$OV


      # Quantiles
      quantile_dif <- stats::quantile(ind.dif[[i]] - dir.dif[[i]], c(0.025, 0.5, 0.975))
      quantile_dir <- stats::quantile(dir.dif[[i]], c(0.025, 0.5, 0.975))
      quantile_ind <- stats::quantile(ind.dif[[i]], c(0.025, 0.5, 0.975))
      quantile_nma <- stats::quantile(nma.dif[[i]], c(0.025, 0.5, 0.975))
      quantiles <- list("difference" = quantile_dif, "direct"=quantile_dir, "indirect"=quantile_ind, "nma"=quantile_nma)


      # GGplots

      source <- c("NMA", "Direct", "Indirect")
      l95 <- c(quantile_nma[1], quantile_dir[1], quantile_ind[1])
      med <- c(quantile_nma[2], quantile_dir[2], quantile_ind[2])
      u95 <- c(quantile_nma[3], quantile_dir[3], quantile_ind[3])
      plotdata <- data.frame(source, l95, med, u95)

      #title <- paste0(names(ind.dif)[i], ": Treatment ", comp[2], " vs Treatment ", comp[1])
      title <- paste0(names(ind.dif)[i], ": ", trt.labs[comp[2]], " vs ", trt.labs[comp[1]])

      gg <-
        ggplot2::ggplot(data=plotdata, ggplot2::aes(x=source, y=med, ymin=l95, ymax=u95)) +
        ggplot2::geom_pointrange() +
        #ggplot2::geom_hline(yintercept=0, lty=2) +
        ggplot2::coord_flip() +  # flip coordinates (puts labels on y axis)
        ggplot2::xlab("") + ggplot2::ylab("Treatment effect (95% CrI)") + ggplot2::ggtitle(title) +
        ggplot2::theme(axis.text = ggplot2::element_text(size=15),
                       axis.title = ggplot2::element_text(size=18),
                       title=ggplot2::element_text(size=18)) +
        ggplot2::theme(plot.margin=ggplot2::unit(c(1,1,1,1),"cm"))


      # Density plots (with shaded area of overlap)
      #MBNMA <- Sowers2005_nma[1:length(Direct)]
      molten <- data.frame(ind.dif[[i]], dir.dif[[i]])
      molten <- reshape2::melt(molten, measure.vars=names(molten))
      names(molten) <- c("Estimate", "value")
      linetypes <- c("solid", "dash")
      levels(molten$Estimate) <- c("Indirect", "Direct")

      dens <- ggplot2::ggplot(molten, ggplot2::aes(x=molten$value, linetype=molten$Estimate, fill=molten$Estimate)) +
        ggplot2::geom_density(alpha=0.2) +
        ggplot2::xlab(title) +
        ggplot2::ylab("Posterior density") +
        #scale_linetype_manual(values=factor(molten$Estimate)) +
        #scale_fill_manual(name="Evidence") +
        ggplot2::theme(strip.text.x = ggplot2::element_text(size=12)) +
        ggplot2::theme(axis.text = ggplot2::element_text(size=12),
                       axis.title = ggplot2::element_text(size=14))


      # Add plots for overlap and forest in return
      nodesplit <- list("comparison"= c(trt.labs[comp[2]], trt.labs[comp[1]]),
                        "parameter"=paste("d.", i,
                                          paste("[",comp[1], ",", comp[2], "]", sep=""),
                                          sep=""),
                        "overlap matrix"=overlap.mat,
                        "p.values"=p.values, "quantiles"=quantiles,
                        "forest.plot"=gg, "density.plot"=dens,
                        "direct"=result.dir, "indirect"=result.ind)

      #nodesplit.comparison[[length(nodesplit.comparison)+1]] <- nodesplit
      nodesplit.comparison[[names(ind.dif)[i]]] <- nodesplit

    }

    nodesplit.results[[paste("d",comp[1],comp[2], sep=".")]] <- nodesplit.comparison

  }

  class(nodesplit.results) <- "mb.nodesplit"

  return(nodesplit.results)
}





#' Drops arms with comp treatments to generate dataset for indirect MBNMA
#' @noRd
drop.comp <- function(data, comp) {
  x <- 1
  for (i in seq_along(unique(data$studyID))) {
    # Separate single study
    subset <- subset(data, data$studyID==unique(data$studyID)[i])

    if (all(comp %in% subset$treatment)) {
      # Remove study from dataset
      data <- subset(data, data$studyID!=unique(data$studyID)[i])

      if (subset$narm[1]>2) {
        subset <- subset(subset, subset$treatment!=comp[abs(x)+1])
        x <- x-1

        # Reinsert ammended study
        data <- rbind(data, subset)
      }
    }
  }
  return(data)
}




# FUNCTION IS DEPRACATED!!!
# check.path <- function(data, dropdata,
#                        comparisons=c(comparisons[i,1], comparisons[i,2]),
#                        path, graph) {
#   # Identify if there is still an indirect pathway not involving direct evidence studies
#
#   # Comparisons is a numeric vector of length 2 indicating the treatment comparison on which to be split
#
#   # Then check for new paths...
#
#   path.fail <- TRUE
#
#   while(as.logical(is.finite(igraph::shortest.paths(
#     graph, comparisons[1], comparisons[2]))) == TRUE) {
#
#     # Delete current path from graph to check for next shortest path
#     for (edge in 1:length(path)-1) {
#       del.index <- which(apply(igraph::as_edgelist(graph), 1,
#                                function(x) identical(x, as.character(path[edge:(edge+1)]))))
#       graph <- igraph::as.undirected(igraph::delete_edges(graph, del.index))
#     }
#
#     # Check if the indirect path is not included in any study that has been dropped (for node-splitting)
#     # And if so break from the function with a pass
#     if (any(apply(dropdata,1, function(x) all(path %in% x[["design"]]))) == FALSE) {
#       path.fail <- FALSE
#       break()
#       #} else {
#       #  path <- as.numeric(igraph::shortest_paths(igraph::as.undirected(graph),
#       #                                            comparisons[1], comparisons[2],
#       #                                            weights=NA
#       #  )[["vpath"]][[1]])
#       #}
#
#       # ADD THIS SECTION IF GETTING WARNINGS FROM IGRAPH
#     } else if (as.logical(is.finite(igraph::shortest.paths(
#       graph, comparisons[1], comparisons[2], weights=NULL))) == TRUE){
#
#       path <- as.numeric(igraph::shortest_paths(igraph::as.undirected(graph),
#                                                 comparisons[1], comparisons[2],
#                                                 weights=NA
#       )[["vpath"]][[1]])
#     } else {break()}
#
#   }
#
#   return(list("path.fail"=path.fail, "path"=path))
# }




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


