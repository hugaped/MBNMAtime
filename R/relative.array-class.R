##############################################
#### Functions for class("relative.array") ####
##############################################


#' Print posterior medians (95% credible intervals) for table of relative effects/mean
#' differences between treatments/classes
#'
#' @inheritParams predict.mbnma
#' @param ... further arguments passed to `knitr::kable`
#'
#' @export
print.relative.array <- function(x, digits=1, ...) {

  xmat <- x$relarray

  outmat <- matrix(nrow=nrow(xmat), ncol=ncol(xmat))
  # dimnames(outmat)[[1]] <- dimnames(xmat)[[1]]
  # dimnames(outmat)[[2]] <- dimnames(xmat)[[2]]

  for (i in 1:nrow(xmat)) {
    for (k in 1:ncol(xmat)) {
      if (!is.na(xmat[i,k,1])) {
        outmat[i,k] <- neatCrI(quantile(xmat[i,k,], probs=c(0.025, 0.5, 0.975)), digits = digits)
      }
    }
  }
  diag(outmat) <- dimnames(xmat)[[1]]

  cat(crayon::bold(paste0("========================================\nTreatment comparisons at time = ", x$time, "\n========================================\n")))
  cat("\n")
  #knitr::kable(outmat, ...)

  write.table(format(outmat, justify="centre"), row.names = FALSE, col.names = FALSE, quote=FALSE)
}


