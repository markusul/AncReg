#' @importFrom Rdpack reprompt
#' @importFrom stats complete.cases lm pnorm
#' @importFrom methods is

holm.corr <- function(pv, cut = TRUE){
  # function to perform Bonferroni-Holm correction on matrix with possibility to ignore some entries
  # Input
  # pv (numeric, matrix): p-values, unused p-values should be 1
  # cut (boolean): should correction cut p-values at 1
  # Output
  # numeric matrix: corrected p-values

  np <- length(pv) # total number of p-values
  ntest <- np - sum(rownames(pv) %in% colnames(pv)) # total number of tests, i.e., without self-effects
  i <- c(seq_len(ntest), rep(ntest, np - ntest)) # used for correction factor
  o <- order(pv) # ranking of p-values
  ro <- order(o) # ranking of rankings to sort back
  if (cut) {
    # correct by factor according to ordering, cut at 1, sort back
    pv[] <- pmin(1, cummax((ntest + 1L - i) * pv[o]))[ro]
  }
  else {
    # correct by factor according to ordering, sort back
    pv[] <- cummax((ntest + 1L - i) * pv[o])[ro]
  }
  pv
}
