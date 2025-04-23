#' Recursive Ancestor Detection
#'
#' This function recursively detects all ancestors of a given set of variables in a matrix.
#' Adds ancestors of ancestors to the output matrix.
#' @param pmat A boolean matrix indicating whether a connection was detected.
#' @return A boolean matrix indicating whether a connection was detected or constructed.
#' @export
#' @seealso \code{\link{AncReg}}, \code{\link{summary_graph}}
#' @examples
#' # random DAGS for simulation
#' set.seed(1234)
#'
#' p <- 5 #number of nodes
#' DAG <- pcalg::randomDAG(p, prob = 0.5)
#'
#' B <- matrix(0, p, p) # represent DAG as matrix
#' for (i in 2:p){
#'   for(j in 1:(i-1)){
#'     # store edge weights
#'     B[i,j] <- max(0, DAG@edgeData@data[[paste(j,"|",i, sep="")]]$weight)
#'   }
#' }
#' colnames(B) <- rownames(B) <- LETTERS[1:p]
#'
#' # edge effects to adjecency matrix
#' B <- B != 0
#'
#' # transform adjacency matrix to ancestral matrix
#' recAncestor(B)
recAncestor <- function(pmat) {
  # get variable names
  preds <- colnames(pmat)
  targets <- rownames(pmat)
  ancmat <- pmat # matrix shape for output
  for (i in 1:length(targets)){
    # loop through targets
    # which predictors have been checked, if for some no ancestors where search, do not use these
    tested <- setdiff(preds, targets)
    an <- names(which(pmat[i,])) # names of detected ancestors
    while (length(setdiff(an, tested)) > 0) {
      # if not all ancestors have been checked
      for (k in setdiff(an, tested)) {
        # loop through unchecked
        j <- which(targets == k) # find position of ancestor
        an <- unique(c(an, names(which(pmat[j,])))) # add ancestors of ancestors
        tested <- c(tested, k) # declare it as checked
      }
    }
    ancmat[i, which(preds %in% an)] <- TRUE # add additional ancestors to output matrix
  }
  return(ancmat)
}
