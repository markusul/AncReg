#' Summary of AncReg
#'
#' Summarize the results of AncReg. For models with degree = 0 only the instantaneous graph is returned
#' and for models with degree > 0 the summary graph is returned as well.
#' @param object output from AncReg()
#' @param alpha significance level for determin whether a connection is significant
#' @param verbose should information be printed?
#' @param corr should multiplicity correction be applied?
#' @param ... Further arguments passed to or from other methods.
#' @return A list containing:
#' If \code{degree = 0}:
#' \item{p.val}{A numeric matrix of p-values for the instantaneous graph}
#' \item{graph}{A boolean matrix indicating whether one variable affects another instantaneously}
#' \item{alpha}{The significance level to avoid cycles}
#' If \code{degree > 0}:
#' \item{inst.p.val}{A numeric matrix of p-values for the instantaneous graph}
#' \item{inst.graph}{A boolean matrix indicating whether one variable affects another instantaneously}
#' \item{inst.alpha}{The significance level to avoid cycles}
#' \item{sum.p.val}{A numeric matrix of p-values for the summary graph}
#' \item{sum.graph}{A boolean matrix indicating whether one variable affects another}
#' @method summary AncReg
#' @seealso \code{\link{AncReg}}, \code{\link{instant_graph}}, \code{\link{summary_graph}},
#' \code{\link{instant_p.val}}, \code{\link{summary_p.val}}
#' @export
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
#' # solution in terms of noise
#' Bprime <- MASS::ginv(diag(p) - B)
#'
#' n <- 500
#' N <- matrix(rexp(n * p), ncol = p)
#' X <- t(Bprime %*% t(N))
#' colnames(X) <- LETTERS[1:p]
#'
#' # fit ancestor regression
#' fit <- AncReg(X)
#' # collect ancestral p-values and graph
#' res <- summary(fit, alpha = 1)
#' res
summary.AncReg <- function(object, alpha = 0.05, verbose = FALSE, corr = TRUE, ...){
  inst.p.val <- instant_p.val(object)
  inst.graph <- instant_graph(object, alpha = alpha, verbose = verbose, corr = corr)

  last_name <- colnames(object$p.val)
  last_name <- last_name[length(last_name)]
  el <- strsplit(last_name, '')[[1]]
  if(el[length(el)] == "0"){
    res <- list(p.val = inst.p.val, graph = inst.graph$rec.ancs, alpha = inst.graph$alpha)
    class(res) <- "summary.AncReg"
    return(res)
  }

  sum.p.val <- summary_p.val(object)
  sum.graph <- summary_graph(object, alpha = alpha, corr = corr)

  res <- list(inst.p.val = inst.p.val, inst.graph = inst.graph$rec.ancs,
              inst.alpha = inst.graph$alpha,
              sum.p.val = sum.p.val, sum.graph = sum.graph)
  class(res) <- "summary.AncReg"
  return(res)
}

#' P-values for summary graph
#'
#' Collect p-values for summary graph.
#' @param lin.anc output from AncReg()
#' @return A numeric matrix of p-values for the summary graph
#' @seealso \code{\link{AncReg}}
#' @export
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
#' # solution in terms of noise
#' Bprime <- MASS::ginv(diag(p) - B)
#'
#' n <- 500
#' N <- matrix(rexp(n * p), ncol = p)
#' X <- t(Bprime %*% t(N))
#' colnames(X) <- LETTERS[1:p]
#'
#' # fit ancestor regression
#' fit <- AncReg(X)
#'
#' # collect summary p-values
#' summary_p.val(fit)
summary_p.val <- function(lin.anc){
  # check input format
  if(!is(lin.anc, "AncReg")) {
    stop("lin.anc must be output from AncReg()")
  }
  pv <- lin.anc$p.val

  # get variables names
  preds <- colnames(pv)
  targets <- rownames(pv)
  p <- sum(grepl("\\.0", preds)) # number of predictors per time step
  npv <- ncol(pv)/p # number of time steps
  spv <- pv[,1:p] # matrix shape for summary p-value
  colnames(spv) <- gsub("\\.0.*","",colnames(spv)) # no time indicator
  spv[] <- 1
  for(i in 1:length(targets)){
    # loop through targets
    for (j in 1:p){
      # loop through predictors
      if(targets[i] == colnames(spv)[j]) next() # omit self-effects
      # get p-values at all time steps for that predictor
      pij <- pv[i, seq(j, j + p*(npv - 1), length.out = npv)]
      # combine p-values
      spv[i,j] <- min(min(sort(pij) * npv /(1:npv)) * sum(1/(1:npv)), 1)
    }
  }
  return(spv)
}

#' Summary graph
#'
#' Construct summary graph from p-values and significance level.
#' Recursively constructs all ancestral connections by adding ancestors of ancestors.
#' @param lin.anc output from AncReg()
#' @param alpha significance level
#' @param corr should multiplicity correction be applied?
#' @return A boolean matrix indicating whether one variable affects another
#' @seealso \code{\link{AncReg}}, \code{\link{summary_p.val}}
#' @export
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
#' # solution in terms of noise
#' Bprime <- MASS::ginv(diag(p) - B)
#'
#' n <- 500
#' N <- matrix(rexp(n * p), ncol = p)
#' X <- t(Bprime %*% t(N))
#' colnames(X) <- LETTERS[1:p]
#'
#' # fit ancestor regression
#' fit <- AncReg(X)
#'
#' # generate summary graph
#' summary_graph(fit, alpha = 0.1)
summary_graph <- function(lin.anc, alpha = 0.05, corr = TRUE) {
  spv <- summary_p.val(lin.anc) # get summary p-values
  if(corr){
    # apply Bonferroni-Holm if needed
    pv.corr <- holm.corr(spv)
  } else {
    pv.corr <- spv
  }
  pmat <- pv.corr < alpha # convert to boolean
  # construct recursively
  return(recAncestor(pmat))
}

#' P-values for instantaneous graph
#'
#' Collect p-values for instantaneous graph.
#' @param lin.anc output from AncReg()
#' @return A numeric matrix of p-values for the instantaneous graph
#' @seealso \code{\link{AncReg}}
#' @export
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
#' # solution in terms of noise
#' Bprime <- MASS::ginv(diag(p) - B)
#'
#' n <- 500
#' N <- matrix(rexp(n * p), ncol = p)
#' X <- t(Bprime %*% t(N))
#' colnames(X) <- LETTERS[1:p]
#'
#' # fit ancestor regression
#' fit <- AncReg(X)
#'
#' # collect instantaneous p-values
#' instant_p.val(fit)
instant_p.val <- function(lin.anc){
  # check input format
  if(!is(lin.anc, "AncReg")) {
    stop("lin.anc must be output from AncReg()")
  }
  pv <- lin.anc$p.val

  # get variable names
  preds <- colnames(pv)
  targets <- rownames(pv)

  p <- sum(grepl("\\.0", preds)) # number of instant predictors
  if (p > 0){
    pv <- pv[,1:p] # instant p-values
    colnames(pv) <- gsub("\\.0.*","",colnames(pv)) # shorten column names
  }

  for(ta in targets) pv[ta, ta] <- 1 # no self-effects
  return(pv)
}

#' Instantaneous graph
#'
#' Construct instantaneous graph from p-values and significance level.
#' Recursively constructs all ancestral connections by adding ancestors of ancestors.
#' @param lin.anc output from AncReg()
#' @param alpha significance level
#' @param verbose should information be printed?
#' @param corr should multiplicity correction be applied?
#' @return A list containing:
#' \item{rec.ancs}{A boolean matrix indicating whether one variable affects another instantaneously}
#' \item{alpha}{The significance level to avoid cycles}
#' @seealso \code{\link{AncReg}}, \code{\link{instant_p.val}}
#' @export
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
#' # solution in terms of noise
#' Bprime <- MASS::ginv(diag(p) - B)
#'
#' n <- 500
#' N <- matrix(rexp(n * p), ncol = p)
#' X <- t(Bprime %*% t(N))
#' colnames(X) <- LETTERS[1:p]
#'
#' # fit ancestor regression
#' fit <- AncReg(X)
#'
#' # generate instantaneous graph
#' instant_graph(fit, alpha = 0.01, verbose = TRUE)
instant_graph <- function(lin.anc, alpha = 0.05, verbose = FALSE, corr = TRUE){
  pv <- instant_p.val(lin.anc)
  targets <- rownames(pv)

  # apply Bonferroni-Holm if needed
  if (corr){
    pv.corr <- holm.corr(pv)
  } else {
    pv.corr <- pv
  }

  anc1 <- pv.corr < alpha # detected instantaneous ancestors
  anc <- recAncestor(anc1) # reconstruct further
  if(sum(sapply(targets, function(ta) anc[ta , ta])) == 0){
    # if there are no cycles
    pv.corr[] <- anc # format for output matrix
    return(list(rec.ancs = pv.corr > 0, alpha = alpha))
  } else {
    loop.vars <- names(which(sapply(targets, function(ta) anc[ta , ta]))) # find variables that cause loops
    pvs.mat <- pv.corr[loop.vars, loop.vars] # get submatrix of p-values
    pvs.sub <- pvs.mat[pvs.mat < alpha] # get significant p-values
    new.alpha <- max(pvs.sub) # new significance level
    if(verbose) print(paste("Try decreasing alpha from", round(alpha, 3), "to", round(new.alpha, 3)))

    # create graph with submatrix and new significance level
    pvs.mat <- list(p.val = pvs.mat)
    class(pvs.mat) <- "AncReg"
    out <- instant_graph(pvs.mat, new.alpha, verbose = verbose, corr = FALSE)
    anc1[loop.vars, loop.vars] <- out$rec.ancs == 1 # adapt detected ancestors
    # construct recursively again
    return(list(rec.ancs = recAncestor(anc1), alpha = out$alpha))
  }
}
