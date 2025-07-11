#' Ancestor Regression
#'
#' This function performs ancestor regression for linear structural equation models \insertCite{schultheiss2024ancestorregressionstructuralvector}{AncReg} and
#' vector autoregressive models \insertCite{ancestor}{AncReg} with explicit error control for false discovery, at least asymptomatically.
#' @references
#'   \insertAllCited{}
#' @param x A named numeric matrix containing the observational data.
#' @param degree An integer specifying the order of the SVAR process to be considered. Default is 0 for no time series.
#' @param targets A character vector specifying the variables whose ancestors should be estimated. Default is all variables.
#' @param f A function specifying the non-linearity used in the ancestor regression. Default is a cubic function.
#' @return An object of class "AncReg" containing:
#'   \item{z.val}{A numeric matrix of test statistics.}
#'   \item{p.val}{A numeric matrix of p-values.}
#' @seealso \code{\link{summary.AncReg}}, \code{\link{instant_graph}}, \code{\link{summary_graph}},
#'   \code{\link{instant_p.val}}, \code{\link{summary_p.val}}
#' @export
#' @examples
#' ##### simulated example for inear structural equation models
#'
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
#' n <- 5000
#' N <- matrix(rexp(n * p), ncol = p)
#' X <- t(Bprime %*% t(N))
#' colnames(X) <- LETTERS[1:p]
#'
#' # fit ancestor regression
#' fit <- AncReg(X)
#' # collect ancestral p-values and graph
#' res <- summary(fit)
#' res
#'
#' #compare true and estimated ancestral graph
#' trueGraph <- igraph::graph_from_adjacency_matrix(recAncestor(B != 0))
#' ancGraph <- igraph::graph_from_adjacency_matrix(res$graph)
#'
#' oldpar <- par(mfrow = c(1, 2))
#' plot(trueGraph, main = 'true ancestral graph', vertex.size = 30)
#' plot(ancGraph, main = 'Ancestor Regression', vertex.size = 30)
#'
#' ##### SVAR-example with geyser timeseries
#' geyser <- MASS::geyser
#' # shift waiting such that it is waiting after erruption
#' geyser2 <- data.frame(waiting = geyser$waiting[-1], duration = geyser$duration[-nrow(geyser)])
#'
#' # fit ancestor regression with 6 lags considered
#' fit2 <- AncReg(as.matrix(geyser2), degree = 6)
#' res2 <- summary(fit2)
#' res2
#'
#' # visualize instantaneous ancestry
#' instGraph <- igraph::graph_from_adjacency_matrix(res2$inst.graph)
#' plot(instGraph, edge.label = round(diag(res2$inst.p.val[1:2, 2:1]), 2),
#'      main = 'instantaneous effects', vertex.size = 90)
#'
#' # visualize summary of lagged ancestry
#' sumGraph <- igraph::graph_from_adjacency_matrix(res2$sum.graph)
#' plot(sumGraph, edge.label = round(diag(res2$sum.p.val[1:2, 2:1]), 2),
#'      main = 'summary graph', vertex.size = 90)
#' par(oldpar)
AncReg <- function(x, degree = 0, targets = colnames(x), f = function(x) x^3){
  # input control
  if(!is.matrix(x) || !is.numeric(x)){
    stop("x must be a numeric matrix")
  }
  if(is.null(colnames(x))){
    stop("x must have column names")
  }
  if(!is.numeric(degree) || length(degree) != 1 || degree < 0){
    stop("degree must be a non-negative integer")
  }
  if(!is.character(targets) || length(targets) == 0){
    stop("targets must be a non-empty character vector")
  }
  if(!all(targets %in% colnames(x))){
    stop("all targets must be columns of x")
  }
  if(min(dim(x)) < 2){
    stop("x must have at least 2 rows and 2 columns")
  }
  if(!is.function(f)){
    stop("f must be a function")
  }
  if(!is.numeric(f(x[, targets[1]])) || length(f(x[, targets[1]])) != nrow(x)){
    stop("f must be a function that returns a numeric vector of the same length as x")
  }



  cols <- colnames(x) # all variables
  ind <- which(cols %in% targets) # considered columns
  p <- ncol(x) # number of variables
  xt <- matrix(NA, nrow = nrow(x), ncol = p * (degree +1)) # matrix to store data with lagged predictors
  for(j in 1:p){
    xtj <- tsutils::lagmatrix(x[,j], 0:degree) # all lags for given variable
    xt[,j + p * (0:degree)] <- xtj # store to combined matrix
  }
  xt <- xt[complete.cases(xt),] # ignore incomplete rows
  colnames(xt) <- paste(rep(cols, degree + 1), ".",
                        rep(0:degree, each = p), sep ="") # column names with degree information

  if(degree > 0){
    # residuals from VAR process
    u <- xt[,1:p] - xt[,-c(1:p)] %*%
      solve(crossprod(xt[,-c(1:p)])) %*% crossprod(xt[,-c(1:p)], xt[,1:p])
  } else {
    # original data if no time series considered
    u <- xt
  }
  n <- nrow(u) # number of available full observations
  # matrix to store test-statistics
  z.val <- matrix(NA, nrow = length(targets), ncol = p * (degree + 1),
                  dimnames = list(targets, colnames(xt)))

  for (s in 0:degree){
    # loop over lags
    if(degree > 0){
      # project out data s + 1 to s + degree steps before
      us <- xt[(s+1):n,1:p] - xt[1:(n-s),-c(1:p)] %*%
        solve(crossprod(xt[1:(n-s),-c(1:p)]))  %*%
        crossprod(xt[1:(n-s),-c(1:p)], xt[(s+1):n,1:p])
    } else {
      # original data if no time series considered
      us <- xt
    }

    colnames(us) <- colnames(x)
    wiz <- 1:p + s*p # position to store the z-value

    for (i in 1:length(targets)){
      # loop over target
      tari <- f(us[,targets[i]]) # apply non-linearity to target
      su <- summary(lm(tari ~u[1: (n - s), ])) # fit linear model
      z.val[i, wiz] <- su$coefficients[-1,3] # store z-values
    }
  }

  p.val <- pnorm(abs(z.val), lower.tail = FALSE)*2 # p-values
  res <- list(z.val = z.val, p.val = p.val)
  class(res) <- "AncReg"
  return(res)
}
