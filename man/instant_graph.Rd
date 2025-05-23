% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{instant_graph}
\alias{instant_graph}
\title{Instantaneous graph}
\usage{
instant_graph(lin.anc, alpha = 0.05, verbose = FALSE, corr = TRUE)
}
\arguments{
\item{lin.anc}{output from AncReg()}

\item{alpha}{significance level}

\item{verbose}{should information be printed?}

\item{corr}{should multiplicity correction be applied?}
}
\value{
A list containing:
\item{rec.ancs}{A boolean matrix indicating whether one variable affects another instantaneously}
\item{alpha}{The significance level to avoid cycles}
}
\description{
Construct instantaneous graph from p-values and significance level.
Recursively constructs all ancestral connections by adding ancestors of ancestors.
}
\examples{
# random DAGS for simulation
set.seed(1234)

p <- 5 #number of nodes
DAG <- pcalg::randomDAG(p, prob = 0.5)

B <- matrix(0, p, p) # represent DAG as matrix
for (i in 2:p){
  for(j in 1:(i-1)){
    # store edge weights
    B[i,j] <- max(0, DAG@edgeData@data[[paste(j,"|",i, sep="")]]$weight)
  }
}
colnames(B) <- rownames(B) <- LETTERS[1:p]

# solution in terms of noise
Bprime <- MASS::ginv(diag(p) - B)

n <- 500
N <- matrix(rexp(n * p), ncol = p)
X <- t(Bprime \%*\% t(N))
colnames(X) <- LETTERS[1:p]

# fit ancestor regression
fit <- AncReg(X)

# generate instantaneous graph
instant_graph(fit, alpha = 0.01, verbose = TRUE)
}
\seealso{
\code{\link{AncReg}}, \code{\link{instant_p.val}}
}
