# Instantaneous graph

Construct instantaneous graph from p-values and significance level.
Recursively constructs all ancestral connections by adding ancestors of
ancestors.

## Usage

``` r
instant_graph(lin.anc, alpha = 0.05, verbose = FALSE, corr = TRUE)
```

## Arguments

- lin.anc:

  output from AncReg()

- alpha:

  significance level

- verbose:

  should information be printed?

- corr:

  should multiplicity correction be applied?

## Value

A list containing:

- rec.ancs:

  A boolean matrix indicating whether one variable affects another
  instantaneously

- alpha:

  The significance level to avoid cycles

## See also

[`AncReg`](http://www.markus-ulmer.ch/AncReg/reference/AncReg.md),
[`instant_p.val`](http://www.markus-ulmer.ch/AncReg/reference/instant_p.val.md)

## Examples

``` r
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
X <- t(Bprime %*% t(N))
colnames(X) <- LETTERS[1:p]

# fit ancestor regression
fit <- AncReg(X)

# generate instantaneous graph
instant_graph(fit, alpha = 0.01, verbose = TRUE)
#> $rec.ancs
#>       A     B     C     D     E
#> A FALSE FALSE FALSE FALSE FALSE
#> B FALSE FALSE FALSE FALSE FALSE
#> C FALSE FALSE FALSE FALSE FALSE
#> D FALSE FALSE FALSE FALSE FALSE
#> E FALSE FALSE FALSE FALSE FALSE
#> 
#> $alpha
#> [1] 0.01
#> 
```
