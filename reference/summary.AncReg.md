# Summary of AncReg

Summarize the results of AncReg. For models with degree = 0 only the
instantaneous graph is returned and for models with degree \> 0 the
summary graph is returned as well.

## Usage

``` r
# S3 method for class 'AncReg'
summary(object, alpha = 0.05, verbose = FALSE, corr = TRUE, ...)
```

## Arguments

- object:

  output from AncReg()

- alpha:

  significance level for determin whether a connection is significant

- verbose:

  should information be printed?

- corr:

  should multiplicity correction be applied?

- ...:

  Further arguments passed to or from other methods.

## Value

A list containing: If `degree = 0`:

- p.val:

  A numeric matrix of p-values for the instantaneous graph

- graph:

  A boolean matrix indicating whether one variable affects another
  instantaneously

- alpha:

  The significance level to avoid cycles

If `degree > 0`:

- inst.p.val:

  A numeric matrix of p-values for the instantaneous graph

- inst.graph:

  A boolean matrix indicating whether one variable affects another
  instantaneously

- inst.alpha:

  The significance level to avoid cycles

- sum.p.val:

  A numeric matrix of p-values for the summary graph

- sum.graph:

  A boolean matrix indicating whether one variable affects another

## See also

[`AncReg`](http://www.markus-ulmer.ch/AncReg/reference/AncReg.md),
[`instant_graph`](http://www.markus-ulmer.ch/AncReg/reference/instant_graph.md),
[`summary_graph`](http://www.markus-ulmer.ch/AncReg/reference/summary_graph.md),
[`instant_p.val`](http://www.markus-ulmer.ch/AncReg/reference/instant_p.val.md),
[`summary_p.val`](http://www.markus-ulmer.ch/AncReg/reference/summary_p.val.md)

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
# collect ancestral p-values and graph
res <- summary(fit, alpha = 1)
res
#> $p.val
#>           A           B         C         D         E
#> A 1.0000000 0.355083093 0.0493693 0.5107251 0.9777551
#> B 0.6612617 1.000000000 0.5350391 0.1672114 0.5645382
#> C 0.3190314 0.224148437 1.0000000 0.9467940 0.9968308
#> D 0.8740048 0.006213275 0.3530013 1.0000000 0.2938822
#> E 0.4328651 0.798593250 0.4182837 0.1692897 1.0000000
#> 
#> $graph
#>       A     B     C     D     E
#> A FALSE FALSE  TRUE FALSE FALSE
#> B FALSE FALSE FALSE FALSE FALSE
#> C FALSE FALSE FALSE FALSE FALSE
#> D FALSE  TRUE FALSE FALSE FALSE
#> E FALSE FALSE FALSE FALSE FALSE
#> 
#> $alpha
#> [1] 1
#> 
#> attr(,"class")
#> [1] "summary.AncReg"
```
