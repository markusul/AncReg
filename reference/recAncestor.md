# Recursive Ancestor Detection

This function recursively detects all ancestors of a given set of
variables in a matrix. Adds ancestors of ancestors to the output matrix.

## Usage

``` r
recAncestor(pmat)
```

## Arguments

- pmat:

  A boolean matrix indicating whether a connection was detected.

## Value

A boolean matrix indicating whether a connection was detected or
constructed.

## See also

[`AncReg`](http://www.markus-ulmer.ch/AncReg/reference/AncReg.md),
[`summary_graph`](http://www.markus-ulmer.ch/AncReg/reference/summary_graph.md)

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

# edge effects to adjecency matrix
B <- B != 0

# transform adjacency matrix to ancestral matrix
recAncestor(B)
#>       A     B     C     D     E
#> A FALSE FALSE FALSE FALSE FALSE
#> B FALSE FALSE FALSE FALSE FALSE
#> C FALSE  TRUE FALSE FALSE FALSE
#> D FALSE  TRUE FALSE FALSE FALSE
#> E  TRUE  TRUE  TRUE  TRUE FALSE
```
