
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AncReg

<!-- badges: start -->

[![R-CMD-check](https://github.com/markusul/AncReg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/markusul/AncReg/actions/workflows/R-CMD-check.yaml)
![CRAN Downloads
overall](https://cranlogs.r-pkg.org/badges/grand-total/AncReg?color=brightgreen)
<!-- badges: end -->

Ancestor Regression (AncReg) is a package with methods to test for
ancestral connections in linear structural equation models (Christoph
Schultheiss, Ulmer, and Bühlmann (2024)) and structural vector
autoregressive models (C. Schultheiss and Bühlmann (2023)). Ancestor
Regression provides explicit error control for false causal discovery,
at least asymptotically. To have power, however, it relies on
non-Gaussian distributions for power.

## Installation

To install the Ancestor Regression R package from CRAN, just run

``` r
install.packages(AncReg)
```

You can install the development version of AncReg from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("markusul/AncReg")
```

or

``` r
# install.packages('pak')
pak::pkg_install('markusul/AncReg')
```

## linear structural equation models

This is a basic example on how to use Ancestor Regression using some
simulated data.

``` r
library(AncReg)

# random DAGS for simulation
set.seed(42)

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

# fit Ancestor Regression
fit <- AncReg(X)
#> Registered S3 method overwritten by 'quantmod':
#>   method            from
#>   as.zoo.data.frame zoo
#> Registered S3 methods overwritten by 'tsutils':
#>   method          from   
#>   print.nemenyi   greybox
#>   summary.nemenyi greybox
fit
#> $z.val
#>            A.0        B.0        C.0         D.0        E.0
#> A 13.680849514  0.5970868  0.3703936  0.18152226 -0.6038712
#> B -1.120155016 20.4051461  0.5304711  1.85402719 -1.5040844
#> C  0.003835857 -1.2963950 22.8533011  0.71919705 -0.9608690
#> D  2.369113587  0.8629459  0.3454023 22.88105801 -0.5901320
#> E  3.602545269 -0.1764219 -1.7062338  0.06736191 18.6628463
#> 
#> $p.val
#>            A.0          B.0           C.0           D.0          E.0
#> A 1.321398e-42 5.504495e-01  7.110892e-01  8.559577e-01 5.459293e-01
#> B 2.626477e-01 1.505113e-92  5.957853e-01  6.373527e-02 1.325596e-01
#> C 9.969394e-01 1.948394e-01 1.355009e-115  4.720195e-01 3.366180e-01
#> D 1.783078e-02 3.881672e-01  7.297919e-01 7.173985e-116 5.551022e-01
#> E 3.151165e-04 8.599625e-01  8.796457e-02  9.462936e-01 9.929839e-78
#> 
#> attr(,"class")
#> [1] "AncReg"
```

The summary function can be used to collect and organize the p-values.
Additionally it returns estimated ancestral graphs.

``` r
# collect ancestral p-values and graph
res <- summary(fit)
res
#> $p.val
#>              A         B          C          D         E
#> A 1.0000000000 0.5504495 0.71108921 0.85595766 0.5459293
#> B 0.2626477098 1.0000000 0.59578535 0.06373527 0.1325596
#> C 0.9969394362 0.1948394 1.00000000 0.47201951 0.3366180
#> D 0.0178307764 0.3881672 0.72979192 1.00000000 0.5551022
#> E 0.0003151165 0.8599625 0.08796457 0.94629359 1.0000000
#> 
#> $graph
#>       A     B     C     D     E
#> A FALSE FALSE FALSE FALSE FALSE
#> B FALSE FALSE FALSE FALSE FALSE
#> C FALSE FALSE FALSE FALSE FALSE
#> D FALSE FALSE FALSE FALSE FALSE
#> E  TRUE FALSE FALSE FALSE FALSE
#> 
#> $alpha
#> [1] 0.05
#> 
#> attr(,"class")
#> [1] "summary.AncReg"
```

As we know the truth in the simulated model, we can compare the
estimated ancestral graph with the true one.

``` r
#compare true and estimated ancestral graph
trueGraph <- igraph::graph_from_adjacency_matrix(recAncestor(B != 0))
ancGraph <- igraph::graph_from_adjacency_matrix(res$graph)

par(mfrow = c(1, 2))
plot(trueGraph, main = 'true ancestral graph', vertex.size = 30)
plot(ancGraph, main = 'Ancestor Regression', vertex.size = 30)
```

<img src="man/figures/README-lsem-plot-1.png" width="100%" />

## structural vector autoregressive models

We show an example of a SVAR application using the time series of geyser
eruptions. (C. Schultheiss and Bühlmann (2023))

``` r
geyser <- MASS::geyser
# shift waiting such that it is waiting after erruption
geyser2 <- data.frame(waiting = geyser$waiting[-1], duration = geyser$duration[-nrow(geyser)])

# fit ancestor regression with 6 lags considered
fit2 <- AncReg(as.matrix(geyser2), degree = 6)
res2 <- summary(fit2)
res2
#> $inst.p.val
#>            waiting     duration
#> waiting  1.0000000 0.0004811719
#> duration 0.5109396 1.0000000000
#> 
#> $inst.graph
#>          waiting duration
#> waiting    FALSE     TRUE
#> duration   FALSE    FALSE
#> 
#> $inst.alpha
#> [1] 0.05
#> 
#> $sum.p.val
#>            waiting    duration
#> waiting  1.0000000 0.008733271
#> duration 0.1760936 1.000000000
#> 
#> $sum.graph
#>          waiting duration
#> waiting    FALSE     TRUE
#> duration   FALSE    FALSE
#> 
#> attr(,"class")
#> [1] "summary.AncReg"

# visualize instantaneous ancestry
instGraph <- igraph::graph_from_adjacency_matrix(res2$inst.graph)
plot(instGraph, edge.label = round(diag(res2$inst.p.val[1:2, 2:1]), 2), 
     main = 'instantaneous effects', vertex.size = 90)
```

<img src="man/figures/README-svar-1.png" width="100%" />

``` r

# visualize summary of lagged ancestry
sumGraph <- igraph::graph_from_adjacency_matrix(res2$sum.graph)
plot(sumGraph, edge.label = round(diag(res2$sum.p.val[1:2, 2:1]), 2), 
     main = 'summary graph', vertex.size = 90)
```

<img src="man/figures/README-svar-2.png" width="100%" />

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-ancestor" class="csl-entry">

Schultheiss, C, and P Bühlmann. 2023. “Ancestor Regression in Linear
Structural Equation Models.” *Biometrika* 110 (4): 1117–24.
<https://doi.org/10.1093/biomet/asad008>.

</div>

<div id="ref-schultheiss2024ancestorregressionstructuralvector"
class="csl-entry">

Schultheiss, Christoph, Markus Ulmer, and Peter Bühlmann. 2024.
“Ancestor Regression in Structural Vector Autoregressive Models.”
<https://arxiv.org/abs/2403.03778>.

</div>

</div>
