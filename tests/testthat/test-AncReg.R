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

# solution in terms of noise
Bprime <- MASS::ginv(diag(p) - B)

n <- 100
N <- matrix(rexp(n * p), ncol = p)
X <- t(Bprime %*% t(N))
colnames(X) <- LETTERS[1:p]

res <- AncReg(X)

test_that("AncReg works", {
  expect_equal(names(res), c("z.val", "p.val"))
  expect_false(any(unlist(lapply(res, is.na))))
})

# check functionality with custom inputs
targets <- c('A', 'B')
degree <- 10
res2 <- AncReg(X, degree = degree, targets = targets, exp)

test_that("AncReg works with time series", {
  expect_equal(dim(res2$p.val), c(length(targets), ncol(X) * (degree + 1)))
})

# in case of degree = 0 summary and instant is equal
test_that("degree = 0 -> summary = instant", {
  expect_equal(summary_p.val(res), instant_p.val(res))
  expect_equal(summary_graph(res), instant_graph(res)$rec.ancs)

  p.val <- res$p.val
  diag(p.val) <- 1
  expect_true(all(p.val == instant_p.val(res)))
})

test_that("holm correction can only increase p values", {
  expect_true(all(holm.corr(res$p.val) >= res$p.val))
})





