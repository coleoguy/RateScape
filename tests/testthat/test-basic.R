test_that("makeQ produces valid rate matrices", {
  # Mk model
  Q <- makeQ("mk", nstates = 3, params = list(rate = 0.5))
  expect_equal(nrow(Q), 3)
  expect_equal(ncol(Q), 3)
  expect_equal(rowSums(Q), rep(0, 3), tolerance = 1e-10)
  expect_true(all(Q[row(Q) != col(Q)] >= 0))

  # Chromosome model
  Q_chr <- makeQ("chromosome", params = list(
    gain = 0.2, loss = 0.3, polyploidy = 0.01, max_chrom = 10
  ))
  expect_equal(nrow(Q_chr), 10)
  expect_equal(rowSums(Q_chr), rep(0, 10), tolerance = 1e-10)

  # Gain: state i -> i+1

  expect_equal(Q_chr[3, 4], 0.2)
  # Loss: state i -> i-1
  expect_equal(Q_chr[3, 2], 0.3)
  # Polyploidy: state 3 -> state 6
  expect_equal(Q_chr[3, 6], 0.01)
  # No polyploidy if 2*i > max_chrom
  expect_equal(Q_chr[6, 10], 0)  # 2*6=12 > 10, but let's check 6->10 is 0
})

test_that("discretize_gamma produces valid rate categories", {
  rates <- RateScape:::discretize_gamma(1.0, 4)
  expect_equal(length(rates), 4)
  expect_true(all(rates > 0))
  # Mean should be approximately 1
  expect_equal(mean(rates), 1.0, tolerance = 0.1)
})

test_that("stationary_dist sums to 1", {
  Q <- makeQ("mk", nstates = 4, params = list(rate = 0.5))
  pi <- RateScape:::stationary_dist(Q)
  expect_equal(sum(pi), 1.0, tolerance = 1e-10)
  expect_true(all(pi >= 0))
})

test_that("simRateScape produces valid output", {
  library(ape)
  set.seed(42)
  tree <- rtree(20)
  Q <- makeQ("mk", nstates = 3, params = list(rate = 0.5))
  data <- simRateScape(tree, Q)

  expect_equal(length(data), 20)
  expect_true(all(data >= 1 & data <= 3))
  expect_equal(names(data), tree$tip.label)
})

test_that("pruning_likelihood returns finite values", {
  library(ape)
  set.seed(42)
  tree <- rtree(20)
  Q <- makeQ("mk", nstates = 3, params = list(rate = 0.5))
  data <- simRateScape(tree, Q)
  rates <- rep(1.0, Nedge(tree))

  ll <- RateScape:::pruning_likelihood_scaled(tree, data, Q, rates, "equal")
  expect_true(is.finite(ll))
  expect_true(ll < 0)  # Log-likelihood should be negative
})

test_that("higher rates on true-fast branches improve likelihood", {
  library(ape)
  set.seed(42)
  tree <- rtree(30)
  Q <- makeQ("mk", nstates = 4, params = list(rate = 0.3))

  # Simulate with some fast branches
  true_rates <- rep(1, Nedge(tree))
  true_rates[1:5] <- 5.0
  data <- simRateScape(tree, Q, rates = true_rates)

  # Likelihood under homogeneous model
  ll_homo <- RateScape:::pruning_likelihood_scaled(
    tree, data, Q, rep(1, Nedge(tree)), "equal"
  )

  # Likelihood under true rates
  ll_true <- RateScape:::pruning_likelihood_scaled(
    tree, data, Q, true_rates, "equal"
  )

  # True rates should fit better (or at least as well)
  expect_true(ll_true >= ll_homo - 1)  # Allow small tolerance
})
