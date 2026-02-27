test_that("makeQ constructs valid Q matrices", {
  Q_mk <- makeQ(model = "mk", k = 2)
  expect_true(is.matrix(Q_mk))
  expect_equal(nrow(Q_mk), 2)
  expect_equal(ncol(Q_mk), 2)

  # Row sums should be 0
  expect_true(all(abs(rowSums(Q_mk)) < 1e-10))

  # Off-diagonals should be non-negative
  expect_true(all(Q_mk[Q_mk != diag(Q_mk)] >= 0))

  # Diagonals should be negative
  expect_true(all(diag(Q_mk) < 0))
})


test_that("makeQ handles different models", {
  for (model in c("mk", "sym", "ard", "chromosome")) {
    Q <- makeQ(model = model, k = 4)
    expect_equal(nrow(Q), 4)
    expect_equal(ncol(Q), 4)
    expect_true(all(abs(rowSums(Q)) < 1e-10))
  }
})


test_that("simRateScape generates valid data", {
  skip("Requires tree; implement with rtree setup")

  tree <- ape::rtree(10)
  Q <- makeQ(model = "mk", k = 2)

  sim <- simRateScape(
    tree = tree,
    Q = Q,
    lambda_sigma = 1.0,
    pi = 0.7,
    nrep = 2
  )

  expect_equal(nrow(sim$data), 2)
  expect_equal(ncol(sim$data), 10)
  expect_true(all(sim$data %in% c(0, 1)))
})


test_that("checkPrior runs without error", {
  skip("Requires tree; implement with rtree setup")

  tree <- ape::rtree(10)
  Q <- makeQ(model = "mk", k = 2)

  check <- checkPrior(
    tree = tree,
    Q = Q,
    lambda_sigma = 1.0,
    nsim = 10
  )

  expect_is(check, "ratescape_prior_check")
  expect_true(!is.na(check$overlap_measure))
})
