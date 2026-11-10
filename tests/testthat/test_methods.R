source(file.path("..", "..", "R", "methods.R"))
source(file.path("..", "..", "R", "prepare.R"))

test_that("Newcombe interval matches the published worked example", {
  # Newcombe (1998) Table II example (e): 56/70 vs 48/80, difference 0.2000, interval 0.0524 to 0.3339
  y <- c(rep(1, 56), rep(0, 14), rep(1, 48), rep(0, 32))
  arm <- c(rep(1, 70), rep(0, 80))
  r <- risk_difference(y, arm)
  expect_equal(r$rd, 0.2, tolerance = 1e-9)
  expect_equal(r$low, 0.0524, tolerance = 5e-4)
  expect_equal(r$high, 0.3339, tolerance = 5e-4)
})

test_that("standardised difference equals the crude difference when there are no covariates", {
  set.seed(2)
  d <- data.frame(arm = rbinom(4000, 1, 0.5))
  d$y <- rbinom(4000, 1, 0.3 + 0.1 * d$arm)
  expect_equal(standardised_rd(y ~ arm, d, "arm"), mean(d$y[d$arm == 1]) - mean(d$y[d$arm == 0]), tolerance = 1e-8)
})

test_that("adjustment recovers a known treatment effect under a prognostic covariate", {
  set.seed(3)
  n <- 20000
  d <- data.frame(arm = rbinom(n, 1, 0.5), x = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-0.2 + 1.5 * d$x - 0.3 * d$arm))
  truth <- mean(plogis(-0.2 + 1.5 * d$x - 0.3)) - mean(plogis(-0.2 + 1.5 * d$x))
  r <- adjusted_rd(y ~ arm + x, d, "arm", reps = 50)
  expect_lt(abs(r$rd - truth), 0.01)
  expect_true(r$low < truth && truth < r$high)
})

test_that("subgroup interaction test is not significant when the effect is the same everywhere", {
  set.seed(4)
  n <- 6000
  d <- data.frame(arm = rbinom(n, 1, 0.5), g = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- rbinom(n, 1, plogis(-0.3 * d$arm))
  s <- subgroup(d, "y", "arm", "g", character(0))
  expect_equal(nrow(s), 3)
  expect_gt(unique(s$p_interaction), 0.01)
})

test_that("AUC is 0.5 for noise and 1 for perfect ranking", {
  expect_equal(auc(c(0, 0, 1, 1), c(0.1, 0.2, 0.8, 0.9)), 1)
  set.seed(5)
  expect_lt(abs(auc(rbinom(5000, 1, 0.4), runif(5000)) - 0.5), 0.03)
})

test_that("prepare merges the pilot H code into medium dose and codes missing 6 month outcomes as NA", {
  raw <- data.frame(RXASP = c("Y", "N", "N"), RXHEP = c("H", "L", "N"), AGE = 70, SEX = "M", RCONSC = "F",
                    STYPE = c("OTH", "LACS", "TACS"), RDELAY = 10, RSBP = 150, RATRIAL = c("Y", "", "N"),
                    RCT = "Y", RVISINF = "N", RASP3 = "N", RHEP24 = "N", OCCODE = c(2, 9, 0), ID14 = 0, ISC14 = 0,
                    H14 = 0, TRAN14 = 0, STRK14 = 0, TD = c(100, 200, 5), DIED = c(0, 0, 1), CMPLASP = "Y",
                    DNOSTRK = "N", CNTRYNUM = 1, stringsAsFactors = FALSE)
  d <- prepare(raw)
  expect_equal(as.character(d$hep_dose), c("M", "L", "N"))
  expect_equal(d$dd6, c(1L, NA, NA))
  expect_equal(d$time, c(100, 183, 5))
  expect_equal(d$died, c(0L, 0L, 1L))
  expect_true(is.na(d$af[2]))
})
