# Statistical methods used in the reanalysis. Kept free of IST specifics so they can be tested on
# simulated data.

library(survival)

# Newcombe hybrid score interval for a difference of two proportions (method 10, Newcombe 1998)
wilson <- function(x, n, z = qnorm(0.975)) {
  p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(centre - half, centre + half)
}

risk_difference <- function(y, arm) {
  ok <- !is.na(y)
  y <- y[ok]; arm <- arm[ok]
  x1 <- sum(y[arm == 1]); n1 <- sum(arm == 1)
  x0 <- sum(y[arm == 0]); n0 <- sum(arm == 0)
  p1 <- x1 / n1; p0 <- x0 / n0
  w1 <- wilson(x1, n1); w0 <- wilson(x0, n0)
  d <- p1 - p0
  lo <- d - sqrt((p1 - w1[1])^2 + (w0[2] - p0)^2)
  hi <- d + sqrt((w1[2] - p1)^2 + (p0 - w0[1])^2)
  p <- suppressWarnings(chisq.test(matrix(c(x1, n1 - x1, x0, n0 - x0), 2), correct = FALSE)$p.value)
  list(events_treated = x1, n_treated = n1, events_control = x0, n_control = n0,
       risk_treated = p1, risk_control = p0, rd = d, low = lo, high = hi, p = p)
}

# Covariate adjusted risk difference by standardisation (g computation): fit a logistic model,
# predict every patient's risk under treatment and under control, and average.
standardised_rd <- function(formula, data, arm) {
  fit <- glm(formula, data = data, family = binomial)
  d1 <- data; d1[[arm]] <- 1
  d0 <- data; d0[[arm]] <- 0
  mean(predict(fit, d1, type = "response")) - mean(predict(fit, d0, type = "response"))
}

adjusted_rd <- function(formula, data, arm, reps = as.integer(Sys.getenv("REPS", 400)), seed = 1) {
  data <- data[complete.cases(data[, all.vars(formula)]), ]
  est <- standardised_rd(formula, data, arm)
  set.seed(seed)
  boots <- replicate(reps, standardised_rd(formula, data[sample.int(nrow(data), replace = TRUE), ], arm))
  fit <- glm(formula, data = data, family = binomial)
  co <- summary(fit)$coefficients[arm, ]
  list(rd = est, low = unname(quantile(boots, 0.025)), high = unname(quantile(boots, 0.975)),
       se_boot = sd(boots), odds_ratio = exp(co[["Estimate"]]),
       or_low = exp(co[["Estimate"]] - 1.96 * co[["Std. Error"]]),
       or_high = exp(co[["Estimate"]] + 1.96 * co[["Std. Error"]]),
       p = co[["Pr(>|z|)"]], n = nrow(data))
}

# Treatment effect within each level of a baseline factor, plus a likelihood ratio test for interaction
subgroup <- function(data, y, arm, by, covars) {
  data <- data[!is.na(data[[y]]) & !is.na(data[[by]]), ]
  rows <- lapply(levels(droplevels(data[[by]])), function(l) {
    s <- data[data[[by]] == l, ]
    r <- risk_difference(s[[y]], s[[arm]])
    data.frame(subgroup = by, level = l, n = nrow(s), risk_treated = r$risk_treated,
               risk_control = r$risk_control, rd = r$rd, low = r$low, high = r$high)
  })
  base <- paste(c(arm, covars), collapse = " + ")
  f0 <- as.formula(paste(y, "~", base, "+", by))
  f1 <- as.formula(paste(y, "~", base, "+", by, "+", arm, ":", by))
  m0 <- glm(f0, data = data, family = binomial)
  m1 <- glm(f1, data = data, family = binomial)
  p <- anova(m0, m1, test = "LRT")[2, "Pr(>Chi)"]
  out <- do.call(rbind, rows)
  out$p_interaction <- p
  out
}

holm <- function(p) p.adjust(p, method = "holm")

# Kaplan Meier curve at chosen days, and an adjusted Cox model with a proportional hazards check
km_table <- function(time, event, arm, days) {
  fit <- survfit(Surv(time, event) ~ arm)
  s <- summary(fit, times = days, extend = TRUE)
  data.frame(arm = sub("arm=", "", as.character(s$strata)), day = s$time, surv = s$surv,
             low = s$lower, high = s$upper, at_risk = s$n.risk)
}

cox <- function(formula, data, arm) {
  fit <- coxph(formula, data = data, model = TRUE, x = TRUE)
  co <- summary(fit)$conf.int[arm, ]
  z <- cox.zph(fit)
  list(hr = unname(co[1]), low = unname(co[3]), high = unname(co[4]),
       p = unname(summary(fit)$coefficients[arm, "Pr(>|z|)"]),
       ph_test_p_treatment = unname(z$table[arm, "p"]), ph_test_p_global = unname(z$table["GLOBAL", "p"]))
}

# Discrimination and calibration of a logistic prognostic model under K fold cross validation
auc <- function(y, p) {
  r <- rank(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cv_prognostic <- function(formula, data, k = 5, seed = 1) {
  data <- data[complete.cases(data[, all.vars(formula)]), ]
  set.seed(seed)
  fold <- sample(rep(1:k, length.out = nrow(data)))
  p <- numeric(nrow(data))
  for (i in 1:k) {
    fit <- glm(formula, data = data[fold != i, ], family = binomial)
    p[fold == i] <- predict(fit, data[fold == i, ], type = "response")
  }
  y <- data[[all.vars(formula)[1]]]
  lp <- qlogis(p)
  slope <- unname(coef(glm(y ~ lp, family = binomial))[2])
  intercept <- unname(coef(glm(y ~ offset(lp), family = binomial))[1])
  dec <- cut(p, quantile(p, 0:10 / 10), include.lowest = TRUE, labels = FALSE)
  list(auc = auc(y, p), calibration_slope = slope, calibration_in_the_large = intercept, n = nrow(data),
       deciles = data.frame(decile = 1:10, predicted = tapply(p, dec, mean), observed = tapply(y, dec, mean)))
}

# Standardised mean difference for baseline balance (continuous or 0/1)
smd <- function(x, arm) {
  a <- x[arm == 1]; b <- x[arm == 0]
  (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2)
}
