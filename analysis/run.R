# Reanalysis of the International Stroke Trial (IST, 19,435 patients, 1991 to 1996).
# Writes docs/data/results.json for the dashboard and report, and tables to outputs/.
#
#   bash scripts/get_data.sh
#   Rscript analysis/run.R

suppressPackageStartupMessages({ library(jsonlite); library(mice); library(splines) })
source("R/prepare.R")
source("R/methods.R")

raw <- read.csv("data/ist_raw.csv", stringsAsFactors = FALSE)
checks <- data_checks(raw)
d <- prepare(raw)
d$af_c <- addNA(d$af)  # keep patients whose AF status was not recorded, as their own level
dir.create("outputs", showWarnings = FALSE)
dir.create("docs/data", recursive = TRUE, showWarnings = FALSE)

COV <- c("ns(age, 3)", "sex", "consc", "stype", "ns(delay, 3)", "ns(sbp, 3)", "af_c", "infarct_visible", "ct_done")
f_adj <- function(y, arm, other) as.formula(paste(y, "~", paste(c(arm, other, COV), collapse = " + ")))

# ---- baseline balance -------------------------------------------------------------------------------
bal_vars <- list(
  "Age, years (mean)" = d$age, "Female" = as.integer(d$sex == "F"), "Drowsy or unconscious" = as.integer(d$consc != "Alert"),
  "Total anterior circulation stroke" = as.integer(d$stype == "TACS"), "Atrial fibrillation" = as.integer(d$af == "Y"),
  "Hours from stroke to randomisation (mean)" = d$delay, "Systolic BP, mmHg (mean)" = d$sbp,
  "CT scan before randomisation" = d$ct_done, "Aspirin in previous 3 days" = d$asp_before)
baseline <- do.call(rbind, lapply(names(bal_vars), function(v) {
  x <- bal_vars[[v]]
  cont <- grepl("mean", v)
  data.frame(variable = v,
             aspirin = mean(x[d$asp == 1], na.rm = TRUE), no_aspirin = mean(x[d$asp == 0], na.rm = TRUE),
             heparin = mean(x[d$hep == 1], na.rm = TRUE), no_heparin = mean(x[d$hep == 0], na.rm = TRUE),
             smd_aspirin = smd(x, d$asp), smd_heparin = smd(x, d$hep), continuous = cont)
}))
write.csv(baseline, "outputs/baseline.csv", row.names = FALSE)

# ---- primary outcome: dead or dependent at 6 months --------------------------------------------------
primary <- list(
  aspirin = list(crude = risk_difference(d$dd6, d$asp), adjusted = adjusted_rd(f_adj("dd6", "asp", "hep"), d, "asp")),
  heparin = list(crude = risk_difference(d$dd6, d$hep), adjusted = adjusted_rd(f_adj("dd6", "hep", "asp"), d, "hep"))
)
dead6 <- list(aspirin = risk_difference(d$dead6, d$asp), heparin = risk_difference(d$dead6, d$hep))
factorial_p <- anova(glm(f_adj("dd6", "asp", "hep"), data = d, family = binomial),
                     glm(update(f_adj("dd6", "asp", "hep"), . ~ . + asp:hep), data = d, family = binomial),
                     test = "LRT")[2, "Pr(>Chi)"]

# ---- 14 day events: benefit and harm per 1,000 patients --------------------------------------------
ev <- c(death14 = "Death", isc14 = "Recurrent ischaemic stroke", haem14 = "Haemorrhagic stroke",
        bleed14 = "Major extracranial bleed", stroke_death14 = "Death or recurrent stroke")
events14 <- do.call(rbind, lapply(names(ev), function(e) {
  rbind(
    with(risk_difference(d[[e]], d$asp), data.frame(comparison = "Aspirin vs none", event = ev[[e]], treated = risk_treated, control = risk_control, per_1000 = 1000 * rd, low = 1000 * low, high = 1000 * high, p = p)),
    with(risk_difference(d[[e]][d$hep_dose != "M"], d$hep[d$hep_dose != "M"]), data.frame(comparison = "Low dose heparin vs none", event = ev[[e]], treated = risk_treated, control = risk_control, per_1000 = 1000 * rd, low = 1000 * low, high = 1000 * high, p = p)),
    with(risk_difference(d[[e]][d$hep_dose != "L"], d$hep[d$hep_dose != "L"]), data.frame(comparison = "Medium dose heparin vs none", event = ev[[e]], treated = risk_treated, control = risk_control, per_1000 = 1000 * rd, low = 1000 * low, high = 1000 * high, p = p))
  )
}))
rownames(events14) <- NULL
write.csv(events14, "outputs/events_14_days.csv", row.names = FALSE)

# heparin by dose on the primary outcome, and a trend test across dose
dose <- data.frame(dose = c("None", "Low", "Medium"),
                   dd6 = tapply(d$dd6, d$hep_dose, mean, na.rm = TRUE),
                   n = as.vector(table(d$hep_dose)))
d$dose_num <- as.integer(d$hep_dose) - 1
dose_trend_p <- summary(glm(update(f_adj("dd6", "asp", ""), . ~ . + dose_num), data = d, family = binomial))$coefficients["dose_num", 4]

# ---- survival to 6 months ---------------------------------------------------------------------------
days <- c(0, 7, 14, 30, 60, 90, 120, 150, 183)
km <- list(aspirin = km_table(d$time, d$died, d$asp, days), heparin = km_table(d$time, d$died, d$hep, days))
fc <- function(arm, other) as.formula(paste("Surv(time, died) ~", paste(c(arm, other, COV), collapse = " + ")))
cox_res <- list(aspirin = cox(fc("asp", "hep"), d, "asp"), heparin = cox(fc("hep", "asp"), d, "hep"))

# ---- subgroups, with Holm correction across the pre specified factors ------------------------------
d$age_g <- cut(d$age, c(0, 69, 79, 120), labels = c("Under 70", "70 to 79", "80 and over"))
d$delay_g <- cut(d$delay, c(0, 12, 24, 48), labels = c("Up to 12 h", "13 to 24 h", "25 to 48 h"))
d$sbp_g <- cut(d$sbp, c(0, 149, 179, 400), labels = c("Under 150", "150 to 179", "180 and over"))
d$ct_g <- factor(ifelse(d$ct_done == 1, "CT before randomisation", "No CT before"))
groups <- c(age_g = "Age", sex = "Sex", consc = "Consciousness", stype = "Stroke type", af = "Atrial fibrillation",
            delay_g = "Delay", sbp_g = "Systolic BP", ct_g = "CT scan")
adj_simple <- c("age", "sex", "consc", "stype")
run_sub <- function(arm, other) {
  s <- do.call(rbind, lapply(names(groups), function(g) {
    cv <- setdiff(adj_simple, sub("_g$", "", g))
    if (g == "age_g") cv <- setdiff(cv, "age")
    x <- subgroup(d, "dd6", arm, g, c(other, cv))
    x$label <- groups[[g]]
    x
  }))
  p <- unique(s[, c("subgroup", "p_interaction")])
  p$p_holm <- holm(p$p_interaction)
  merge(s, p[, c("subgroup", "p_holm")], by = "subgroup", sort = FALSE)
}
subs <- list(aspirin = run_sub("asp", "hep"), heparin = run_sub("hep", "asp"))
write.csv(subs$aspirin, "outputs/subgroups_aspirin.csv", row.names = FALSE)
write.csv(subs$heparin, "outputs/subgroups_heparin.csv", row.names = FALSE)

# ---- sensitivity analyses on the aspirin primary result -----------------------------------------------
fa <- f_adj("dd6", "asp", "hep")
mi_vars <- c("dd6", "asp", "hep", "age", "sex", "consc", "stype", "delay", "sbp", "af_c", "infarct_visible",
             "ct_done", "death14", "stroke_death14")
mi_d <- d[, mi_vars]
mi_d$dd6 <- factor(mi_d$dd6)
meth <- make.method(mi_d); meth[meth != ""] <- "logreg"
imp <- mice(mi_d, m = 20, method = meth, printFlag = FALSE, seed = 7)
num <- function(x) { x$dd6 <- as.integer(as.character(x$dd6)); x }
mi_rd <- sapply(1:20, function(i) standardised_rd(fa, num(complete(imp, i)), "asp"))
mi_fit <- with(imp, glm(dd6 ~ asp + hep + ns(age, 3) + sex + consc + stype + ns(delay, 3) + ns(sbp, 3) + af_c +
                          infarct_visible + ct_done, family = binomial))
mi_pool <- summary(pool(mi_fit))
mi_p <- mi_pool$p.value[mi_pool$term == "asp"]
mi_or <- exp(mi_pool$estimate[mi_pool$term == "asp"] + c(0, -1.96, 1.96) * mi_pool$std.error[mi_pool$term == "asp"])
worst <- d; worst$dd6[is.na(worst$dd6)] <- ifelse(worst$asp[is.na(worst$dd6)] == 1, 1L, 0L)
best <- d; best$dd6[is.na(best$dd6)] <- ifelse(best$asp[is.na(best$dd6)] == 1, 0L, 1L)
sens <- list(
  list(analysis = "Main: complete cases, adjusted", rd = primary$aspirin$adjusted$rd, low = primary$aspirin$adjusted$low, high = primary$aspirin$adjusted$high, p = primary$aspirin$adjusted$p, n = primary$aspirin$adjusted$n),
  list(analysis = "Missing outcomes imputed (20 imputations)", rd = mean(mi_rd), low = NA, high = NA, p = mi_p, n = nrow(d),
       odds_ratio = mi_or[1], or_low = mi_or[2], or_high = mi_or[3]),
  c(list(analysis = "Missing counted against aspirin"), risk_difference(worst$dd6, worst$asp)[c("rd", "low", "high", "p")], n = nrow(d)),
  c(list(analysis = "Missing counted in favour of aspirin"), risk_difference(best$dd6, best$asp)[c("rd", "low", "high", "p")], n = nrow(d)),
  with(adjusted_rd(fa, d[d$compliant_asp, ], "asp"), list(analysis = "Per protocol: compliant with aspirin allocation", rd = rd, low = low, high = high, p = p, n = n)),
  with(adjusted_rd(fa, d[!d$not_stroke, ], "asp"), list(analysis = "Excluding final diagnosis of not a stroke", rd = rd, low = low, high = high, p = p, n = n)),
  with(adjusted_rd(fa, d[d$ct_done == 1, ], "asp"), list(analysis = "Only patients scanned before randomisation", rd = rd, low = low, high = high, p = p, n = n))
)
mi_between <- var(mi_rd)

# ---- prognostic model for the primary outcome ------------------------------------------------------
prog <- cv_prognostic(as.formula(paste("dd6 ~", paste(COV, collapse = " + "))), d)

# ---- write -------------------------------------------------------------------------------------------
results <- list(
  checks = checks,
  n_analysed_6m = sum(!is.na(d$dd6)),
  baseline = baseline,
  primary = primary, dead6 = dead6, factorial_interaction_p = factorial_p,
  events14 = events14,
  heparin_dose = list(table = dose, trend_p = dose_trend_p),
  km = km, cox = cox_res,
  subgroups = subs,
  sensitivity = sens, mi_between_imputation_sd = sqrt(mi_between),
  prognostic = prog,
  published = list(aspirin_dd6 = c(62.2, 63.5), source = "IST Collaborative Group, Lancet 1997;349:1569-81")
)
write_json(results, "docs/data/results.json", auto_unbox = TRUE, digits = 6, pretty = TRUE, na = "null")

cat(sprintf("Aspirin: dead or dependent %.1f%% vs %.1f%%, crude RD %.2f (%.2f to %.2f) per 100, p = %.3f\n",
            100 * primary$aspirin$crude$risk_treated, 100 * primary$aspirin$crude$risk_control,
            100 * primary$aspirin$crude$rd, 100 * primary$aspirin$crude$low, 100 * primary$aspirin$crude$high, primary$aspirin$crude$p))
cat(sprintf("Aspirin adjusted RD %.2f (%.2f to %.2f) per 100, OR %.3f, p = %.4f\n", 100 * primary$aspirin$adjusted$rd,
            100 * primary$aspirin$adjusted$low, 100 * primary$aspirin$adjusted$high, primary$aspirin$adjusted$odds_ratio, primary$aspirin$adjusted$p))
cat(sprintf("Heparin: %.1f%% vs %.1f%%, adjusted RD %.2f (%.2f to %.2f), p = %.3f\n", 100 * primary$heparin$crude$risk_treated,
            100 * primary$heparin$crude$risk_control, 100 * primary$heparin$adjusted$rd, 100 * primary$heparin$adjusted$low,
            100 * primary$heparin$adjusted$high, primary$heparin$adjusted$p))
cat("Factorial interaction p =", signif(factorial_p, 3), "\n")
print(events14[, c("comparison", "event", "treated", "control", "per_1000", "low", "high", "p")], digits = 3)
print(dose); cat("dose trend p", dose_trend_p, "\n")
str(cox_res)
print(unique(subs$aspirin[, c("label", "p_interaction", "p_holm")]), digits = 3)
print(unique(subs$heparin[, c("label", "p_interaction", "p_holm")]), digits = 3)
print(data.table::rbindlist(sens, fill = TRUE), digits = 3)
str(prog[1:4])
