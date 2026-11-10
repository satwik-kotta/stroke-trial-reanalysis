# Read the IST database, check it, and build the analysis dataset.
#
# One row per randomised patient. Treatment is analysed as randomised (intention to treat).
# Outcomes:
#   dd6      dead or dependent at 6 months (OCCODE 1 or 2), the trial's primary outcome
#   death14  dead within 14 days
#   isc14    recurrent ischaemic stroke within 14 days
#   haem14   haemorrhagic stroke within 14 days
#   bleed14  major non cerebral bleed within 14 days (transfused or fatal)
#   stroke_death14  death or non fatal recurrent stroke within 14 days

HORIZON <- 183  # days, the 6 month follow up

yn <- function(x) ifelse(x == "Y", 1L, ifelse(x == "N", 0L, NA_integer_))

data_checks <- function(raw) {
  occ <- raw$OCCODE
  list(
    patients = nrow(raw),
    columns = ncol(raw),
    duplicate_rows = sum(duplicated(raw)),
    aspirin_arms = as.list(table(raw$RXASP)),
    heparin_arms = as.list(table(raw$RXHEP)),
    pilot_high_dose_coded_H = sum(raw$RXHEP == "H"),
    outcome_6m_missing_code_8_9 = sum(occ %in% c(8, 9)),
    outcome_6m_undocumented_code_0 = sum(occ == 0),
    outcome_6m_known = sum(occ %in% 1:4),
    atrial_fibrillation_blank = sum(raw$RATRIAL == "" | is.na(raw$RATRIAL)),
    final_diagnosis_not_stroke = sum(raw$DNOSTRK == "Y", na.rm = TRUE),
    age_range = range(raw$AGE),
    delay_hours_range = range(raw$RDELAY),
    sbp_range = range(raw$RSBP),
    time_to_death_missing = sum(is.na(raw$TD)),
    died_but_flag_mismatch = sum(raw$DIED == 1 & raw$FDEAD == "N", na.rm = TRUE),
    aspirin_noncompliant = sum(raw$CMPLASP == "N"),
    heparin_noncompliant = sum(raw$CMPLHEP == "N")
  )
}

prepare <- function(raw) {
  occ <- raw$OCCODE
  d <- data.frame(
    asp = as.integer(raw$RXASP == "Y"),
    # the pilot coded medium dose as H; the database notes say H and M are the same dose
    hep_dose = factor(ifelse(raw$RXHEP == "H", "M", raw$RXHEP), levels = c("N", "L", "M")),
    age = raw$AGE,
    sex = factor(raw$SEX, levels = c("M", "F")),
    consc = factor(raw$RCONSC, levels = c("F", "D", "U"), labels = c("Alert", "Drowsy", "Unconscious")),
    stype = factor(ifelse(raw$STYPE == "OTH", "PACS", raw$STYPE), levels = c("PACS", "LACS", "POCS", "TACS")),
    delay = raw$RDELAY,
    sbp = raw$RSBP,
    af = factor(ifelse(raw$RATRIAL %in% c("Y", "N"), raw$RATRIAL, "Not recorded"), levels = c("N", "Y", "Not recorded")),
    ct_done = yn(raw$RCT),
    infarct_visible = yn(raw$RVISINF),
    asp_before = yn(raw$RASP3),
    hep_before = yn(raw$RHEP24),
    dd6 = ifelse(occ %in% 1:2, 1L, ifelse(occ %in% 3:4, 0L, NA_integer_)),
    dead6 = ifelse(occ == 1, 1L, ifelse(occ %in% 2:4, 0L, NA_integer_)),
    death14 = raw$ID14,
    isc14 = raw$ISC14,
    haem14 = raw$H14,
    bleed14 = raw$TRAN14,
    stroke_death14 = as.integer(raw$ID14 == 1 | raw$STRK14 == 1),
    time = pmin(raw$TD, HORIZON),
    died = as.integer(raw$DIED == 1 & raw$TD <= HORIZON),
    compliant_asp = raw$CMPLASP == "Y",
    not_stroke = raw$DNOSTRK == "Y",
    country = raw$CNTRYNUM
  )
  d$hep <- as.integer(d$hep_dose != "N")
  d$af[d$af == "Not recorded"] <- NA
  d$af <- droplevels(d$af)
  d
}
