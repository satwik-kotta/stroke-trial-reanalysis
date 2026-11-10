/*---------------------------------------------------------------------------------------------------
  International Stroke Trial: primary analyses in SAS
  The same primary analyses as analysis/run.R, for checking the R results in SAS.

  To run on SAS OnDemand for Academics (free): run scripts/get_data.sh, upload data/ist_raw.csv to
  your home folder, set the path below, and submit.
---------------------------------------------------------------------------------------------------*/
%let path = ~/ist_raw.csv;

proc import datafile="&path" out=raw dbms=csv replace;
  guessingrows=max;
run;

data ist;
  set raw;
  asp = (RXASP = "Y");
  length hep_dose $1;
  hep_dose = RXHEP;
  if hep_dose = "H" then hep_dose = "M";      /* pilot coded medium dose as H */
  hep = (hep_dose ne "N");
  if OCCODE in (1, 2) then dd6 = 1;
  else if OCCODE in (3, 4) then dd6 = 0;
  else dd6 = .;                                /* codes 0, 8 and 9: outcome not known */
  time = min(TD, 183);
  died = (DIED = 1 and TD <= 183);
  if RATRIAL not in ("Y", "N") then RATRIAL = "U";
  if STYPE = "OTH" then STYPE = "PACS";
run;

/* Check: arms, missing outcomes */
proc freq data=ist;
  tables asp*hep_dose dd6 / missing;
run;

/* Primary outcome: dead or dependent at 6 months, crude risk difference with Newcombe interval */
proc freq data=ist;
  where dd6 ne .;
  tables asp*dd6 / riskdiff(cl=newcombe column=2) chisq;
  tables hep*dd6 / riskdiff(cl=newcombe column=2) chisq;
run;

/* Adjusted for baseline prognosis, both factorial arms in one model */
proc logistic data=ist;
  where dd6 ne .;
  class SEX RCONSC STYPE RATRIAL RVISINF RCT / param=ref;
  effect s_age = spline(AGE / naturalcubic basis=tpf(noint) knotmethod=equal(2));
  effect s_delay = spline(RDELAY / naturalcubic basis=tpf(noint) knotmethod=equal(2));
  effect s_sbp = spline(RSBP / naturalcubic basis=tpf(noint) knotmethod=equal(2));
  model dd6(event="1") = asp hep s_age SEX RCONSC STYPE s_delay s_sbp RATRIAL RVISINF RCT;
  oddsratio asp;
run;

/* Factorial interaction */
proc logistic data=ist;
  where dd6 ne .;
  class SEX RCONSC STYPE RATRIAL / param=ref;
  model dd6(event="1") = asp hep asp*hep AGE SEX RCONSC STYPE RDELAY RSBP RATRIAL;
run;

/* Survival to 6 months */
proc lifetest data=ist plots=survival(cl) timelist=14 30 90 183;
  time time*died(0);
  strata asp;
run;

proc phreg data=ist;
  class SEX RCONSC STYPE RATRIAL / param=ref;
  model time*died(0) = asp hep AGE SEX RCONSC STYPE RDELAY RSBP RATRIAL;
  assess ph / resample;
  hazardratio asp;
run;
