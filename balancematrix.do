clear all

// Set root directory
if "`c(os)'" == "MacOSX" cd "/Users/alvaro/Library/Application Support/Stata/ado/personal/rddsga/"
if "`c(os)'" == "Windows" cd "C:\ado\personal\rddsga"

// Set adopath and discard loaded programs
*adopath + old_ados/balancepscore
discard

// Load dataset
use data/procurement_test2, clear
*use data/midRiskRDDdataset.dta, clear
// Run command
rddsga ///
high_direct p1 p2 p3 size_PRE1 size_PRE2 audited dTR1-dTR3 Year2 Year1 Zone1-Zone3 /// 
if (dis_cutoff2>-4 & dis_cutoff2<4), ///
  psweight(weight4) pscore(ps_flexmodel41) comsup(soporte) logit ///
  namgroup(Low share/High share)

*-------------------------------------------------------------------------------
*-------------------------------------------------------------------------------

capture program drop balancematrix
program define balancematrix, rclass
syntax varlist `touse', matname()
// Count observations in each treatment group
qui count if `touse' & `comsup' & `treatvar'==0
local Ncontrols = `r(N)'
qui count if `touse' & `comsup' & `treatvar'==1
local Ntreated = `r(N)'

// Compute propensity score weighting vector 
qui gen `psweight' = ///
  `Ntreated'/(`Ntreated'+`Ncontrols')/`pscore'*(`treatvar'==1) + ///
  `Ncontrols'/(`Ntreated'+`Ncontrols')/(1-`pscore')*(`treatvar'==0) ///
  if `touse' & `comsup' 

* Compute stats for each covariate 
*-------------------------------------------------------------------------------

local j = 0
foreach var of varlist `covariates' {
  local j=`j'+1

  // Compute and store conditional expectations
  qui reg `var' `T0' `treatvar' [iw=`psweight'] if `touse' & `comsup', noconstant
  local coef`j'_T0 = _b[`T0']
  local coef`j'_T1 = _b[`treatvar']

  // Compute and store mean differences and their p-values
  qui reg `var' `T0' [iw=`psweight'] if `touse' & `comsup'
  matrix m = r(table)
  scalar diff`j'=m[1,1] // mean difference
  local pval`j' = m[4,1] // p-value 

  // Standardized mean difference
  qui summ `var' if `touse' & `comsup'
  local stddiff`j' = (diff`j')/r(sd)
}

* Global stats
*-------------------------------------------------------------------------------

// Mean of absolute standardized mean differences (ie. stddiff + ... + stddiff`k')
/* todo: this begs to be vectorized */
local totaldiff = 0
forvalues j = 1/`numcov' {
  local totaldiff = abs(`stddiff`j'') + `totaldiff' // sum over `j' (covariates)
}
local totaldiff = `totaldiff'/`numcov' // compute mean 

// F-statistic and global p-value
qui reg `varlist' [iw=`psweight'] if `touse' & `comsup' 
local Fstat = e(F)
local pval_global = 1-F(e(df_m),e(df_r),e(F))

di in ye       "**************************************************** "
di in ye       "Propensity score-psweighting "
di in ye       "**************************************************** "

tempname `matname'
matrix `matname' = J(`numcov'+4, 4, .)
matrix colnames `matname' = "Mean `G0'" "Mean `G1'" "StMeanDiff" p-value
matrix rownames `matname' = `covariates' Observations Abs(StMeanDiff) F-statistic p-value

forvalues j = 1/`numcov' {
  matrix `matname'[`j',1] = round(`coef`j'_T0', 10^(-`bdec'))  
  matrix `matname'[`j',2] = round(`coef`j'_T1', 10^(-`bdec'))  
  matrix `matname'[`j',3] = round(`stddiff`j'', 10^(-`bdec'))
  matrix `matname'[`j',4] = round(`pval`j'', 10^(-`bdec')) 
}

matrix `matname'[`numcov'+1,1] = `Ncontrols'
matrix `matname'[`numcov'+1,2] = `Ntreated'
matrix `matname'[`numcov'+2,3] = round(`totaldiff', 10^(-`bdec'))
matrix `matname'[`numcov'+3,4] = round(`Fstat', 10^(-`bdec'))        
matrix `matname'[`numcov'+4,4] = round(`pval_global', 10^(-`bdec'))      

matrix list `matname', noheader
return matrix baltab1 = `matname'

end
