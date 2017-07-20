*! 0.3 Alvaro Carril 20jul2017
program define rddsga, rclass
version 11.1 /* todo: check if this is the real minimum */
syntax varlist(min=2 numeric fv) [if] [in] , [ ///
  subgroup(name) treatment(name) /// importan inputs
	psweight(name) pscore(name) comsup(name) /// newvars
  balance(varlist numeric) showbalance logit /// balancepscore opts
	BWidth(real 0) cutoff(real 0) ///
  vce(string) ///
]

// Mark observations to be used
marksample touse, novarlist

*-------------------------------------------------------------------------------
* Check inputs
*-------------------------------------------------------------------------------

// Check that depvar and assignvar are not factor variables
local fvops = "`s(fvops)'" == "true" | _caller() >= 11 
if `fvops' { 
  local vv: di "version " ///
  string(max(11,_caller())) ", missing: " 
  gettoken first rest : varlist
  gettoken second rest : rest
  _fv_check_depvar `first'
  capture _fv_check_depvar `second'
  if _rc!=0 {
    di as error "assignvar may not be a factor variable"
    exit 198
  }
}

// psweight(): define new propensity score weighting variable or use a tempvar
if "`psweight'" != "" confirm new variable `psweight'
else tempvar psweight

// comsup(): define new common support variable or use a tempvar
if "`comsup'" != "" confirm new variable `comsup'
else tempvar comsup

// pscore(): define new propensity score variable or use a tempvar
if "`pscore'" != "" confirm new variable `pscore'
else tempvar pscore

*-------------------------------------------------------------------------------
* Process inputs
*-------------------------------------------------------------------------------

// Extract outcome variable
local depvar : word 1 of `varlist'

// Extract assignment variable
local assignvar :	word 2 of `varlist'

// Define covariates list
local covariates : list varlist - depvar
local covariates : list covariates - assignvar

// Add c. stub to continuous covariates for factor interactions
foreach var in `covariates' {
  capture _fv_check_depvar `var'
  if _rc != 0 local fv_covariates `fv_covariates' `var'
  else local fv_covariates `fv_covariates' c.`var'
}

// Create complementary subgroup var
tempvar subgroup0
qui gen `subgroup0' = (`subgroup' == 0) if !mi(`subgroup')

// Extract balance variables
if "`balance'" == "" local balance `covariates'
local n_balance `: word count `balance''

// Define model to fit (probit is default)
if "`logit'" != "" local binarymodel logit
else local binarymodel probit

// Create bandwidth condition 
local bwidth abs(`assignvar') < `bwidth'

// Create indicator cutoff variable
tempvar cutoffvar
gen `cutoffvar' = (`assignvar'>`cutoff')

*-------------------------------------------------------------------------------
* Compute balance table matrices
*-------------------------------------------------------------------------------

* Original balance
*-------------------------------------------------------------------------------
balancematrix, matname(oribal)  ///
  touse(`touse') bwidth(`bwidth') balance(`balance') ///
  subgroup(`subgroup') subgroup0(`subgroup0') n_balance(`n_balance')
return add

// Display balance matrix and global stats
if "`showbalance'" != "" {
  matlist oribal, border(rows) format(%9.3g) title("Original balance:")
  di "Obs. in subgroup 0: " oribal_N_G0
  di "Obs. in subgroup 1: " oribal_N_G1
  di "Mean abs(std_diff): " oribal_avgdiff
  di "F-statistic: " oribal_Fstat
  di "Global p-value: " oribal_pval_global
}

* Propensity Score Weighting balance
*-------------------------------------------------------------------------------
balancematrix, matname(pswbal)  ///
  touse(`touse') bwidth(`bwidth') balance(`balance') ///
  psw psweight(`psweight') pscore(`pscore') comsup(`comsup') binarymodel(`binarymodel') ///
	subgroup(`subgroup') subgroup0(`subgroup0') n_balance(`n_balance')
return add

// Display balance matrix and global stats
if "`showbalance'" != "" {
  matlist pswbal, border(rows) format(%9.3g) title("Propensity Score Weighting balance:")
  di "Obs. in subgroup 0: " pswbal_N_G0
  di "Obs. in subgroup 1: " pswbal_N_G1
  di "Mean abs(std_diff): " pswbal_avgdiff
  di "F-statistic: " pswbal_Fstat
  di "Global p-value: " pswbal_pval_global
}

*-------------------------------------------------------------------------------
* rddsga
*-------------------------------------------------------------------------------

// IVREG
/*
ivregress 2sls `depvar' i.`subgroup'#(`fv_covariates' i.gpaoXuceXr c.`assignvar' c.`assignvar'#`cutoffvar') ///
  (i.`subgroup'#1.`treatment' = i.`subgroup'#`cutoffvar') ///
  if `bwidth', vce(`vce') noconstant
*/
/*
*reg `x' `Z1' `Z0' `C`S`i''' `FE'  if `X'>-(`bw1') & `X'<(`bw1'), vce(cluster gpaoXuceXrk)
*reg I_CURaudit `Z1' `Z0' `C`S`i''' `FE'  if -(`bw1')<`assignvar' & `assignvar'<(`bw1'), vce(cluster gpaoXuceXrk)
*/

* Clear any ereturn results and end main program
*-------------------------------------------------------------------------------
ereturn clear
end

*===============================================================================
* Define auxiliary subroutines
*===============================================================================

*-------------------------------------------------------------------------------
* balancematrix: compute balance table matrices and other statistics
*-------------------------------------------------------------------------------
program define balancematrix, rclass
syntax, matname(string) /// important inputs, differ by call
  touse(name) bwidth(string) balance(varlist) /// unchanging inputs
  [psw psweight(name) pscore(name) comsup(name) binarymodel(string)] /// only needed for PSW balance
  subgroup(name) subgroup0(name) n_balance(int) // todo: eliminate these? can be computed by subroutine at low cost

* Create variables specific to PSW matrix
*-------------------------------------------------------------------------------
if "`psw'" != "" { // if psw
  // Fit binary response model
  qui `binarymodel' `subgroup' `balance' if `touse' & `bwidth'

  // Generate pscore variable and clear stored results
  qui predict `pscore' if `touse' & `bwidth'
  ereturn clear

  // Genterate common support varible
  capture drop `comsup'
  if "`comsup'" != "" {
    qui sum `pscore' if `subgroup' == 1
    qui gen `comsup' = ///
      (`pscore' >= `r(min)' & ///
       `pscore' <= `r(max)') if `touse' & `bwidth'
    label var `comsup' "Dummy for obs. in common support"
  }
  else qui gen `comsup' == 1 if `touse' & `bwidth'

  // Count observations in each treatment group
  qui count if `touse' & `bwidth' & `comsup' & `subgroup'==0
  local N_G0 = `r(N)'
  qui count if `touse' & `bwidth' & `comsup' & `subgroup'==1
  local N_G1 = `r(N)'

  // Compute propensity score weighting vector
  cap drop `psweight'
  qui gen `psweight' = ///
    `N_G1'/(`N_G1'+`N_G0')/`pscore'*(`subgroup'==1) + ///
    `N_G0'/(`N_G1'+`N_G0')/(1-`pscore')*(`subgroup'==0) ///
    if `touse' & `bwidth' & `comsup' 
} // end if psw

* Count obs. in each treatment group if not PSW matrix
*-------------------------------------------------------------------------------
else { // if nopsw
  qui count if `touse' & `bwidth' & `subgroup'==0
  local N_G0 = `r(N)'
  qui count if `touse' & `bwidth' & `subgroup'==1
  local N_G1 = `r(N)'
} // end if nopsw

* Compute stats specific for each covariate 
*-------------------------------------------------------------------------------
local j = 0
foreach var of varlist `balance' {
  local ++j

  // Compute and store conditional expectations
  if "`psw'" == "" qui reg `var' `subgroup0' `subgroup' if `touse' & `bwidth', noconstant /* */
  else qui reg `var' `subgroup0' `subgroup' [iw=`psweight'] if `touse' & `bwidth' & `comsup', noconstant
  local coef`j'_G0 = _b[`subgroup0']
  local coef`j'_G1 = _b[`subgroup']

  // Compute and store mean differences and their p-values
  if "`psw'" == "" qui reg `var' `subgroup0' if `touse' & `bwidth'
  else qui reg `var' `subgroup0' [iw=`psweight'] if `touse' & `bwidth' & `comsup'
  matrix m = r(table)
  scalar diff`j'=m[1,1] // mean difference
  local pval`j' = m[4,1] // p-value 

  // Standardized mean difference
  if "`psw'" == "" qui summ `var' if `touse' & `bwidth'
  else qui summ `var' if `touse' & `bwidth' & `comsup'
  local stddiff`j' = (diff`j')/r(sd)
}

* Compute global stats
*-------------------------------------------------------------------------------
// Mean of absolute standardized mean differences (ie. stddiff + ... + stddiff`k')
/* todo: this begs to be vectorized */
local avgdiff = 0
forvalues j = 1/`n_balance' {
  local avgdiff = abs(`stddiff`j'') + `avgdiff' // sum over `j' (balance)
}
local avgdiff = `avgdiff'/`n_balance' // compute mean 

// F-statistic and global p-value
if "`psw'" == "" qui reg `subgroup' `balance' if `touse' & `bwidth'
else qui reg `subgroup' `balance' [iw=`psweight'] if `touse' & `bwidth' & `comsup' 
local Fstat = e(F)
local pval_global = 1-F(e(df_m),e(df_r),e(F))

* Create balance matrix
*-------------------------------------------------------------------------------
// Matrix parameters
matrix `matname' = J(`n_balance', 4, .)
matrix colnames `matname' = mean_G0 mean_G1 std_diff p-value
matrix rownames `matname' = `balance'

// Add per-covariate values 
forvalues j = 1/`n_balance' {
  matrix `matname'[`j',1] = `coef`j'_G0'
  matrix `matname'[`j',2] = `coef`j'_G1'
  matrix `matname'[`j',3] = `stddiff`j''
  matrix `matname'[`j',4] = `pval`j''
}

// Return matrix and other scalars
scalar `matname'_N_G0 = `N_G0'
scalar `matname'_N_G1 = `N_G1'
scalar `matname'_avgdiff = `avgdiff'
scalar `matname'_Fstat = `Fstat'
scalar `matname'_pval_global = `pval_global'

return matrix `matname' = `matname', copy
return scalar `matname'_avgdiff = `avgdiff'
return scalar `matname'_Fstat = `Fstat'
return scalar `matname'_pvalue = `pval_global'
return scalar `matname'_N_G1 = `N_G1'
return scalar `matname'_N_G0 = `N_G0'

end

********************************************************************************

/* 
CHANGE LOG
0.3
  - Standardize syntax to merge with original rddsga.ado
0.2
  - Implement balancematrix as separate subroutine
  - Standardize balancematrix output
0.1
	- First working version, independent of project
	- Remove any LaTeX output
	- Modify some option names and internal locals

KNOWN ISSUES/BUGS:
  - Global stats don't agree with the ones computed by original balancepscore
    ~ computed mean in differences is same; r(sd) is different, maybe due to
      differences in treatment groups? check if variable.
  - Per-covariate stats don't agree with original balancepscore
    ~ In original balance this was due to different usage of `touse'; original
      ado includes obs. with missing values in depvar (and balance?)

TODOS AND IDEAS:
  - Create subroutine of matlist formatting for display of balancematrix output
  - Implement matrix manipulation in Mata
  - Get rid of subgroup0 hack
*/
