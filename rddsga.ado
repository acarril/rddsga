*! 0.3 Alvaro Carril 19jul2017
program define rddsga, rclass
version 11.1 /* todo: check if this is the real minimum */
syntax varlist(min=2 numeric) [if] [in] [ , ///
	psweight(name) pscore(name) comsup(name) balvars(varlist numeric) showbalance logit /// balancepscore opts
	cutoff(name)  /// rddsga opts
]

*-------------------------------------------------------------------------------
* Check inputs
*-------------------------------------------------------------------------------

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

// Mark observations to be used
marksample touse, novarlist

// Define model to fit (probit is default)
if "`logit'" != "" local binarymodel logit
else local binarymodel probit

// Extract treatment variable and create complementary t0 tempvar
local treatvar :	word 2 of `varlist'
tempvar t0
qui gen `t0' = (`treatvar' == 0) if !mi(`treatvar')

// Extract balance variables
if "`balvars'" == "" local balvars `varlist'
local balvars : list balvars - treatvar // remove treatvar if present
local n_balvars `: word count `balvars''

*-------------------------------------------------------------------------------
* Compute balance table matrices
*-------------------------------------------------------------------------------

* Original balance
*-------------------------------------------------------------------------------
balancematrix, matname(oribal)  ///
  touse(`touse') balvars(`balvars') ///
  treatvar(`treatvar') t0(`t0') n_balvars(`n_balvars')
return add

// Display balance matrix and global stats
if "`showbalance'" != "" {
  matlist oribal, border(rows) format(%9.3g) title("Original balance:")
  di "Obs. in T0: " oribal_Ncontrols
  di "Obs. in T1: " oribal_Ntreated
  di "Mean abs(std_diff) = " oribal_avgdiff
  di "F-statistic: " oribal_Fstat
  di "Global p-value: " oribal_pval_global
}

* Propensity Score Weighting balance
*-------------------------------------------------------------------------------
balancematrix, matname(pswbal)  ///
  touse(`touse') balvars(`balvars') ///
  psw psweight(`psweight') pscore(`pscore') comsup(`comsup') binarymodel(`binarymodel') ///
	treatvar(`treatvar') t0(`t0') n_balvars(`n_balvars')
return add

// Display balance matrix and global stats
if "`showbalance'" != "" {
  matlist pswbal, border(rows) format(%9.3g) title("Propensity Score Weighting balance:")
  di "Obs. in T0: " pswbal_Ncontrols
  di "Obs. in T1: " pswbal_Ntreated
  di "Mean abs(std_diff) = " pswbal_avgdiff
  di "F-statistic: " pswbal_Fstat
  di "Global p-value: " pswbal_pval_global
}

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
  touse(name) balvars(varlist) /// unchanging inputs
  [psw psweight(name) pscore(name) comsup(name) binarymodel(string)] /// only needed for PSW balance
	treatvar(name) t0(name) n_balvars(int) // todo: eliminate these? can be computed by subroutine at low cost

* Create variables specific to PSW matrix
*-------------------------------------------------------------------------------
if "`psw'" != "" { // if psw
  // Fit binary response model
  qui `binarymodel' `treatvar' `balvars' if `touse'

  // Generate pscore variable and clear stored results
  qui predict `pscore' if `touse'
  ereturn clear

  // Genterate common support varible
  capture drop `comsup'
  if "`comsup'" != "" {
    qui sum `pscore' if `treatvar' == 1
    qui gen `comsup' = ///
      (`pscore' >= `r(min)' & ///
       `pscore' <= `r(max)') if `touse'
    label var `comsup' "Dummy for obs. in common support"
  }
  else qui gen `comsup' == 1 if `touse'

  // Count observations in each treatment group
  qui count if `touse' & `comsup' & `treatvar'==0
  local Ncontrols = `r(N)'
  qui count if `touse' & `comsup' & `treatvar'==1
  local Ntreated = `r(N)'

  // Compute propensity score weighting vector
  cap drop `psweight'
  qui gen `psweight' = ///
    `Ntreated'/(`Ntreated'+`Ncontrols')/`pscore'*(`treatvar'==1) + ///
    `Ncontrols'/(`Ntreated'+`Ncontrols')/(1-`pscore')*(`treatvar'==0) ///
    if `touse' & `comsup' 
} // end if psw



* Count obs. in each treatment group if not PSW matrix
*-------------------------------------------------------------------------------
else { // if nopsw
  qui count if `touse' & `treatvar'==0
  local Ncontrols = `r(N)'
  qui count if `touse' & `treatvar'==1
  local Ntreated = `r(N)'
} // end if nopsw

* Compute stats specific for each covariate 
*-------------------------------------------------------------------------------
local j = 0
foreach var of varlist `balvars' {
  local ++j

  // Compute and store conditional expectations
  if "`psw'" == "" qui reg `var' `t0' `treatvar' if `touse', noconstant /* */
  else qui reg `var' `t0' `treatvar' [iw=`psweight'] if `touse' & `comsup', noconstant
  local coef`j'_T0 = _b[`t0']
  local coef`j'_T1 = _b[`treatvar']

  // Compute and store mean differences and their p-values
  if "`psw'" == "" qui reg `var' `t0' if `touse'
  else qui reg `var' `t0' [iw=`psweight'] if `touse' & `comsup'
  matrix m = r(table)
  scalar diff`j'=m[1,1] // mean difference
  local pval`j' = m[4,1] // p-value 

  // Standardized mean difference
  if "`psw'" == "" qui summ `var' if `touse'
  else qui summ `var' if `touse' & `comsup'
  local stddiff`j' = (diff`j')/r(sd)
}

* Compute global stats
*-------------------------------------------------------------------------------
// Mean of absolute standardized mean differences (ie. stddiff + ... + stddiff`k')
/* todo: this begs to be vectorized */
local avgdiff = 0
forvalues j = 1/`n_balvars' {
  local avgdiff = abs(`stddiff`j'') + `avgdiff' // sum over `j' (balvars)
}
local avgdiff = `avgdiff'/`n_balvars' // compute mean 

// F-statistic and global p-value
if "`psw'" == "" qui reg `treatvar' `balvars' if `touse'
else qui reg `treatvar' `balvars' [iw=`psweight'] if `touse' & `comsup' 
local Fstat = e(F)
local pval_global = 1-F(e(df_m),e(df_r),e(F))

* Create balance matrix
*-------------------------------------------------------------------------------
// Matrix parameters
matrix `matname' = J(`n_balvars', 4, .)
matrix colnames `matname' = mean_T0 mean_T1 std_diff p-value
matrix rownames `matname' = `balvars'

// Add per-covariate values 
forvalues j = 1/`n_balvars' {
  matrix `matname'[`j',1] = `coef`j'_T0'
  matrix `matname'[`j',2] = `coef`j'_T1'
  matrix `matname'[`j',3] = `stddiff`j''
  matrix `matname'[`j',4] = `pval`j''
}

// Return matrix and other scalars
scalar `matname'_Ncontrols = `Ncontrols'
scalar `matname'_Ntreated = `Ntreated'
scalar `matname'_avgdiff = `avgdiff'
scalar `matname'_Fstat = `Fstat'
scalar `matname'_pval_global = `pval_global'

return matrix `matname' = `matname', copy
return scalar `matname'_avgdiff = `avgdiff'
return scalar `matname'_Fstat = `Fstat'
return scalar `matname'_pvalue = `pval_global'
return scalar `matname'_N_T1 = `Ntreated'
return scalar `matname'_N_T0 = `Ncontrols'

end

********************************************************************************

/* 
CHANGE LOG
0.3
  - Integrate with rddsga.ado v2.0 original program
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
      ado includes obs. with missing values in depvar (and balvars?)

TODOS AND IDEAS:
  - Create subroutine of matlist formatting for display of balancematrix output
  - Implement matrix manipulation in Mata
  - Get rid of t0 hack for control units
*/
