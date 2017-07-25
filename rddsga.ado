*! 0.6 Alvaro Carril 25jul2017
program define rddsga, rclass
version 11.1
syntax varlist(min=2 numeric fv) [if] [in] , ///
  SGroup(name) BWidth(real) [ Treatment(name) Cutoff(real 0) /// important inputs
  	PSWeight(name) PSCore(name) COMsup(name) noCOMsupaux /// newvars
    BALance(varlist numeric) DIBALance probit /// balancepscore opts
    IVreg REDUCEDform FIRSTstage vce(string) QUADratic ] // model opts

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
    di as error "assignvar {bf:`second'} may not be a factor variable"
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

// Mark observations to be used
marksample touse, novarlist

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

// Create complementary sgroup var
tempvar sgroup0
qui gen `sgroup0' = (`sgroup' == 0) if !mi(`sgroup')

// Extract balance variables
if "`balance'" == "" local balance `covariates'
local n_balance `: word count `balance''

// Define model to fit (logit is default)
if "`probit'" != "" local binarymodel probit
else local binarymodel logit

// Create bandwidth condition 
local bwidthtab `bwidth'
local bwidth abs(`assignvar') < `bwidth'

// Create indicator cutoff variable
tempvar cutoffvar
gen `cutoffvar' = (`assignvar'>`cutoff')
lab var `cutoffvar' "Treatment"

// Compute spline options
if "`quadratic'" != "" {
  tempvar assignXcutoff
  gen `assignXcutoff' = `assignvar'*`cutoffvar'
  local quad c.`assignvar'#c.`assignvar' c.`assignXcutoff'#c.`assignXcutoff'
}

*-------------------------------------------------------------------------------
* Compute balance table matrices
*-------------------------------------------------------------------------------

* Original balance
*-------------------------------------------------------------------------------
balancematrix, matname(oribal)  ///
  touse(`touse') bwidth(`bwidth') balance(`balance') ///
  sgroup(`sgroup') sgroup0(`sgroup0') n_balance(`n_balance')
return add

// Display balance matrix and global stats
if "`dibalance'" != "" {
  matlist oribal, border(rows) format(%9.3g) title("Unweighted balance:")
  di "Obs. in subgroup 0: " oribal_N_G0
  di "Obs. in subgroup 1: " oribal_N_G1
  di "Mean abs(std_diff): " oribal_avgdiff
  di "F-statistic: " oribal_Fstat
  di "Global p-value: " oribal_pval_global
}

* Propensity Score Weighting balance
*-------------------------------------------------------------------------------
balancematrix, matname(pswbal)  ///
  psw psweight(`psweight') touse(`touse') bwidth(`bwidth') balance(`balance') ///
  pscore(`pscore') comsup(`comsup') comsupaux(`comsupaux') binarymodel(`binarymodel') ///
	sgroup(`sgroup') sgroup0(`sgroup0') n_balance(`n_balance') 
return add

// Display balance matrix and global stats
if "`dibalance'" != "" {
  matlist pswbal, border(rows) format(%9.3g) title("Propensity Score Weighting balance:")
  di "Obs. in subgroup 0: " pswbal_N_G0
  di "Obs. in subgroup 1: " pswbal_N_G1
  di "Mean abs(std_diff): " pswbal_avgdiff
  di "F-statistic: " pswbal_Fstat
  di "Global p-value: " pswbal_pval_global
}

*-------------------------------------------------------------------------------
* Model
*-------------------------------------------------------------------------------

// Create dummy _nl_1 variable for nlcomhack
gen _nl_1 = 1
label var _nl_1 "Difference"

// Variable and value labels

label define sgroup 0 "Subgroup 0" 1 "Subgroup 1"
label values `sgroup' sgroup

label define cutoffvar 0 "Cutoff 0" 1 "Cutoff 1"
label values `cutoffvar' cutoffvar
label variable `cutoffvar' "Cutoff"

label define treatment 0 "Control" 1 "Treated"
label values `treatment' treatment

* First stage
*-------------------------------------------------------------------------------

if "`firststage'" != "" {
  // Original
  qui reg `treatment' i.`sgroup'#1.`cutoffvar' ///
    i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#`cutoffvar' `quad') ///
    if `touse' & `bwidth', vce(`vce') noconstant
  estimates title: "Unweighted first stage"
  estimates store noW_firststage
  
  // PSW
  qui reg `treatment' i.`sgroup'#1.`cutoffvar' ///
    i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#`cutoffvar' `quad') ///
    [pw=`psweight'] if `touse' & `bwidth', vce(`vce') noconstant
  estimates title: "PSW first stage"
  estimates store PSW_firststage
  
  // Output with esttab if installed; if not, default to estimates table 
  capture which estout
  if _rc!=111 {
    qui estadd local bwidthtab `bwidthtab'
    esttab noW_firststage PSW_firststage, ///
      title("First stage:") nonumbers mtitles("Unweighted" "PSW") ///
      keep(*`sgroup'#1.`cutoffvar') b(3) label ///
      se(3) star(* 0.10 ** 0.05 *** 0.01) ///
      stats(N N_clust rmse bwidthtab, fmt(0 0 3) label(Observations Clusters RMSE Bandwidth))
  }
  else {
    estimates table noW_firststage PSW_firststage, ///
      b(%9.3g) se(%9.3g) keep(i.`sgroup'#1.`cutoffvar') ///
      stats(N) varlabel title("First stage:") fvlabel
  }
}

* Reduced form
*-------------------------------------------------------------------------------
if "`reducedform'" != "" {
  // Original
  qui reg `depvar' i.`sgroup'#1.`cutoffvar' ///
    i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#`cutoffvar' `quad') ///
    if `touse' & `bwidth', vce(`vce') noconstant
  estimates title: "Unweighted reduced form"
  estimates store noW_reducedform

  // PSW
  qui reg `depvar' i.`sgroup'#1.`cutoffvar' ///
    i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#`cutoffvar' `quad') ///
    [pw=`psweight'] if `touse' & `bwidth', vce(`vce') noconstant
  estimates title: "PSW reduced form"
  estimates store PSW_reducedform

  // Output with esttab if installed; if not, default to estimates table 
  capture which estout
  if _rc!=111 {
    qui estadd local bwidthtab `bwidthtab'
    esttab noW_reducedform PSW_reducedform, ///
      title("Reduced form:") nonumbers mtitles("Unweighted" "PSW") ///
      keep(*`sgroup'#1.`cutoffvar') b(3) label ///
      se(3) star(* 0.10 ** 0.05 *** 0.01) ///
      stats(N N_clust rmse bwidthtab, fmt(0 0 3) label(Observations Clusters RMSE Bandwidth))
  }
  else {
    estimates table noW_reducedform PSW_reducedform, ///
      b(%9.3g) se(%9.3g) keep(i.`sgroup'#1.`cutoffvar') ///
      stats(N) varlabel title("Reduced form:") fvlabel
  }
}

* Instrumental variables
*-------------------------------------------------------------------------------
if "`ivreg'" != "" {
  // Original
  qui ivregress 2sls `depvar' _nl_1 ///
    i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#`cutoffvar' `quad') ///
    (i.`sgroup'#1.`treatment' = i.`sgroup'#`cutoffvar') ///
    if `touse' & `bwidth', vce(`vce') noconstant
  nlcomhack `sgroup' `treatment'
  estimates title: "Unweighted IVREG"
  estimates store noW_ivreg
  
  // PSW
  qui ivregress 2sls `depvar' _nl_1 ///
    i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#`cutoffvar' `quad') ///
    (i.`sgroup'#1.`treatment' = i.`sgroup'#`cutoffvar') /// (exogenous = endogenous)
    [pw=`psweight'] if `touse' & `bwidth', vce(`vce') noconstant
  nlcomhack `sgroup' `treatment'
  estimates title: "PSW IVREG"
  estimates store PSW_ivreg

  // Output with esttab if installed; if not, default to estimates table 
  capture which estout
  if _rc!=111 {
    qui estadd local bwidthtab `bwidthtab'
    esttab noW_ivreg PSW_ivreg, ///
      title("IV regression:") nonumbers mtitles("Unweighted" "PSW") ///
      keep(*`sgroup'#1.`treatment' _nl_1) label ///
      b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
      stats(N N_clust rmse bwidthtab, ///
        fmt(0 0 3 0) labels(Observations Clusters RMSE Bandwidth))
  }
  else{
    estimates table noW_ivreg PSW_ivreg, ///
      b(%9.3g) se(%9.3g) keep(i.`sgroup'#1.`treatment') ///
      stats(N) varlabel title("IV regression:") fvlabel
  }
}

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
  [psw psweight(name) pscore(name) comsup(name) comsupaux(string) binarymodel(string)] /// only needed for PSW balance
  sgroup(name) sgroup0(name) n_balance(int) // todo: eliminate these? can be computed by subroutine at low cost

* Create variables specific to PSW matrix
*-------------------------------------------------------------------------------
if "`psw'" != "" { // if psw
  // Fit binary response model
  qui `binarymodel' `sgroup' `balance' if `touse' & `bwidth'

  // Generate pscore variable and clear stored results
  qui predict double `pscore' if `touse' & `bwidth' & !mi(`sgroup')
  ereturn clear

  // Compute common support area by default; if not, equal comsup to 1
  if "`comsupaux'" != "nocomsupaux" {
    qui sum `pscore' if `sgroup' == 1 /* todo: check why this is like that */
    qui gen `comsup' = ///
      (`pscore' >= `r(min)' & ///
       `pscore' <= `r(max)')
    label var `comsup' "Dummy for obs. in common support"
  }
  else qui gen `comsup' = 1 if `touse' & `bwidth' & !mi(`sgroup')

  // Count observations in each treatment group
  qui count if `touse' & `bwidth' & `comsup' & `sgroup'==0
  local N_G0 = `r(N)'
  qui count if `touse' & `bwidth' & `comsup' & `sgroup'==1
  local N_G1 = `r(N)'

  // Compute propensity score weighting vector
  cap drop `psweight'
  qui gen `psweight' = ///
    `N_G1'/(`N_G1'+`N_G0')/`pscore'*(`sgroup'==1) + ///
    `N_G0'/(`N_G1'+`N_G0')/(1-`pscore')*(`sgroup'==0) ///
    if `touse' & `bwidth' & `comsup' 
} // end if psw

* Count obs. in each treatment group if not PSW matrix
*-------------------------------------------------------------------------------
else { // if nopsw
  qui count if `touse' & `bwidth' & `sgroup'==0
  local N_G0 = `r(N)'
  qui count if `touse' & `bwidth' & `sgroup'==1
  local N_G1 = `r(N)'
} // end if nopsw

* Compute stats specific for each covariate 
*-------------------------------------------------------------------------------
local j = 0
foreach var of varlist `balance' {
  local ++j

  // Compute and store conditional expectations
  if "`psw'" == "" qui reg `var' `sgroup0' `sgroup' if `touse' & `bwidth', noconstant /* */
  else qui reg `var' `sgroup0' `sgroup' [iw=`psweight'] if `touse' & `bwidth' & `comsup', noconstant
  local coef`j'_G0 = _b[`sgroup0']
  local coef`j'_G1 = _b[`sgroup']

  // Compute and store mean differences and their p-values
  if "`psw'" == "" qui reg `var' `sgroup0' if `touse' & `bwidth'
  else qui reg `var' `sgroup0' [iw=`psweight'] if `touse' & `bwidth' & `comsup'
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
if "`psw'" == "" qui reg `sgroup' `balance' if `touse' & `bwidth'
else qui reg `sgroup' `balance' [iw=`psweight'] if `touse' & `bwidth' & `comsup' 
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


*-------------------------------------------------------------------------------
* nlcomhack: hack b and V matrices to inlude nlcom results
*-------------------------------------------------------------------------------
program nlcomhack, eclass
  tempname b V 
  matrix `b' = e(b)
  matrix `V' = e(V)
  qui nlcom _b[1.`1'#1.`2'] - _b[0.`1'#1.`2']
  matrix `b'[1,3] = r(b)
  matrix `V'[3,3] = r(V)
  ereturn repost b = `b' V = `V'
end

********************************************************************************

/* 
CHANGE LOG
0.6
  - Implement nlcom hack to show difference as additional coefficient 
0.5
  - Fist working version with IVREG, reduced form and first stage equations
  - Implement output reporting with estimates table and estout
  - Default binarymodel is logit
0.4
  - First working version with IVREG equation
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
  - Should we use pweights or iweights? iw don't work with ivregress.

TODOS AND IDEAS:
  - Create subroutine of matlist formatting for display of balancematrix output
  - Implement matrix manipulation in Mata
  - Get rid of sgroup0 hack
  - Allow that groupvar is not necessarily an indicator variable
  - Is it possible to allow for N subgroups?
*/
