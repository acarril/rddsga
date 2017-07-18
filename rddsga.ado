*! 0.2 Alvaro Carril 17jul2017
program define rddsga, rclass
version 11.1 /* todo: check if this is the real minimum */
syntax varlist(min=2 numeric) [if] [in] [ , ///
	psweight(name) pscore(name) comsup(name) logit ///
	namgroup(string) bdec(int 3) ///
]

*-------------------------------------------------------------------------------
* Check inputs
*-------------------------------------------------------------------------------

// psweight(): define new propensity score weighting variable or use a tempvar
if `"`psweight'"' != `""' confirm new variable `psweight'
else tempvar psweight

// comsup(): define new common support variable or use a tempvar
if `"`comsup'"' != `""' confirm new variable `comsup'
else tempvar comsup

// pscore(): define new propensity score variable or use a tempvar
if `"`pscore'"' != `""' confirm new variable `pscore'
else tempvar pscore

*-------------------------------------------------------------------------------
* Process inputs
*-------------------------------------------------------------------------------

// Mark observations to be used
marksample touse

// Define model to fit (probit is default)
if `"`logit'"' != `""' local binarymodel logit
else local binarymodel probit

// Extract treatment variable and create complementary t0 tempvar
local treatvar :	word 1 of `varlist'
tempvar t0
qui gen `t0' = (`treatvar' == 0) if !mi(`treatvar')

// Extract covariates
local covariates : list varlist - treatvar
local numcov `: word count `covariates''

// Treatment groups names
if `"`namgroup'"' != `""'  {
	local pos=strpos("`namgroup'","/")
	local G0=substr("`namgroup'",1,`pos'-1)
	local G1=substr("`namgroup'",`pos'+1,.)
}
else {
	local G0="G0"
	local G1="G1"
}

*-------------------------------------------------------------------------------
* Compute balance table matrices
*-------------------------------------------------------------------------------

* Original balance
*-------------------------------------------------------------------------------
balancematrix, matname(orbal) nopsw ///
  touse(`touse') comsup(`comsup') treatvar(`treatvar') pscore(`pscore') ///
  psweight(weight5) covariates(`covariates') treatvar(`treatvar') numcov(`numcov') ///
  t0(`t0') bdec(`bdec') binarymodel(`binarymodel') 
return add

* Propensity Score Weighting balance
*-------------------------------------------------------------------------------
balancematrix, matname(pwsbal) /// 
  touse(`touse') comsup(`comsup') treatvar(`treatvar') pscore(`pscore') ///
	psweight(weight5) covariates(`covariates') treatvar(`treatvar')	numcov(`numcov') ///
	t0(`t0') bdec(`bdec') binarymodel(`binarymodel')
return add

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
capture program drop balancematrix
program define balancematrix, rclass
syntax, matname(string) psweight(name) comsup(name) /// important inputs, differ by call
	touse(name) treatvar(name) pscore(name) covariates(varlist) bdec(int) /// unchanging inputs
	treatvar(name) t0(name) numcov(int) /// todo: eliminate these; can be computed by subroutine at low cost
  [nopsw] binarymodel(string)

* Create variables specific to PSW matrix
*-------------------------------------------------------------------------------

if "`psw'" == "" { // if psw
  // Fit binary response model
  qui `binarymodel' `treatvar' `covariates' if `touse'

  // Generate pscore variable and clear stored results
  tempvar pscore
  qui predict double `pscore' if `touse'
  label var `pscore' "Estimated propensity score"
  ereturn clear // Clear e() stored results

  // Genterate common support varible
  capture drop `comsup'
  if `"`comsup'"' != `""' {
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

* Stats for each covariate 
*-------------------------------------------------------------------------------
local j = 0
foreach var of varlist `covariates' {
  local j=`j'+1

  // Compute and store conditional expectations
  if "`psw'" != "" qui reg `var' `t0' `treatvar' if `touse', noconstant
  else qui reg `var' `t0' `treatvar' [iw=`psweight'] if `touse' & `comsup', noconstant
  local coef`j'_T0 = _b[`t0']
  local coef`j'_T1 = _b[`treatvar']

  // Compute and store mean differences and their p-values
  if "`psw'" != "" qui reg `var' `t0' if `touse'
  else qui reg `var' `t0' [iw=`psweight'] if `touse' & `comsup'
  matrix m = r(table)
  scalar diff`j'=m[1,1] // mean difference
  local pval`j' = m[4,1] // p-value 

  // Standardized mean difference
  if "`psw'" != "" qui summ `var' if `touse'
  else qui summ `var' if `touse' & `comsup'
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
if "`psw'" != "" qui reg `varlist' if `touse'
else qui reg `varlist' [iw=`psweight'] if `touse' & `comsup' 
local Fstat = e(F)
local pval_global = 1-F(e(df_m),e(df_r),e(F))

* Balance matrix
*-------------------------------------------------------------------------------

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

matrix list `matname'
return matrix `matname' = `matname'

end

********************************************************************************

/* 
CHANGE LOG
0.2
	- Implement balancematrix as separate subroutine
0.1
	- First working version, independent of project
	- Remove any LaTeX output
	- Modify some option names and internal locals

TODOS (AND IDEAS TO MAKE RDDSGA EVEN COOLER)
- Create sub-program with loop that defines balance matrices
- Implement matrix manipulation in Mata
*/
