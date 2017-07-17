*! 0.1 Alvaro Carril 17jul2017
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

// Extract treatment variable and create complementary T0 tempvar
local treatvar :	word 1 of `varlist'
tempvar T0
qui gen `T0' = (`treatvar' == 0) if !mi(`treatvar')

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
* Propensity Score
*-------------------------------------------------------------------------------

// Fit binary response model
capture drop comsup /* todo: don't drop automatically, user-generated name */
qui `binarymodel' `treatvar' `covariates' if `touse'

// Generate pscore variable and clear stored results
tempvar pscore
qui predict double `pscore' if `touse'
label var `pscore' "Estimated propensity score"
ereturn clear // Clear e() stored results

*-------------------------------------------------------------------------------
* Common Support
*-------------------------------------------------------------------------------

// Genterate common support varible
if `"`comsup'"' != `""' {
	qui sum `pscore' if `treatvar' == 1
	qui gen `comsup' = ///
		(`pscore' >= `r(min)' & ///
		 `pscore' <= `r(max)') if `touse'
	label var `comsup' "Dummy for obs. in common support"
}
else qui gen `comsup' == 1 if `touse'

*-------------------------------------------------------------------------------
* Original Balance
*-------------------------------------------------------------------------------

// Count observations in each treatment group
qui count if `touse' & `treatvar'==0
local Ncontrols = `r(N)'
qui count if `touse' & `treatvar'==1
local Ntreated = `r(N)'

* Compute stats for each covariate 
*-------------------------------------------------------------------------------

local j=0
foreach var of varlist `covariates' {
	local j=`j'+1

	// Compute and store conditional expectations
	qui reg `var' `T0' `treatvar' if `touse', noconstant
	local coef`j'_T0 = _b[`T0']
	local coef`j'_T1 = _b[`treatvar']

	// Compute and store mean differences and their p-values
	qui reg `var' `T0' if `touse'
	matrix m=r(table)
	scalar diff`j' = m[1,1] // mean difference
	local pval`j' = m[4,1] // p-value 

	// Standardized mean difference
	qui summ `var' if `touse'
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
qui reg `varlist' if `touse'
local Fstat = e(F)
local pval_global = 1-F(e(df_m),e(df_r),e(F))

di in ye       "**************************************************** "
di in ye	     "ORIGINAL BALANCE "
di in ye	     "**************************************************** "


tempname orbal
matrix `orbal' = J(`numcov'+4,4,.)

local j=0                              
foreach var of varlist `covariates' {
	local j=`j'+1  
	matrix `orbal'[`j',1] = round(`coef`j'_T0', 10^(-`bdec'))	
	matrix `orbal'[`j',2] = round(`coef`j'_T1', 10^(-`bdec'))
	matrix `orbal'[`j',3] = round(`stddiff`j'', 10^(-`bdec'))
	matrix `orbal'[`j',4] = round(`pval`j'', 10^(-`bdec'))
	local rown3 "`rown3' `var'"
}

matrix `orbal'[`numcov'+1,1] = `Ncontrols'
matrix `orbal'[`numcov'+1,2] = `Ntreated'
matrix `orbal'[`numcov'+2,3] = round(`totaldiff', 10^(-`bdec'))
matrix `orbal'[`numcov'+3,4] = round(`Fstat', 10^(-`bdec'))
matrix `orbal'[`numcov'+4,4] = round(`pval_global', 10^(-`bdec'))

matrix colnames `orbal' = "Mean `G0'" "Mean `G1'" "StMeanDiff" p-value 
matrix rownames `orbal' = `rown3' Observations Abs(StMeanDiff) F-statistic p-value

local form ", noheader"
*XXX RD: where is format coming from? how does one specify it as non-missing? XXX
if "`format'" != "" {
	local form "`form' `format'"
}
matrix list `orbal' `form'

if "`matrix'" != "" {
	matrix `matrix' = `orbal'
}

return matrix orbal = `orbal'

*-------------------------------------------------------------------------------
* Propensity Score Weighting
*-------------------------------------------------------------------------------

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
di in ye	     "Propensity score-psweighting "
di in ye	     "**************************************************** "


tempname balimp
matrix `balimp' = J(`numcov'+4,4,.)

local j=0                              
foreach var of varlist `covariates' {
	local j=`j'+1  
	matrix `balimp'[`j',1] = round(`coef`j'_T0', 10^(-`bdec'))	
	matrix `balimp'[`j',2] = round(`coef`j'_T1', 10^(-`bdec'))	
	matrix `balimp'[`j',3] = round(`stddiff`j'', 10^(-`bdec'))
	matrix `balimp'[`j',4] = round(`pval`j'', 10^(-`bdec'))	
	local rown4 "`rown4' `var'"
}

matrix `balimp'[`numcov'+1,1] = `Ncontrols'
matrix `balimp'[`numcov'+1,2] = `Ntreated'
matrix `balimp'[`numcov'+2,3] = round(`totaldiff', 10^(-`bdec'))
matrix `balimp'[`numcov'+3,4] = round(`Fstat', 10^(-`bdec'))				
matrix `balimp'[`numcov'+4,4] = round(`pval_global', 10^(-`bdec'))			

matrix colnames `balimp' = "Mean `G0'" "Mean `G1'" "StMeanDiff" p-value
matrix rownames `balimp' = `rown4' Observations Abs(StMeanDiff) F-statistic p-value

local form ", noheader"
if "`format'" != "" {
	local form "`form' `format'"
}
matrix list `balimp' `form'
	
if "`matrix'" != "" {
	matrix `matrix' = `balimp'
}

return matrix balimp = `balimp'

ereturn clear

end

*-------------------------------------------------------------------------------
* Define auxiliary programs
*-------------------------------------------------------------------------------



********************************************************************************

/* 
CHANGE LOG
0.1
	- First working version, independent of project
	- Remove any LaTeX output
	- Modify some option names and internal locals

TODOS (AND IDEAS TO MAKE RDDSGA EVEN COOLER)
- Create sub-program with loop that defines balance matrices
- Implement matrix manipulation in Mata
*/
