*! 0.1 Alvaro Carril 14jul2017
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

/*******/
/* NEW */
/*******/

// Fit binary response model
capture drop comsup /* todo: don't drop automatically, user-generated name */
qui `binarymodel' `treatvar' `covariates' if `touse'

// Generate pscore variable and clear results
tempvar pscore
qui predict double `pscore' if `touse'
label var `pscore' "Estimated propensity score"
ereturn clear // Clear e() stored results

/* NEW */
/*******/

// Region of common support
if `"`comsup'"' != `""'  {
	// Genterate common support varible
	qui sum `pscore' if `treatvar' == 1
	qui gen `comsup' = (`pscore' >= `r(min)' & ///
											`pscore' <= `r(max)') if `touse'
	label var `comsup' "Dummy for obs. in common support"
}
else qui gen `comsup' == 1 if `touse'

// Count observations in each sample
qui count if `touse' & `treatvar'==0
local N0 = `r(N)'
qui count if `touse' & `treatvar'==1
local N1 = `r(N)'


*** Original balance

local j=0
foreach var of varlist `covariates' {
	local j=`j'+1

	*get mean group 1 and mean group 2
	qui reg `var'  `T0' `treatvar'  if `touse', noconstant
	mat m=r(table)
	/* mean G0  */ 	local m`j'1:  di m[1,1]
	/* mean G1 */ 	local m`j'2:  di m[1,2]


	qui reg `var'  `T0'  if `touse'
	mat m=r(table)
	scalar dif_`var'=m[1,1]
	/* p-value */ 	local m`j'4:  di  m[4,1]


	qui tabstat `var'  if `touse', stat(sd) save
	matrix overall= r(StatTotal)'

	local stdiff_`j'=(dif_`var')/overall[1,1]
	local m`j'3:  di (dif_`var')/overall[1,1]
}

*** abs sd mean difference (stdiff1+stdiff2+...+stdiffk), i.e, k shows the sum of the standard mean difference
local k=0
forvalue j=1/`numcov' {
	foreach x  of numlist `stdiff_`j'' {
		local k=abs(`x')+`k'
	}
}

** Then, we show the mean of the standard mean difference (the sum divided by the number of covariates):
local l=`numcov'+1
local m`l'3:di `k'/`numcov'

*** F statistics

qui reg `varlist' if `touse'
local l=`l'+1
local m`l'4: di e(F)

** p-value
local l=`l'+1
local m`l'4: di 1-F(e(df_m),e(df_r),e(F))

**************************
*** Propensity-score Weighting // empieza con touse y com
**************************


**Observations:
preserve
qui keep if  `touse'  & `comsup' & `treatvar'==0
local Npsweight0=_N
restore
preserve
qui keep if `touse' & `comsup' & `treatvar'==1
local Npsweight1=_N
restore

* compute psweights

qui sum `treatvar' if `touse' & `comsup' & `treatvar'==1
local hd=r(N)
qui sum high_direct if `touse' & `comsup' & `treatvar'==0
local ld=r(N)

*** gen the psweight for each observation using non-conditional probability and conditional probability.
qui gen `psweight' = `hd'/(`hd'+`ld')/`pscore'*(`treatvar'==1) + `ld'/(`hd'+`ld')/(1-`pscore')*(`treatvar'==0) if `touse' & `comsup' 


local j=0
foreach var of varlist `covariates' {
	local j=`j'+1

	qui reg `var' `T0' `treatvar' [iw=`psweight'] if `touse' & `comsup'  , noconstant
	mat m=r(table)

	/* Low direct  */ 	local Weight`j'_1:  di  m[1,1]
	/* High direct  */ 	local Weight`j'_2:  di  m[1,2]


	qui reg `var'  `T0' [iw=`psweight'] if `touse' & `comsup' 
	mat m=r(table)
	scalar dif_`var'=m[1,1]
	/* p-value */ 	local Weight`j'_4:  di m[4,1]

	qui tabstat `var' if `touse' & `comsup'  , stat(sd) save
	matrix overall= r(StatTotal)'
	local stdiff_`j'=(dif_`var')/overall[1,1]
	local Weight`j'_3:  di  (dif_`var')/overall[1,1]
}


*** abs sd mean difference 
local k=0
forvalue j=1/`numcov' {
	foreach x  of numlist `stdiff_`j'' {
		local k=abs(`x')+`k'
	}
}

** Compute
local l=`numcov'+1
local Weight`l'_3:di `k'/`numcov'

*** global F-STATISTIC and P-VALUE 
qui reg `varlist'  [iw=`psweight'] if `touse' & `comsup' 

local l=`l'+1
local Weight`l'_4: di  e(F)
local l=`l'+1
local Weight`l'_4: di  1-F(e(df_m),e(df_r),e(F))


di in ye       "**************************************************** "
di in ye	     "ORIGINAL BALANCE "
di in ye	     "**************************************************** "


tempname orbal
matrix `orbal' = J(`numcov'+4,4,.)

local j=0                              
foreach var of varlist `covariates' {
	local j=`j'+1  
	matrix `orbal'[`j',1] = round(`m`j'1',10^(-`bdec'))	
	matrix `orbal'[`j',2] = round(`m`j'2',10^(-`bdec'))
	matrix `orbal'[`j',3] = round(`m`j'3',10^(-`bdec'))
	matrix `orbal'[`j',4] = round(`m`j'4',10^(-`bdec'))
	local rown3 "`rown3' `var'"
}

matrix `orbal'[`numcov'+1,1] = `N0'
matrix `orbal'[`numcov'+1,2] = `N1'			
local l=`numcov'+1
matrix `orbal'[`numcov'+2,3] = round(`m`l'3',10^(-`bdec'))
local l=`l'+1		
matrix `orbal'[`numcov'+3,4] = round(`m`l'4',10^(-`bdec'))
local l=`l'+1			
matrix `orbal'[`numcov'+4,4] = round(`m`l'4',10^(-`bdec'))


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

di in ye       "**************************************************** "
di in ye	     "Propensity score-psweighting "
di in ye	     "**************************************************** "


tempname balimp
matrix `balimp' = J(`numcov'+4,4,.)

*** Complete the matrix with the values respective:
local j=0                              
foreach var of varlist `covariates' {
	local j=`j'+1  
	matrix `balimp'[`j',1] = round(`Weight`j'_1', 10^(-`bdec'))	
	matrix `balimp'[`j',2] = round(`Weight`j'_2', 10^(-`bdec'))	
	matrix `balimp'[`j',3] = round(`Weight`j'_3', 10^(-`bdec'))
	matrix `balimp'[`j',4] = round(`Weight`j'_4', 10^(-`bdec'))	
	
	local rown4 "`rown4' `var'"
}

matrix `balimp'[`numcov'+1,1] = `Npsweight0'
matrix `balimp'[`numcov'+1,2] = `Npsweight1'			
local l=`numcov'+1
matrix `balimp'[`numcov'+2,3] = round(`Weight`l'_3',10^(-`bdec'))
local l=`l'+1		
matrix `balimp'[`numcov'+3,4] = round(`Weight`l'_4',10^(-`bdec'))			
local l=`l'+1			
matrix `balimp'[`numcov'+4,4] = round(`Weight`l'_4',10^(-`bdec'))			


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


*** RETURN

eret clear 

*** Produce the comsup variable to be generated after running the command
/* if `"`comsup'"' != "" {
	qui g comsup = `COMSUP'
	label var comsup "Dummy for obs. in common support"
} */

end

********************************************************************************

/* 
CHANGE LOG
0.1
	- First working version, independent of project
	- Remove any LaTeX output
	- Modify some option names and internal locals

TODOS (AND IDEAS TO MAKE RDDSGA EVEN COOLER)
- Create sub-program with loop that defines balance matrices
*/
