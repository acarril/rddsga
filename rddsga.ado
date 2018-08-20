*! 1.0.3 Andre Cazor 03012018
program define rddsga, eclass
version 11.1
syntax varlist(min=2 numeric fv) [if] [in] , ///
  SGroup(name) BWidth(real) [ fuzzy(name) Cutoff(real 0) Kernel(string) /// important inputs
  	IPSWeight(name) PSCore(name) comsup /// newvars
    BALance(varlist numeric) DIBALance probit /// balancepscore opts
    IVregress REDUCEDform FIRSTstage vce(string) p(int 1) /// model opts
    noBOOTstrap bsreps(real 50) FIXEDbootstrap BLOCKbootstrap(string) NORMal  noipsw weights(string) ] // bootstrap options

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


if ("`weights'"~="") {
tempvar oweights
g `oweights'=`weights'
}
else {
tempvar oweights
g `oweights'=1
}


/*
// Check that fuzzy() is specified if ivregress is specified
if "`ivregress'" != "" & "`fuzzy'" == "" {
  di as error "fuzzy() must be specified with ivregress"
  exit 198
}*/

// ipsweight(): define new propensity score weighting variable or use a tempvar
if "`ipsweight'" != "" confirm new variable `ipsweight'
else tempvar ipsweight


// comsup(): define tempvar for no common support (default option)
if "`comsup'" == "" {
tempvar NOCOMSUP
g `NOCOMSUP'=1
}
// pscore(): define new propensity score variable or use a tempvar
if "`pscore'" != "" confirm new variable `pscore'
else tempvar pscore



// Issue warning if no covariates and no vars in balance when propensity score weighting is used:
if ("`ipsw'"!="noipsw")  {
if `: list sizeof varlist'<=2 & `: list sizeof balance'==0 {
  di as error "either {it:indepvars} or {bf:balance()} must be specified"
  exit 198
}
}

//  When Propensity score weighting is not used, request covariates only if dibalance is specified:
else {
if "`dibalance'" != ""  & `: list sizeof varlist'<=2 & `: list sizeof balance'==0 {
  di as error "either {it:indepvars} or {bf:balance()} must be specified"
  exit 198
}
}



// Issue warning if options bsreps and normal are specified along with nobootstrap
if "`bootstrap'" == "nobootstrap" & (`bsreps' != 50 | "`normal'" != "") {
  di as text "Warning: options " as result "bsreps" as text " and " as result "normal" ///
    as text " are irrelevant if " as result "nobootstrap" as text " is specified"
}




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
if "`balance'" == "" & "`ipsw'" != "noipsw"  local balance `covariates'
local n_balance `: word count `balance''

if "`balance'" == "" & "`dibalance'" != "" & "`ipsw'" == "noipsw"  local balance `covariates'
local n_balance `: word count `balance''

// Define model to fit (logit is default)
if "`probit'" != "" local binarymodel probit
else local binarymodel logit

// Create bandwidth condition 
local bwidthtab `bwidth'
local bwidth abs(`assignvar') < `bwidth'

// Create indicator cutoff variable
*tempvar cutoffvar
*gen _cutoff = (`assignvar'>`cutoff')
*lab var _cutoff "fuzzy"
confirm new variable _cutoff
gen _cutoff = (`assignvar'>`cutoff')

// Compute spline options
forval i=1/`p' {
tempvar assignvar_`i'
qui gen `assignvar_`i''=`assignvar'^`i'
tempvar gassignvar_`i'
qui gen `gassignvar_`i''=`assignvar_`i''*_cutoff
local polynm `polynm' i.`sgroup'#c.`assignvar_`i'' i.`sgroup'#c.`gassignvar_`i''
}


tempvar kwt
// create weights for kernel 

* default 
if "`kernel'"=="" {
local kernel = "uni" 
}

  if ("`kernel'"=="epanechnikov" | "`kernel'"=="epa") {
    local kernel_type = "Epanechnikov"
    qui g double `kwt'=max(0,3/4*(`bwidthtab'^2-abs(`2')^2))*`oweights'
  }
  else if ("`kernel'"=="triangular" | "`kernel'"=="tri") {
      local kernel_type = "Triangular"
      qui g double `kwt'=max(0,`bwidthtab'-abs(`2'))*`oweights'
  }
  else {
    local kernel_type = "Uniform"
    qui g double `kwt'=(-`bwidthtab'<=(`2') & `2'<`bwidthtab')*`oweights'
  }
       

// Create weight local for balance
if "`ipsw'" == "noipsw" local weight = ""
else local weight "[pw=`ipsweight']"


*** Count number of observations when noipsw is specified and there are not variables to balance:

if `: list sizeof balance'==0 {

 qui count if `touse' & `bwidth' & `sgroup'==0
 local N_G0 = `r(N)'
 qui count if `touse' & `bwidth' & `sgroup'==1
 local N_G1 = `r(N)'
 
scalar unw_N_G1 = `N_G1'
scalar unw_N_G0 = `N_G0' 


}

*-------------------------------------------------------------------------------
* Compute balance table matrices
*-------------------------------------------------------------------------------

* Original balance
*-------------------------------------------------------------------------------

if `: list sizeof balance'!=0 {

// Compute balanace matrix 
balancematrix, matname(unw)  ///
 weights(`oweights') touse(`touse') bwidth(`bwidth') balance(`balance') ///
  sgroup(`sgroup') sgroup0(`sgroup0') n_balance(`n_balance')
  

// Store balance matrix and computed balance stats
matrix unw = e(unw)
foreach s in unw_N_G0 unw_N_G1 unw_pvalue unw_Fstat unw_avgdiff {
  scalar `s' = e(`s')
}
// Display balance matrix and balance stats
if "`dibalance'" != "" {
  di _newline as result "Unweighted"
  matlist unw, ///
    border(rows) format(%9.3g) noblank
  di "Obs. in subgroup 0: " unw_N_G0
  di "Obs. in subgroup 1: " unw_N_G1
  di "Mean abs(std_diff): " unw_avgdiff
  di "F-statistic: " unw_Fstat
  di "Global p-value: " unw_pval_global
}

if "`ipsw'" != "noipsw" {

* Propensity Score Weighting balance
*-------------------------------------------------------------------------------
// Compute balanace matrix 

balancematrix, matname(ipsw)  ///
  psw ipsweight(`ipsweight') weights(`oweights') touse(`touse') bwidth(`bwidth') balance(`balance') ///
  pscore(`pscore') comsup nocomsup(`NOCOMSUP') binarymodel(`binarymodel') ///
	sgroup(`sgroup') sgroup0(`sgroup0') n_balance(`n_balance') 
	
// Store balance matrix and computed balance stats
matrix ipsw = e(ipsw)
foreach s in ipsw_N_G0 ipswN_G1 ipsw_pvalue ipsw_Fstat ipsw_avgdiff {
  scalar `s' = e(`s')
}
// Display balance matrix and balance stats
if "`dibalance'" != "" {
  di _newline as result "Inverse Propensity Score Weighting"
  matlist ipsw, ///
  border(rows) format(%9.3g) noblank
  di "Obs. in subgroup 0: " ipsw_N_G0
  di "Obs. in subgroup 1: " ipsw_N_G1
  di "Mean abs(std_diff): " ipsw_avgdiff
  di "F-statistic: " ipsw_Fstat
  di "Global p-value: " ipsw_pval_global
}
}
}



*-------------------------------------------------------------------------------
* Estimation
*-------------------------------------------------------------------------------


tempvar kernelipsw


// Create weight local for regression
if "`ipsw'" == "noipsw" { 
local weight = "[pw=`kwt']"
}
else { 
qui g double `kernelipsw'=`ipsweight'*`kwt'
local weight "[pw=`kernelipsw']"
}




* First stage
*-------------------------------------------------------------------------------
if "`firststage'" != "" {
    mat IndIV=[0]
  // Regression
    qui reg `fuzzy' i.`sgroup'#1._cutoff i.`sgroup' ///
    i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#_cutoff)  `polynm' ///
    `weight' if `touse' & `bwidth', vce(`vce')
   
  ** Escalar used as a indicator to compute bootstrap. 
    mat fixed=[0] 

if "`blockbootstrap'"!="" {
ereturn local blockbootstrap `blockbootstrap'
}
    
  // Compute bootstrapped variance-covariance matrix and post results
  if "`bootstrap'" != "nobootstrap" myboo `sgroup' _cutoff `bsreps'
  // If no bootstrap, trim b and V to show only RD estimates
  else epost `sgroup' _cutoff 
}

* Reduced form
*-------------------------------------------------------------------------------
if "`reducedform'" != "" {
  // Regression
    mat IndIV=[0]
    qui reg `depvar' i.`sgroup'#1._cutoff i.`sgroup' ///
    i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#_cutoff) `polynm' ///
    `weight' if `touse' & `bwidth', vce(`vce')
  
  ** Escalar used as a indicator to compute bootstrap. 
    mat fixed=[0]
   
if "`blockbootstrap'"!="" {
ereturn local blockbootstrap `blockbootstrap'
}

  // Compute bootstrapped variance-covariance matrix and post results
  if "`bootstrap'" != "nobootstrap" myboo `sgroup' _cutoff `bsreps'
  // If no bootstrap, trim b and V to show only RD estimates
  else epost `sgroup' _cutoff 
}

* Instrumental variables
*-------------------------------------------------------------------------------
if "`ivregress'" != "" {
  // Regression

  mat IndIV=[1]
 
 qui reg `fuzzy' i.`sgroup'#1._cutoff i.`sgroup' ///
   i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#_cutoff) `polynm' ///
   `weight' if `touse' & `bwidth', vce(`vce')
   
   
 local coeffFSg0: di _b[0.`sgroup'#1._cutoff]
 local coeffFSg1: di _b[1.`sgroup'#1._cutoff]
 
 ** If RDD is fuzzy and we use a fixed bootstrap, so we save the first stage coefficient 
   mat FS=[`coeffFSg0',`coeffFSg1']

   qui reg `depvar' i.`sgroup'#1._cutoff i.`sgroup' ///
    i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#_cutoff) `polynm' ///
    `weight' if `touse' & `bwidth', vce(`vce')
	
local RFline `e(cmdline)'
  
 qui ivregress 2sls `depvar' i.`sgroup' ///
   i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#_cutoff) `polynm' ///
    (i.`sgroup'#1.`fuzzy' = i.`sgroup'#1._cutoff) ///
    `weight' if `touse' & `bwidth', vce(`vce')

** Escalar used as a indicator to compute bootstrap. 0 if bootstrap is computed using IV estimation
mat fixed=[0]

if "`blockbootstrap'" != "" {
ereturn local blockbootstrap `blockbootstrap'
}


if "`fixedbootstrap'" != "" {  

** If RDD is fuzzy and we use fixed bootstrap, so we use reduced form estimation to compute bootstrap in the IV
ereturn local cmdline `RFline'  

** Escalar used as a indicator to compute bootstrap. 1 if bootstrap is computed using reduced form estimation
** and keeping first stage fixed. 
 mat fixed=[1]
    }

  // Compute bootstrapped variance-covariance matrix and post results
  if "`bootstrap'" != "nobootstrap" myboo `sgroup' `fuzzy' `bsreps'  

  // If no bootstrap, trim b and V to show only RD estimates
  else epost `sgroup' `fuzzy'
*  mat cumulative = e(cumulative)
*  ereturn matrix cumulative = cumulative

}


* Post balance results
*-------------------------------------------------------------------------------
// Post global balance stats
if `: list sizeof balance'!=0 & "`ipsw'" != "noipsw"   {

foreach w in unw ipsw {
  foreach s in N_G0 N_G1 pvalue Fstat avgdiff {
    ereturn scalar `w'_`s' = `w'_`s'
  }

}
// Post balance matrices
ereturn matrix ipsw ipsw
ereturn matrix unw unw

}
if "`ipsw'" == "noipsw" & `: list sizeof balance'!=0 {

local w unw
  foreach s in N_G0 N_G1 pvalue Fstat avgdiff {
    ereturn scalar `w'_`s' = `w'_`s'
  }

// Post balance matrices
ereturn matrix unw unw
}
if "`ipsw'" == "noipsw" & `: list sizeof balance'==0  {

ereturn scalar unw_N_G1 = `N_G1'
ereturn scalar unw_N_G0 = `N_G0'

}

*-------------------------------------------------------------------------------
* Results
*-------------------------------------------------------------------------------

* Post and display estimation results
*-------------------------------------------------------------------------------
if "`ivregress'" != "" | "`reducedform'" != "" | "`firststage'" != "" {

  // Post abridged b and V matrices
  ereturn repost b=b V=V, resize
  // Display estimates by subgroup
*  di as result "Subgroup estimates"
*  ereturn display
  // Display difference of subgroup estimates 
*  di _newline as result "Difference estimate"
  if "`ivregress'" == "" {
*   di as text "_nl_1 = _b[1.`sgroup'#1._cutoff] - _b[0.`sgroup'#1._cutoff]" _continue
    qui nlcom _b[1.`sgroup'#1._cutoff] - _b[0.`sgroup'#1._cutoff]
  }
  else {
*   di as text "_nl_1 = _b[1.`sgroup'#1.`fuzzy'] - _b[0.`sgroup'#1.`fuzzy']" _continue
    qui nlcom _b[1.`sgroup'#1.`fuzzy'] - _b[0.`sgroup'#1.`fuzzy']
  } 

  * Compute and store subgroup estimates 
  *-------------------------------------------------------------------------------
  if "`ivregress'" == "" scalar df = e(df_r)
  else scalar df = e(df_m)

// Estimates by subgroup
  forvalues g = 0/1 {
    // Coefficient
    matrix e_b = e(b)
    scalar b_g`g' = e_b[1,`=`g'+1']
    // Standard error
    matrix e_V = e(V)
    scalar se_g`g' = sqrt(e_V[`=`g'+1',`=`g'+1'])
    // t-stat 
    scalar t_g`g' = b_g`g'/se_g`g'
    // P-value
    scalar p_g`g' = ttail(df, abs(t_g`g'))*2
*    scalar p_g`g' = e(p_g`g')
    // Confidence interval
    scalar ci_lb_g`g' = b_g`g' + invttail(df, 0.975)*se_g`g'
    scalar ci_ub_g`g' = b_g`g' + invttail(df, 0.025)*se_g`g'
*    scalar ci_ub_g`g'_norm = e(ci_ub_g`g')
*    scalar ci_lb_g`g' = e(ci_lb_g`g')
  }

  * Compute and store difference estimates 
  *-------------------------------------------------------------------------------
  // Coefficient
  matrix e_b_diff = r(b)
  scalar b_diff = e_b_diff[1,1]
  // Standard error
  matrix e_V_diff = r(V)
  scalar se_diff = sqrt(e_V_diff[1,1])
  // t-stat 
  scalar t_diff = b_diff/se_diff
  // P>|t|
  scalar p_diff = 2*(1-normal(abs(t_diff))) // norm and boot are the same 
  // Confidence interval
  scalar ci_lb_diff = b_diff + invttail(`=r(N)-df', 0.975)*se_diff // fix: not exactly the same as nlcom
  scalar ci_ub_diff = b_diff + invttail(`=r(N)-df', 0.025)*se_diff // fix: not exactly the same as nlcom


  * Display estimation results
  *-------------------------------------------------------------------------------
  // Normal based
  if "`normal'" != "" {
    di as text "{hline 13}{c TT}{hline 64}"
    di as text %12s abbrev("`depvar'",12) " {c |}" ///
      _col(15) "{ralign 11:Coef.}" ///
      _col(26) "{ralign 12:Std. Err.}" ///
      _col(38) "{ralign 8:t }" /// notice extra space
      _col(46) "{ralign 8:P>|t|}" ///
      _col(54) "{ralign 25:[95% Conf. Interval]}" 
    di as text "{hline 13}{c +}{hline 64}"
    di as text "Subgroup" _col(14) "{c |}"
    forvalues g = 0/1 {
      display as text %12s abbrev("`g'",12) " {c |}" ///
        as result ///
        "  " %9.0g b_g`g' ///
        "  " %9.0g se_g`g' ///
        "    " %5.2f t_g`g' ///
        "   " %5.3f p_g`g' ///
        "    " %9.0g ci_lb_g`g' ///
        "   " %9.0g ci_ub_g`g'
    }
    di as text "{hline 13}{c +}{hline 64}"
    display as text "Difference   {c |}" ///
      as result ///
      "  " %9.0g b_diff ///
      "  " %9.0g se_diff ///
      "    " %5.2f t_diff ///
      "   " %5.3f p_diff ///
      "    " %9.0g ci_lb_diff ///
      "   " %9.0g ci_ub_diff
    di as text "{hline 13}{c BT}{hline 64}"
  }
  // Empirical 
  else {
    di as text "{hline 13}{c TT}{hline 64}"
    di as text %12s abbrev("`depvar'",12) " {c |}" ///
      _col(15) "{ralign 11:Coef.}" ///
      _col(26) "{ralign 12:Std. Err.}" ///
      _col(38) "{ralign 8:z }" /// notice extra space
      _col(46) "{ralign 8:P>|z|}" ///
      _col(58) "{ralign 25:[95% Conf. Interval] (P)}" 
    di as text "{hline 13}{c +}{hline 64}"
    di as text "Subgroup" _col(14) "{c |}"
    forvalues g = 0/1 {
      display as text %12s abbrev("`g'",12) " {c |}" ///
        as result ///
        "  " %9.0g b_g`g' ///
        "  " %9.0g se_g`g' ///
        "    " %5.2f t_g`g' ///
        "   " %5.3f p_g`g' ///
        "    " %9.0g ci_lb_g`g' ///
        "   " %9.0g ci_ub_g`g'
    }
    di as text "{hline 13}{c +}{hline 64}"
    display as text "Difference   {c |}" ///
      as result ///
      "  " %9.0g b_diff ///
      "  " %9.0g se_diff ///
      "    " %5.2f t_diff ///
      "   " %5.3f p_diff ///
      "    " %9.0g ci_lb_diff ///
      "   " %9.0g ci_ub_diff
    di as text "{hline 13}{c BT}{hline 64}"
  }
}

* End
*-------------------------------------------------------------------------------
cap drop _cutoff
end

*===============================================================================
* Define auxiliary subroutines
*===============================================================================

*-------------------------------------------------------------------------------
* epost: post matrices in e(b) and e(V); leave other ereturn results unchanged
*-------------------------------------------------------------------------------
program epost, eclass
  // Store results: scalars
  local scalars: e(scalars)
  foreach scalar of local scalars {
    local `scalar' = e(`scalar')
  }
  // Store results: macros
  local macros: e(macros)
  foreach macro of local macros {
    local `macro' = e(`macro')
  }
  // Store results: matrices (drop V_modelbased; b and V are computed below)
  local matrices: e(matrices)
  // Store results: functions
  tempvar esample
  gen `esample' = e(sample)
  // b and V matrices
  matrix b = e(b)
  matrix V = e(V)
  matrix b = b[1, "0.`1'#1.`2'".."1.`1'#1.`2'"]
  matrix V = V["0.`1'#1.`2'".."1.`1'#1.`2'", "0.`1'#1.`2'".."1.`1'#1.`2'"]
  ereturn post, esample(`esample')
  // Post results: scalars
  foreach scalar of local scalars {
    ereturn scalar `scalar' = ``scalar''
  }
  // Post results: macros
  foreach macro of local macros {
    ereturn local `macro' ``macro''
  }
end

*-------------------------------------------------------------------------------
* myboo: compute bootstrapped variance-covariance matrix & adjust ereturn results
*-------------------------------------------------------------------------------
program define myboo, eclass

  // Store results: scalars
  local scalars: e(scalars)

  
  foreach scalar of local scalars {
    local `scalar' = e(`scalar')
  }
  // Store results: macros
  local macros: e(macros)
  foreach macro of local macros {
    local `macro' = e(`macro')
  }
  

  // Store results: matrices (drop V_modelbased; b and V are computed below)
  local matrices: e(matrices)
  // Store results: functions
  tempvar esample
  gen `esample' = e(sample)
  * Extract b submatrix with subgroup coefficients
  matrix b = e(b)
  matrix b = b[1, "0.`1'#1.`2'".."1.`1'#1.`2'"]
  matrix colnames b = 0.`1'#1.`2' 1.`1'#1.`2'
 
  // Start bootstrap 
  di "" // empty line on purpose
  _dots 0, title(Bootstrap replications) reps(`3')
  cap mat drop cumulative // more elegant solution?
  forvalues i=1/`3' {
    preserve
    bsample, strata(`e(blockbootstrap)') // sample w/ replacement; default sample size is _N
    qui `e(cmdline)' // use full regression specification left out by reg
    tempname this_run
    ** If RDD is not fuzzy or we do bootstrap both first stage and reduced form 
    if IndIV[1,1]==0 | fixed[1,1]==0  {
      mat `this_run' = (_b[0.`1'#1.`2'], _b[1.`1'#1.`2'])
    }
    	 
		
    ** If RDD is fuzzy and we do only bootstrap in the reduced form keeping the first stage fixed.
if IndIV[1,1]==1 & fixed[1,1]==1 {  
    local CoeffIVg0_`i'= _b[0.`1'#1._cutoff]/FS[1,1]
    local CoeffIVg1_`i'= _b[1.`1'#1._cutoff]/FS[1,2]
    mat `this_run' = (`CoeffIVg0_`i'',`CoeffIVg1_`i'')
    }  
    
    mat cumulative = nullmat(cumulative) \ `this_run'
    restore
    _dots `i' 0
  }
  
 
  
  di _newline
  // Compute variance-covariance matrix 
  /* This procedure was achieved with the variance mata function, but could be 
  computed with cross() or crossdev() mata functions. */
  cap mat drop V
  mata: cumulative = st_matrix("cumulative")
  mata: st_matrix("V", variance(cumulative)) // see help mf_mean
  /*
  // New computation
  mata: cumulative = st_matrix("cumulative")
  mata: st_matrix("means", mean(cumulative))
  matrix means = means'
  matrix U = J(1,`3',1)
  matrix means = means * U 
  matrix cumulative = cumulative'
  matrix V = cumulative - means
  matrix V = V*V'
  matrix V = V/`=`3'-1'
  */
  // Add names to rows and columns of V 
  mat rownames V = 0.`1'#1.`2' 1.`1'#1.`2'
  mat colnames V = 0.`1'#1.`2' 1.`1'#1.`2'
  // Return 
  ereturn post, esample(`esample')
  // Post results: scalars
  foreach scalar of local scalars {
    ereturn scalar `scalar' = ``scalar''
  }
  ereturn scalar N_reps = `3'
  ereturn scalar level = 95
  // Empirical p-values 
  cap scalar drop bscoef
  forvalues g = 0/1 {
    local count = 0 // reset count for second subgroup 
    forvalues i = 1/`3' {
      scalar bscoef = cumulative[`i',`=`g'+1']
      if abs(bscoef) >= abs(b[1,`=`g'+1']) local count = `count'+1
    }
    scalar pval`g' = (1+`count') / (`3' + 1)
*    di "pval`g' = (1+`count') / (`3' + 1)"
    ereturn scalar pval`g' = pval`g'
  }
  // Empirical confidence intervals
  svmat cumulative, names(_subgroup)
  forvalues g = 0/1 {
    qui centile _subgroup`=`g'+1', centile(2.5 97.5)
    drop _subgroup`=`g'+1'
    scalar lb_g`g' = r(c_1)
    ereturn scalar lb_g`g' = lb_g`g'
    scalar ub_g`g' = r(c_2)
    ereturn scalar ub_g`g' = ub_g`g'
  }
  
  // Post results: macros
  foreach macro of local macros {
    if "`macro'" == "clustvar" continue // skip this macro as it doesn't apply
    ereturn local `macro' ``macro''
  }
  ereturn local vcetype "Bootstrap"
  ereturn local vce "bootstrap"
  ereturn local prefix "bootstrap"
*  ereturn matrix cumulative = cumulative
  // Drop auxiliary matrices 
*  cap mat drop cumulative
*  cap mat drop  means
*  cap mat drop U 
end

*-------------------------------------------------------------------------------
* nlcompost: modify b matrix after nlcom
*-------------------------------------------------------------------------------
program nlcompost, eclass
  matrix b = e(b)
  matrix colnames b = Difference 
  ereturn repost b = b, rename // renames V matrix as well
end

*-------------------------------------------------------------------------------
* nlcomhack: hack b and V matrices to inlude nlcom results
*-------------------------------------------------------------------------------
program nlcomhack, eclass
  tempname b V nlcom_V
  matrix `b' = e(b)
  matrix `V' = e(V)
  local i = colnumb(`b', "_nl_1")
  qui nlcom _b[1.`1'#1.`2'] - _b[0.`1'#1.`2']
  matrix `nlcom_V' = r(V) // for some reason this is necessary
  matrix `b'[1,`i'] = r(b)
  matrix `V'[`i',`i'] = `nlcom_V'[1,1] // ...and this
  ereturn repost b = `b' V = `V'
end


*-------------------------------------------------------------------------------
* balancematrix: compute balance table matrices and other statistics
*-------------------------------------------------------------------------------
program define balancematrix, eclass
syntax, matname(string) /// important inputs, differ by call
  touse(name) weights(string) bwidth(string) balance(varlist) /// unchanging inputs
  [psw ipsweight(name) pscore(name) comsup nocomsup(name) binarymodel(string)] /// only needed for PSW balance
  sgroup(name) sgroup0(name) n_balance(int) // todo: eliminate these? can be computed by subroutine at low cost


* Create variables specific to PSW matrix
*-------------------------------------------------------------------------------
if "`psw'" != "" { // if psw
  // Fit binary response model
 qui cap drop comsup
 qui `binarymodel' `sgroup' `balance' [pw=`weights'] if `touse' & `bwidth'


  // Generate pscore variable and clear stored results
  qui predict double `pscore' if `touse' & `bwidth' & !mi(`sgroup')
  ereturn clear

  // No compute common support area as default (create a aux variable)
if "`nocomsup'" != "" {
  tempvar COMSUP
  qui gen `COMSUP' = 1 if `touse' & `bwidth' & !mi(`sgroup')
}  
else {
qui sum `pscore' if `sgroup' == 1 /* todo: check why this is like that */
tempvar COMSUP
    qui gen `COMSUP' = ///
      (`pscore' >= r(min) & ///
       `pscore' <= r(max))
     
   label var `COMSUP' "Dummy for obs. in common support"
	
   qui g comsup = `COMSUP'

   label var comsup "Dummy for obs. in common support"
}
  // Count observations in each fuzzy group
  qui count if `touse' & `bwidth' & `COMSUP' & `sgroup'==0 & !mi(`pscore')
  local N_G0 = `r(N)'
  qui count if `touse' & `bwidth' & `COMSUP' & `sgroup'==1 & !mi(`pscore')
  local N_G1 = `r(N)'


  // Compute propensity score weighting vector
  cap drop `ipsweight'
  qui gen `ipsweight' = ///
   ( `N_G1'/(`N_G1'+`N_G0')/`pscore'*(`sgroup'==1) + ///
    `N_G0'/(`N_G1'+`N_G0')/(1-`pscore')*(`sgroup'==0)) ///
    if `touse' & `bwidth' & `COMSUP' & !mi(`sgroup')

	tempvar nweights
	qui gen `nweights'=`ipsweight'*`weights'	
    
} // end if psw


* Count obs. in each fuzzy group if not PSW matrix
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
  if "`psw'" == "" qui reg `var' `sgroup0' `sgroup' [iw=`weights'] if `touse' & `bwidth', noconstant /* */
  else qui reg `var' `sgroup0' `sgroup' [iw=`nweights'] if `touse' & `bwidth' & `COMSUP', noconstant
  local coef`j'_G0 = _b[`sgroup0']
  local coef`j'_G1 = _b[`sgroup']
  
  // Compute and store mean differences and their p-values
  if "`psw'" == "" qui reg `var' `sgroup0' [iw=`weights'] if `touse' & `bwidth'
  else qui reg `var' `sgroup0' [iw=`nweights'] if `touse' & `bwidth' & `COMSUP'
  matrix m = r(table)
  scalar diff`j'=m[1,1] // mean difference
  local pval`j' = m[4,1] // p-value 

  // Standardized mean difference
  if "`psw'" == "" qui summ `var' if `touse' & `bwidth' & !mi(`sgroup')
  else qui summ `var' if `touse' & `bwidth' & `COMSUP' & !mi(`sgroup')
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
if "`psw'" == "" qui reg `sgroup' `balance' [iw=`weights']  if `touse' & `bwidth'
else qui reg `sgroup' `balance' [iw=`nweights'] if `touse' & `bwidth' & `COMSUP' 
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

ereturn matrix `matname' = `matname', copy
ereturn scalar `matname'_avgdiff = `avgdiff'
ereturn scalar `matname'_Fstat = `Fstat'
ereturn scalar `matname'_pvalue = `pval_global'
ereturn scalar `matname'_N_G1 = `N_G1'
ereturn scalar `matname'_N_G0 = `N_G0'





end

********************************************************************************

/* 
CHANGE LOG
1.0
  - Compute block bootstrapped variance-covariance matrix
  - Fix First Stage to compute bootstrap in IV
0.9
  - Compute bootstrapped variance-covariance matrix
  - Make program (and subprograms) e-class
  - Allow issuing no model
0.8
  - Add synthetic dataset for examples
0.7 
  - First alpha version ready for full usage
  - Implement nlcom hack to all models, detect diff coef position automatically
0.6
  - Implement nlcom hack to show difference as additional coefficient in ivreg
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
  - Should we use pweights or iweights? iw don't work with ivregress.

TODOS AND IDEAS:
  - Create subroutine of matlist formatting for display of balancematrix output
  - Implement matrix manipulation in Mata
  - Get rid of sgroup0 hack
  - Allow that groupvar is not necessarily an indicator variable
  - Is it possible to allow for N subgroups?
*/
