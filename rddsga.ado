program define rddsga, rclass byable(recall)
version 11.1
syntax varlist [if] [in] ,  weight(varlist) treat(varlist) bw(string) [c(real 0) x(varlist) covs(varlist) fe(varlist) cluster(string) ivreg FirstStage ReducedForm all  spline(string) addnamtex(string) latex dir(string) bdec(int 3)  labelt0(string)]

tokenize `varlist'


/* retrieve the Outcome */
local Y  `1' 

/* retrieve the distance to the cutoff */
local X  `2' 

*** SPLINE
if "`spline'" == "" { 
*** we define three splines. The default is linear, linear, liner, but it could be included quadratic regression.
local S "L" "L" "L"

forval i=1/3 {
local S`i': word `i' of `S'
}

}

*** Interaction with the treatment
tokenize `treat'
/* retrieve the treatment indicator */
local T  `1'

************* LINEAR

/* indicator */
tempvar g
g `g'=(`X'>=`c')

/* interaction */
tempvar gX
g `gX'=`g'*`X'

** define local containing cutoff, instrument,  interaction between instrument and cutoff and the covariates. 
local VL `x' `g' `X' `gX' `covs'

tempvar T0
qui gen `T0'=(`T'==0) if `T'!=.

* generating interaction terms between treatment and each of the covariates 
foreach var of varlist `VL' {
tempvar T1`var'
cap g `T1`var''=`T'*`var'
tempvar T0`var'
cap g `T0`var''=`T0'*`var'


local VARL `VARL' `T1`var'' `T0`var''
}


*************** QUADRATIC
tempvar Xsq
qui g `Xsq'=`X'^2

tempvar gXsq
qui g `gXsq'=`gX'^2

local VQ `Xsq'  `gXsq' 

foreach var of varlist `VQ' {
tempvar T1`var'
qui cap g `T1`var''=`T'*`var'
tempvar T0`var'
qui cap g `T0`var''=`T0'*`var'

local VARQ1 `VARQ1' `T0`var'' `T1`var''
}
local VARQ `VARL' `VARQ1' 



tempvar touse
g `touse'=0

qui replace `touse'=1 `if' `in'



/* Bandwidth  */
if "`bw'" == "" { 
   qui ssc install rd
   tokenize `varlist'
 rd `1' `2'
 tempname bw
 qui matrix `bw'=(e(w50),e(w), e(w200))
 qui matrix list `bw'
}

******************************* weight ******

 tokenize weight
 local numw `: word count `numw''
 local col=colsof(`bw')
 
 

**** WEIGHT
forval i=1/`col' {
tokenize `weight'
local w`i' ``i''
}

**FIXED EFFECT
foreach var of varlist `fe' {

tempvar T0`var'
cap g `T0`var''=`T0'*`var'

tempvar T`var'
cap g `T`var''=`T'*`var'

local FE `FE' i.`T0`var'' i.`T`var''
	
}

*** FIRST STAGE

if `"`FirstStage'"' != `""' | `"`all'"' != `""' { 

forval i=1/`col' {
local S`i': word `i' of `spline'



tokenize `VAR`S`i'''
local numvar`S`i'' `: word count `VAR`S`i''''


*** independent variable 
** 0
local X0 `1'
** 1 
local X1 `2'

** Instrument
*** 0
local Z0 `3'
*** 1
local Z1 `4'

*** save other variables
forvalue j=5/`numvar`S`i''' {
local C`S`i'' `C`S`i''' ``j''
}

local bw`i'=`bw'[1,`i']


qui xi: reg `x' `Z0' `Z1' `C`S`i''' `FE'  if `X'>-(`bw`i'') & `X'<(`bw`i''), vce(cluster `cluster')

mat m=r(table)
mat VAR=e(V)

	/*Diff mean*/    local C`i'R5:  di %12.`bdec'fc m[1,1]-m[1,2]
	/* s.e. for diff*/ local C`i'R6: di %12.`bdec'fc (sqrt(m[2,1]^2+m[2,2]^2-2*VAR[2,1])) // VAR[2,1] corresponds to the covariance. 
						local C`i'R6: di trim("`C`i'R6'")
	/*(s.e. diff) */    local C`i'R6="(`C`i'R6')"
        /* p-val test */ local pv`i'R5: di %12.`bdec'fc 2*(1-t(e(df_r),abs((m[1,1]-m[1,2])/(sqrt(m[2,1]^2+m[2,2]^2-2*VAR[2,1])))))
        /*coef H  */     local C`i'R1: di %12.`bdec'fc m[1,1]
	/*s.e. H  */     local C`i'R2: di %12.`bdec'fc m[2,1]
				local C`i'R2: di trim("`C`i'R2'")
	/*(s.e. H) */    local C`i'R2="(`C`i'R2')"
	/*p-val H */     local pv`i'R1: di %12.`bdec'fc m[4,1]
	/*coef L  */     local C`i'R3: di %12.`bdec'fc m[1,2]
	/*s.e. L  */     local C`i'R4: di %12.`bdec'fc m[2,2]
				local C`i'R4: di trim("`C`i'R4'")
	/*(s.e. L) */    local C`i'R4="(`C`i'R4')"
	/*p-val L */     local pv`i'R3: di %12.`bdec'fc m[4,2]
	/*bandwidth */   local C`i'R7="$\pm`bw`i''$"		
	/*numObs */      local C`i'R8: di %12.0fc e(N)
	if "`S`i''"=="L" {
	/*spline*/  	 local C`i'R10="Linear"	
				}		
	else {
	/*spline*/  	 local C`i'R10="Quadr."	

			}
if `"`fe'"' != "" {			
	/*fe */      	local C`i'R11="Yes"
	/*control */    local C`i'R12="Yes"
	}
else {			
	/*fe */      	local C`i'R11="No"
	/*control */    local C`i'R12="No"
	}
	
*stars
	forvalues j = 1(2)5 {
	if `pv`i'R`j''<=0.10 & `pv`i'R`j''>0.05 {
	local C`i'R`j'="$`C`i'R`j''^{*}$"
	}
	
	else if `pv`i'R`j''<=0.05 & `pv`i'R`j''>0.01 {
	local C`i'R`j'="$`C`i'R`j''^{**}$"
	}
	else if `pv`i'R`j''<=0.01 { 
	local C`i'R`j'="$`C`i'R`j''^{***}$"
	}	
	}
	
	
qui xi: reg `x' `Z0' `Z1' `C`S`i''' `FE' [iw=`w`i''] if `X'>-(`bw`i'') & `X'<(`bw`i''), vce(cluster `cluster')
mat m=r(table)
mat VAR=e(V)

	/*Diff mean*/    local C`i'R5_w:  di %12.`bdec'fc m[1,1]-m[1,2]
	/* s.e. for diff*/ local C`i'R6_w: di %12.`bdec'fc (sqrt(m[2,1]^2+m[2,2]^2-2*VAR[2,1]))
						local C`i'R6: di trim("`C`i'R6_w'")
	/*(s.e. diff) */    local C`i'R6_w="(`C`i'R6_w')"
        /* p-val test */ local pv`i'R5_w: di %12.`bdec'fc 2*(1-t(e(df_r),abs((m[1,1]-m[1,2])/(sqrt(m[2,1]^2+m[2,2]^2-2*VAR[2,1])))))
        /*coef H  */     local C`i'R1_w: di %12.`bdec'fc m[1,1]
	/*s.e. H  */     local C`i'R2_w: di %12.`bdec'fc m[2,1]
				local C`i'R2_w: di trim("`C`i'R2_w'")
	/*(s.e. H) */    local C`i'R2_w="(`C`i'R2_w')"
	/*p-val H */     local pv`i'R1_w: di %12.`bdec'fc m[4,1]
	/*coef L  */     local C`i'R3_w: di %12.`bdec'fc m[1,2]
	/*s.e. L  */     local C`i'R4_w: di %12.`bdec'fc m[2,2]
				local C`i'R4_w: di trim("`C`i'R4_w'")
	/*(s.e. L) */    local C`i'R4_w="(`C`i'R4_w')"
	/*p-val L */     local pv`i'R3_w: di %12.`bdec'fc m[4,2]
	/*bandwidth */   local C`i'R7_w="$\pm`bw`i''$"		
	/*numObs */      local C`i'R8_w: di %12.0fc e(N)

	*stars
	forval j=1(2)5 {
	if `pv`i'R`j'_w'<=0.10 & `pv`i'R`j'_w'>0.05 {
	local C`i'R`j'_w="$`C`i'R`j'_w'^{*}$"
	}
	else if `pv`i'R`j'_w'<=0.05 & `pv`i'R`j'_w'>0.01 {
	local C`i'R`j'_w="$`C`i'R`j'_w'^{**}$"
	}
	else if `pv`i'R`j'_w'<=0.01 { 
	local C`i'R`j'_w="$`C`i'R`j'_w'^{***}$"
	}
}
}

if "`latex'" != "" { 
********************************************************************************
	
*** Row titles
local Zvar: variable label `treat'

	local lbl_1 "1\{RI $\ge$ cutoff\} $\times$ `Zvar'" 
	local lbl_2 ""
	local lbl_3 "1\{RI $\ge$ cutoff\} $\times$ `labelt0'" 
	local lbl_4 ""
	local lbl_5 "Difference Estimate"
	local lbl_6 ""
	local lbl_7 "Bandwidth"
	local lbl_8 "Observations"
	local lbl_9 "R-squared"
 
	local lbl_10 "Spline"
	local lbl_11 "Stratum fixed effects"
	local lbl_12 "Additional controls"
	
********************************************************************************
* Create TeX file
********************************************************************************	
texdoc init $TEX/FirstStage_`treat'_`addnamtex'.tex, replace
tex \\ \\ [-1.5ex]
tex \hline\hline \\ [-1.5ex]
tex {} & (1) & (2) & (3)   \\
tex [1ex] \\ [-1.5ex]
tex & \multicolumn{3}{c}{\quad Panel A: Nonweighted}\\\\
forvalue t=1/8 {
tex `lbl_`t''  & `C1R`t''  & `C2R`t'' & `C3R`t''  \\
}
tex \\
tex & \multicolumn{3}{c}{\quad Panel B: Propensity score-weighted}\\\\
forvalue t=1/8 {
tex `lbl_`t''  & `C1R`t'_w'  & `C2R`t'_w' & `C3R`t'_w'  \\
}
tex [1ex] \hline \\ [-1.5ex]
forvalue t=10/12 {
tex `lbl_`t''  & `C1R`t''  & `C2R`t'' & `C3R`t''  \\
}
tex [1ex] \hline\hline \\ [-1.5ex]
texdoc close

}
}


*** ReducedForm

if `"`ReducedForm'"' != `""' | `"`all'"' != `""' { 

forval i=1/`col' {
local S`i': word `i' of `spline'



tokenize `VAR`S`i'''
local numvar`S`i'' `: word count `VAR`S`i''''


*** xendent variable 
** 0
local X0 `1'
** 1 
local X1 `2'

** Instrument
*** 0
local Z0 `3'
*** 1
local Z1 `4'

*** save other variables
forvalue j=5/`numvar`S`i''' {
local C`S`i'' `C`S`i''' ``j''
}

local bw`i'=`bw'[1,`i']


qui xi: reg `Y' `Z0' `Z1' `C`S`i''' `FE'  if `X'>-(`bw`i'') & `X'<(`bw`i''), vce(cluster `cluster')

mat m=r(table)
mat VAR=e(V)

	/*Diff mean*/    local C`i'R5:  di %12.`bdec'fc m[1,1]-m[1,2]
	/* s.e. for diff*/ local C`i'R6: di %12.`bdec'fc (sqrt(m[2,1]^2+m[2,2]^2-2*VAR[2,1]))
						local C`i'R6: di trim("`C`i'R6'")
	/*(s.e. diff) */    local C`i'R6="(`C`i'R6')"
        /* p-val test */ local pv`i'R5: di %12.`bdec'fc 2*(1-t(e(df_r),abs((m[1,1]-m[1,2])/(sqrt(m[2,1]^2+m[2,2]^2-2*VAR[2,1])))))
        /*coef H  */     local C`i'R1: di %12.`bdec'fc m[1,1]
	/*s.e. H  */     local C`i'R2: di %12.`bdec'fc m[2,1]
				local C`i'R2: di trim("`C`i'R2'")
	/*(s.e. H) */    local C`i'R2="(`C`i'R2')"
	/*p-val H */     local pv`i'R1: di %12.`bdec'fc m[4,1]
	/*coef L  */     local C`i'R3: di %12.`bdec'fc m[1,2]
	/*s.e. L  */     local C`i'R4: di %12.`bdec'fc m[2,2]
				local C`i'R4: di trim("`C`i'R4'")
	/*(s.e. L) */    local C`i'R4="(`C`i'R4')"
	/*p-val L */     local pv`i'R3: di %12.`bdec'fc m[4,2]
	/*bandwidth */   local C`i'R7="$\pm`bw`i''$"		
	/*numObs */      local C`i'R8: di %12.0fc e(N)
	if "`S`i''"=="L" {
	/*spline*/  	 local C`i'R10="Linear"	
				}		
	else {
	/*spline*/  	 local C`i'R10="Quadr."	

			}
			
if `"`fe'"' != "" {			
	/*fe */      	local C`i'R11="Yes"
	/*control */    local C`i'R12="Yes"
	}
else {			
	/*fe */      	local C`i'R11="No"
	/*control */    local C`i'R12="No"
	}
*stars
	forvalues j = 1(2)5 {
	if `pv`i'R`j''<=0.10 & `pv`i'R`j''>0.05 {
	local C`i'R`j'="$`C`i'R`j''^{*}$"
	}
	
	else if `pv`i'R`j''<=0.05 & `pv`i'R`j''>0.01 {
	local C`i'R`j'="$`C`i'R`j''^{**}$"
	}
	else if `pv`i'R`j''<=0.01 { 
	local C`i'R`j'="$`C`i'R`j''^{***}$"
	}	
	}
		
	

qui xi: reg `Y' `Z0' `Z1' `C`S`i''' `FE' [iw=`w`i''] if `X'>-(`bw`i'') & `X'<(`bw`i''), vce(cluster `cluster')
mat m=r(table)
mat VAR=e(V)

	/*Diff mean*/    local C`i'R5_w:  di %12.`bdec'fc m[1,1]-m[1,2]
	/* s.e. for diff*/ local C`i'R6_w: di %12.`bdec'fc (sqrt(m[2,1]^2+m[2,2]^2-2*VAR[2,1]))
						local C`i'R6: di trim("`C`i'R6_w'")
	/*(s.e. diff) */    local C`i'R6_w="(`C`i'R6_w')"
        /* p-val test */ local pv`i'R5_w: di %12.`bdec'fc 2*(1-t(e(df_r),abs((m[1,1]-m[1,2])/(sqrt(m[2,1]^2+m[2,2]^2-2*VAR[2,1])))))
        /*coef H  */     local C`i'R1_w: di %12.`bdec'fc m[1,1]
	/*s.e. H  */     local C`i'R2_w: di %12.`bdec'fc m[2,1]
				local C`i'R2_w: di trim("`C`i'R2_w'")
	/*(s.e. H) */    local C`i'R2_w="(`C`i'R2_w')"
	/*p-val H */     local pv`i'R1_w: di %12.`bdec'fc m[4,1]
	/*coef L  */     local C`i'R3_w: di %12.`bdec'fc m[1,2]
	/*s.e. L  */     local C`i'R4_w: di %12.`bdec'fc m[2,2]
				local C`i'R4_w: di trim("`C`i'R4_w'")
	/*(s.e. L) */    local C`i'R4_w="(`C`i'R4_w')"
	/*p-val L */     local pv`i'R3_w: di %12.`bdec'fc m[4,2]
	/*bandwidth */   local C`i'R7_w="$\pm`bw`i''$"		
	/*numObs */      local C`i'R8_w: di %12.0fc e(N)	
	
*stars
	forval j=1(2)5 {
	if `pv`i'R`j'_w'<=0.10 & `pv`i'R`j'_w'>0.05 {
	local C`i'R`j'_w="$`C`i'R`j'_w'^{*}$"
	}
	else if `pv`i'R`j'_w'<=0.05 & `pv`i'R`j'_w'>0.01 {
	local C`i'R`j'_w="$`C`i'R`j'_w'^{**}$"
	}
	else if `pv`i'R`j'_w'<=0.01 { 
	local C`i'R`j'_w="$`C`i'R`j'_w'^{***}$"
	}
}
	
}

if "`latex'" != "" { 
********************************************************************************
	
*** Row titles
local Zvar: variable label `treat'

	local lbl_1 "1\{RI $\ge$ cutoff\} $\times$ `Zvar'" 
	local lbl_2 ""
	local lbl_3 "1\{RI $\ge$ cutoff\} $\times$ `labelt0'" 
	local lbl_4 ""
	local lbl_5 "Difference Estimate"
	local lbl_6 ""
	local lbl_7 "Bandwidth"
	local lbl_8 "Observations"
	local lbl_9 "R-squared"
 
	local lbl_10 "Spline"
	local lbl_11 "Stratum fixed effects"
	local lbl_12 "Additional controls"
	
********************************************************************************
* Create TeX file
********************************************************************************	
texdoc init $TEX/ReducedForm_`treat'_`addnamtex'.tex, replace
tex \\ \\ [-1.5ex]
tex \hline\hline \\ [-1.5ex]
tex {} & (1) & (2) & (3)   \\
tex [1ex] \\ [-1.5ex]
tex & \multicolumn{3}{c}{\quad Panel A: Nonweighted}\\\\
forvalue t=1/8 {
tex `lbl_`t''  & `C1R`t''  & `C2R`t'' & `C3R`t''  \\
}
tex \\
tex & \multicolumn{3}{c}{\quad Panel B: Propensity score-weighted}\\\\
forvalue t=1/8 {
tex `lbl_`t''  & `C1R`t'_w'  & `C2R`t'_w' & `C3R`t'_w'  \\
}
tex [1ex] \hline \\ [-1.5ex]
forvalue t=10/12 {
tex `lbl_`t''  & `C1R`t''  & `C2R`t'' & `C3R`t''  \\
}
tex [1ex] \hline\hline \\ [-1.5ex]
texdoc close

}
}




*** IVREG

if `"`ivreg'"' != `""' | `"`all'"' != `""' { 

forval i=1/`col' {
local S`i': word `i' of `spline'



tokenize `VAR`S`i'''
local numvar`S`i'' `: word count `VAR`S`i''''


*** xendent variable 
** 0
local X0 `1'
** 1 
local X1 `2'

** Instrument
*** 0
local Z0 `3'
*** 1
local Z1 `4'

*** save other variables
forvalue j=5/`numvar`S`i''' {
local C`S`i'' `C`S`i''' ``j''
}

local bw`i'=`bw'[1,`i']



 qui xi: ivreg `Y' `C`S`i''' `FE' (`X0' `X1' = `Z0' `Z1') if `X'>-(`bw`i'') & `X'<(`bw`i''), cluster(`cluster')

mat m=r(table)
mat VAR=e(V)

	/*Diff mean*/    local C`i'R5:  di %12.`bdec'fc m[1,1]-m[1,2]
	/* s.e. for diff*/ local C`i'R6: di %12.`bdec'fc (sqrt(m[2,1]^2+m[2,2]^2-2*VAR[2,1]))
						local C`i'R6: di trim("`C`i'R6'")
	/*(s.e. diff) */    local C`i'R6="(`C`i'R6')"
        /* p-val test */ local pv`i'R5: di %12.`bdec'fc 2*(1-t(e(df_r),abs((m[1,1]-m[1,2])/(sqrt(m[2,1]^2+m[2,2]^2-2*VAR[2,1])))))
        /*coef H  */     local C`i'R1: di %12.`bdec'fc m[1,1]
	/*s.e. H  */     local C`i'R2: di %12.`bdec'fc m[2,1]
				local C`i'R2: di trim("`C`i'R2'")
	/*(s.e. H) */    local C`i'R2="(`C`i'R2')"
	/*p-val H */     local pv`i'R1: di %12.`bdec'fc m[4,1]
	/*coef L  */     local C`i'R3: di %12.`bdec'fc m[1,2]
	/*s.e. L  */     local C`i'R4: di %12.`bdec'fc m[2,2]
				local C`i'R4: di trim("`C`i'R4'")
	/*(s.e. L) */    local C`i'R4="(`C`i'R4')"
	/*p-val L */     local pv`i'R3: di %12.`bdec'fc m[4,2]
	/*bandwidth */   local C`i'R7="$\pm`bw`i''$"		
	/*numObs */      local C`i'R8: di %12.0fc e(N)
	if "`S`i''"=="L" {
	/*spline*/  	 local C`i'R10="Linear"	
				}		
	else {
	/*spline*/  	 local C`i'R10="Quadr."	

			}
if `"`fe'"' != "" {			
	/*fe */      	local C`i'R11="Yes"
	/*control */    local C`i'R12="Yes"
	}
else {			
	/*fe */      	local C`i'R11="No"
	/*control */    local C`i'R12="No"
	}
	
*stars
	forvalues j = 1(2)5 {
	if `pv`i'R`j''<=0.10 & `pv`i'R`j''>0.05 {
	local C`i'R`j'="$`C`i'R`j''^{*}$"
	}
	
	else if `pv`i'R`j''<=0.05 & `pv`i'R`j''>0.01 {
	local C`i'R`j'="$`C`i'R`j''^{**}$"
	}
	else if `pv`i'R`j''<=0.01 { 
	local C`i'R`j'="$`C`i'R`j''^{***}$"
	}	
	}
	
	

 qui xi: ivreg `Y' `C`S`i''' `FE' (`X0' `X1' = `Z0' `Z1') [iw=`w`i''] if `X'>-(`bw`i'') & `X'<(`bw`i''), cluster(`cluster')
mat m=r(table)
mat VAR=e(V)

	/*Diff mean*/    local C`i'R5_w:  di %12.`bdec'fc m[1,1]-m[1,2]
	/* s.e. for diff*/ local C`i'R6_w: di %12.`bdec'fc (sqrt(m[2,1]^2+m[2,2]^2-2*VAR[2,1]))
						local C`i'R6: di trim("`C`i'R6_w'")
	/*(s.e. diff) */    local C`i'R6_w="(`C`i'R6_w')"
        /* p-val test */ local pv`i'R5_w: di %12.`bdec'fc 2*(1-t(e(df_r),abs((m[1,1]-m[1,2])/(sqrt(m[2,1]^2+m[2,2]^2-2*VAR[2,1])))))
        /*coef H  */     local C`i'R1_w: di %12.`bdec'fc m[1,1]
	/*s.e. H  */     local C`i'R2_w: di %12.`bdec'fc m[2,1]
				local C`i'R2_w: di trim("`C`i'R2_w'")
	/*(s.e. H) */    local C`i'R2_w="(`C`i'R2_w')"
	/*p-val H */     local pv`i'R1_w: di %12.`bdec'fc m[4,1]
	/*coef L  */     local C`i'R3_w: di %12.`bdec'fc m[1,2]
	/*s.e. L  */     local C`i'R4_w: di %12.`bdec'fc m[2,2]
				local C`i'R4_w: di trim("`C`i'R4_w'")
	/*(s.e. L) */    local C`i'R4_w="(`C`i'R4_w')"
	/*p-val L */     local pv`i'R3_w: di %12.`bdec'fc m[4,2]
	/*bandwidth */   local C`i'R7_w="$\pm`bw`i''$"		
	/*numObs */      local C`i'R8_w: di %12.0fc e(N)	
	
*stars
	forval j=1(2)5 {
	if `pv`i'R`j'_w'<=0.10 & `pv`i'R`j'_w'>0.05 {
	local C`i'R`j'_w="$`C`i'R`j'_w'^{*}$"
	}
	else if `pv`i'R`j'_w'<=0.05 & `pv`i'R`j'_w'>0.01 {
	local C`i'R`j'_w="$`C`i'R`j'_w'^{**}$"
	}
	else if `pv`i'R`j'_w'<=0.01 { 
	local C`i'R`j'_w="$`C`i'R`j'_w'^{***}$"
	}
}	
	
}

if "`latex'" != "" { 
********************************************************************************
	
*** Row titles

local Xvar: variable label `x'
local Zvar: variable label `treat'

	local lbl_1 "`Xvar' $\times$ `Zvar'" 
	local lbl_2 ""
	local lbl_3 "`Xvar' $\times$ `labelt0'" 
	local lbl_4 ""
	local lbl_5 "Difference Estimate"
	local lbl_6 ""
	local lbl_7 "Bandwidth"
	local lbl_8 "Observations"
	local lbl_9 "R-squared"

	local lbl_10 "Spline"
	local lbl_11 "Stratum fixed effects"
	local lbl_12 "Additional controls"
	
********************************************************************************
* Create TeX file
********************************************************************************	
texdoc init $TEX/SubgroupRDD_`treat'_`addnamtex'.tex, replace
tex \\ \\ [-1.5ex]
tex \hline\hline \\ [-1.5ex]
tex {} & (1) & (2) & (3)   \\
tex [1ex] \\ [-1.5ex]
tex & \multicolumn{3}{c}{\quad Panel A: Nonweighted}\\\\
forvalue t=1/8 {
tex `lbl_`t''  & `C1R`t''  & `C2R`t'' & `C3R`t''  \\
}
tex \\
tex & \multicolumn{3}{c}{\quad Panel B: Propensity score-weighted}\\\\
forvalue t=1/8 {
tex `lbl_`t''  & `C1R`t'_w'  & `C2R`t'_w' & `C3R`t'_w'  \\
}
tex [1ex] \hline \\ [-1.5ex]
forvalue t=10/12 {
tex `lbl_`t''  & `C1R`t''  & `C2R`t'' & `C3R`t''  \\
}
tex [1ex] \hline\hline \\ [-1.5ex]
texdoc close

}
}



end
