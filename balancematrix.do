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

capture program drop balancematrix
program define balancematrix, rclass
/* 1:matrix_out 2:numcov 3:matrix_in */
tempname `1'
matrix `1' = J(`numcov'+4,4,.)
matrix colnames `1' = "Mean `G0'" "Mean `G1'" "StMeanDiff" p-value 
matrix rownames `1' = `rown3' Observations Abs(StMeanDiff) F-statistic p-value
                         
forvalues j = 1/`numcov' {
  local j=`j'+1  
  matrix `1'[`j',1] = round(`m`j'1',10^(-`bdec')) 
  matrix `1'[`j',2] = round(`m`j'2',10^(-`bdec'))
  matrix `1'[`j',3] = round(`m`j'3',10^(-`bdec'))
  matrix `1'[`j',4] = round(`m`j'4',10^(-`bdec'))
  local rown3 "`rown3' `var'"
}

matrix `1'[`numcov'+1,1] = `Ncontrols'
matrix `1'[`numcov'+1,2] = `Ntreated'     
local l=`numcov'+1
matrix `1'[`numcov'+2,3] = round(`m`l'3',10^(-`bdec'))
local l=`l'+1   
matrix `1'[`numcov'+3,4] = round(`m`l'4',10^(-`bdec'))
local l=`l'+1     
matrix `1'[`numcov'+4,4] = round(`m`l'4',10^(-`bdec'))

end
