*------------------------------------------------------------------------------*
* Test programs with procurement dataset
* Alvaro Carril
*------------------------------------------------------------------------------*

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
  psweight(peso) pscore(ps_flexmodel41) comsup(soporte) bdec(5) logit
