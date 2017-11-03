*------------------------------------------------------------------------------*
* Test programs with procurement dataset
* Alvaro Carril
*------------------------------------------------------------------------------*

clear all
set seed 1102

// Set root directory
if "`c(os)'" == "MacOSX" cd "/Users/alvaro/Library/Application Support/Stata/ado/personal/rddsga/"
if "`c(os)'" == "Windows" cd "C:\ado\personal\rddsga"

// Set adopath and discard loaded programs
*adopath + old_ados/balancepscore
discard

// Load dataset
use data/procurement_test3, clear
// Obtain bandwidths
*qui rd sh_licitacion high_direct dis_cutoff2
*local w50 = e(w50)
*local w = e(w)
*local w200 = e(w200)

// Run command
rddsga ///
sh_licitacion dis_cutoff /// outcome a-nd assignvar
p1 p2 p3 sh_licitacionPRE2 sh_directoPRE1 sh_licitacionPRE1 sh_directoPRE2 size_PRE1 size_PRE2 I_PREaudit i.gpaoXuceXr /// covariates
 , ///
  ivreg  ///
  balance(p1 p2 p3 size_PRE1 size_PRE2 audited dTR1-dTR3 Year2 Year1 Zone1-Zone3) ///
  bwidth(10) sgroup(high_direct) cutoff(0) treatment(I_CURaudit) ///
  bsreps(10) // ipsweight(peso) pscore(ps_flexmodel4) comsup(comsup) ///
 