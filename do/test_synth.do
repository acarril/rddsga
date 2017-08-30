*---------------------------------------------------------------*
* Test programs with synthetic dataset
* Alvaro Carril
*------------------------------------------------------------------------------*

clear all
set seed 1102
set matsize 11000

// Set root directory
if "`c(os)'" == "MacOSX" cd "/Users/alvaro/Library/Application Support/Stata/ado/personal/rddsga/"
if "`c(os)'" == "Windows" cd "C:\ado\personal\rddsga"

// Set adopath and discard loaded programs
*adopath + old_ados/balancepscore
discard

// Load dataset
use rddsga_synth, clear

// Examples
rddsga Y Z, balance(X1) sgroup(G) bwidth(10) dibal
rddsga Y Z, balance(X1 X2) sgroup(G) bwidth(10) ipsweight(ipsw)
rddsga Y Z, balance(X1 X2) sgroup(G) bwidth(10) reduced
rddsga Y Z X1 X2, sgroup(G) bwidth(6) ivreg bsreps(10) treatment(T) noipsw
rddsga Y Z X1 X2, sgroup(G) bwidth(6) ivreg bsreps(10) treatment(T)
