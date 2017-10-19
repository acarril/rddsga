*------------------------------------------------------------------------------*
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

rddsga Y X, balance(W1) sgroup(G) bwidth(10) dibal
rddsga Y X, balance(W1 W2) sgroup(G) bwidth(10) ipsweight(ipsw)
rddsga Y X, balance(W1 W2) sgroup(G) bwidth(10) reduced
rddsga Y X W1 W2, sgroup(G) bwidth(6) ivreg bsreps(100) treatment(D) noipsw
rddsga Y X W1 W2, sgroup(G) bwidth(6) ivreg bsreps(100) treatment(D)
