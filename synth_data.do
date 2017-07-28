*-------------------------------------------------------------------------------
* Generate synthetic dataset for rddsga
* Alvaro Carril
*-------------------------------------------------------------------------------
clear all
discard
set seed 112
* Generate data
*-------------------------------------------------------------------------------
// Define values
local N 10000
local x0 .2
local x1 .5
// Set observations
set obs `N'
// Create running variable (noramlized to [-100, 100])
gen runvar = rnormal()
qui summ runvar
replace runvar = 200/(r(max)-r(min))*(runvar-r(max))+100
// Create treatment indicator from running variable
gen Z = (runvar > 0)
// Generate subgroup indicator
gen G = round(runiform())
// Covariates
gen X1 = .
replace X1 = rnormal() if G
replace X1 = rnormal(.7,0.8) if !G
gen X2 = .
replace X2 = rnormal(-1, 1.5) if G
replace X2 = rnormal(-1.2, 1.2) if !G
gen Y = .
replace Y = 1 + 10*X1 + X2 + .1*Z + rnormal() if G
replace Y = 0 + 1*X1  + X2 - .1*Z + rnormal() if !G
// Save synthetic dataset
compress
*saveold rddsga_synth, replace version(11)
*saveold data/rddsga_synth, replace version(11)

* Estimation
*-------------------------------------------------------------------------------
// Compute optimal bandwidth
rd Y runvar
// Fit unweighted and PSW model
rddsga Y runvar, balance(X1 X2) sgroup(G) bwidth(10) dibal

// Fit unweighted and PSW model
rddsga Y runvar, balance(X1 X2) sgroup(G) bwidth(10) psweight(PSW)

rd Y runvar
rddsga Y runvar, balance(X1 X2) sgroup(G) bwidth(6.095) reduced dibal
rddsga Y runvar X2, balance(X1 X2) sgroup(G) bwidth(6.095) reduced
