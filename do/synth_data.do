*-------------------------------------------------------------------------------
* Generate synthetic dataset for rddsga
* Alvaro Carril
*-------------------------------------------------------------------------------
// Set root directory
if "`c(os)'" == "MacOSX" cd "/Users/alvaro/Library/Application Support/Stata/ado/personal/rddsga/"
if "`c(os)'" == "Windows" cd "C:\ado\personal\rddsga"
// Initial set up 
clear all
discard
set seed 112
// Define values
local N 10000
local x0 .2
local x1 .5
// Set observations
set obs `N'
// Create running variable (noramlized to [-100, 100])
gen X = rnormal()
qui summ X
replace X = 200/(r(max)-r(min))*(X-r(max))+100
// Create treatment indicator from running variable
gen D = (X + runiform(-1,1) > 0)
// Generate subgroup indicator
gen G = round(runiform())
// Covariates
gen W1 = .
replace W1 = rnormal() if G
replace W1 = rnormal(.7,0.8) if !G 
gen W2 = .
replace W2 = rnormal() if G
replace W2= rnormal(.7,0.8) if !G 
*rd W1 X, bw(10) // se cumple para bw=10, balanceado el covariates en ambos lados del cutoff
rd W2 X, bw(10) // tambien se cumple 
// Create pscore
logit G W1 W2
predict pscore
// Define outcome variable
local bwidth abs(X)<=100
gen Y = .
replace Y = 0.5 -  4*D+0.5*(W1-W2) + rnormal() if `bwidth' & pscore>0.7 & !G 
replace Y = 0.4*D+ rnormal() +0.5*(W1-W2) if `bwidth' & pscore<=0.7 & !G 
replace Y = 1 +    3.5*D  + rnormal() + 0.5*(W1-W2) if `bwidth' & pscore<0.4 &  G 
replace Y = 0.5+-0.7*D+ rnormal() + 0.5*(W1-W2) if `bwidth' & pscore>=0.4 &  G 
*rd Y D X
// Tidy up and save dataset
keep Y X D G W*
order Y X D G W*
lab var Y "Outcome"
lab var X "Running variable"
lab var D "Dreatment"
lab var G "Subgroup"
lab var W1 "Covariate 1"
lab var W2 "Covariate 2"
// Test rddsga on data
*rddsga Y X, sgroup(G) reduced bw(10) dibalance balance(W1 W2) psweight(weight) quad

*drop _est_*
// Save
saveold rddsga_synth, replace
// Run plots do-file
*run do/plots
