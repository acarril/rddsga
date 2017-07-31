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
gen Z = rnormal()
qui summ Z
replace Z = 200/(r(max)-r(min))*(Z-r(max))+100
// Create treatment indicator from running variable
gen T = (Z > 0)
// Generate subgroup indicator
gen G = round(runiform())
// Covariates
gen X1 = .
replace X1 = rnormal() if G
replace X1 = rnormal(.7,0.8) if !G 
gen X2 = .
replace X2 = rnormal() if G
replace X2= rnormal(.7,0.8) if !G 
*rd X1 Z, bw(10) // se cumple para bw=10, balanceado el covariates en ambos lados del cutoff
rd X2 Z, bw(10) // tambien se cumple 
// Create pscore
logit G X1 X2
predict pscore
// Define outcome variable
local bwidth abs(Z)<=100
gen Y = .
replace Y = 1 -  2.5*T + rnormal() if `bwidth' & pscore>0.6 & !G 
replace Y = 1 +  0.1*T + rnormal() if `bwidth' & pscore<0.6 & !G 
replace Y = 0 +    4*T + rnormal() if `bwidth' & pscore<0.3 &  G 
replace Y = 0 + 0.05*T + rnormal() if `bwidth' & pscore>0.3 &  G 
*rd Y T Z
// Tidy up and save dataset
keep Y Z T G X*
order Y Z T G X*
lab var Y "Outcome"
lab var Z "Running variable"
lab var T "Treatment"
lab var G "Subgroup"
lab var X1 "Covariate 1"
lab var X2 "Covariate 2"
// Test rddsga on data
rddsga Y Z, sgroup(G) reduced bw(10) dibalance balance(X1 X2) psweight(weight) quad
// Save
saveold data/rddsga_synth, replace version(11)
saveold rddsga_synth, replace version(11)
// Run plots do-file
run do/plots
