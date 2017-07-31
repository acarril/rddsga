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
replace X2 = rnormal() if G
replace X2= rnormal(.7,0.8) if !G 
*rd X1 runvar, bw(10) // se cumple para bw=10, balanceado el covariates en ambos lados del cutoff
*rd X2 runvar, bw(10) // tambien se cumple 
// Create pscore
logit G X1 X2
predict pscore
// Define outcome variable
local bwidth abs(runvar)<10
gen Y = .
replace Y = 1 +    1*Z + rnormal() if  G & (pscore<0.3 & G) & `bwidth'
replace Y = 1 + 0.05*Z + rnormal() if  G & (pscore>0.3 & G) & `bwidth'
replace Y = 1 -      Z + rnormal() if !G &       pscore>0.6 & `bwidth'
replace Y = 1          + rnormal() if  (pscore<0.4 & !G)    & `bwidth'
// Tidy up and save dataset
keep Y runvar Z G X*
order Y runvar Z G X*
lab var Y "Outcome"
lab var runvar "Running variable"
lab var Z "Treatment"
lab var G "Subgroup"
lab var X1 "Covariate 1"
lab var X2 "Covariate 2"
saveold data/rddsga_synth, replace version(11)
saveold rddsga_synth, replace version(11)
// Test rddsga on data
rddsga Y runvar, sgroup(G) reduced bw(10) dibalance balance(X1 X2) psweight(weight) quad
