* Generate synthetic dataset
*** Andre
clear all
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
gen X = .
replace X = rnormal() if G
replace X = rnormal(.7,0.8) if !G 
gen Y = .
replace Y = 1 + .6*X + 2*Z + rnormal() if G
replace Y = 0 + .4*X - 2*Z + rnormal() if !G

// Estimation
*reg Y X Z##G
*rdrobust Y runvar, covs(G)
rddsga Y runvar, balance(X) sgroup(G) bwidth(10) reduced dibal
