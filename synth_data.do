* Generate synthetic dataset

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
gen X = rnormal()
gen Y = .
replace Y = 1 + .6*X + 2*Z + rnormal() if G
replace Y = 0 + .4*X - 2*Z + rnormal() if !G

// Estimation
reg Y X T##G
rdrobust Y runvar, covs(G)