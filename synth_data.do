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
