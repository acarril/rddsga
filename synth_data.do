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
replace X2 = rnormal() if G
replace X2= rnormal(.7,0.8) if !G 

rd X1 runvar, bw(10) // se cumple para bw=10, balanceado el covariates en ambos lados del cutoff
rd X2 runvar, bw(10) // tambien se cumple 

logit G X1 X2
predict pscore 

gen Y = .
replace Y = 1  + 1*Z + rnormal() if G &  (pscore<0.3 & G)   & abs(runvar)<10
replace Y = 1  + 0.05*Z + rnormal() if G & (pscore>0.3 & G) & abs(runvar)<10
replace Y = 1  - Z + rnormal() if    !G  & pscore>0.6 & abs(runvar)<10
replace Y = 1  + rnormal() if  (pscore<0.4 & !G) & abs(runvar)<10

rddsga Y runvar, sgroup(G) reduced bw(10) dibalance balance(X1 X2) psweight(weight) quad

gen bin=floor(runvar/2)*2+1


keep if abs(runvar)<=10

*** Bin mean (all)

rd Y runvar, bw(10)



		bysort bin: egen bin_Y = mean(Y) 

		graph twoway  (scatter bin_Y  bin,  msymbol(O) mcolor(gray) msize(medium)) ///
		(qfit Y runvar if  runvar>=-10 & runvar<=0, range(-10 0) lcolor(green) lpattern(solid)) /// 
		(qfit Y runvar if  runvar>=0 & runvar<=10, range(0 10) lcolor(green) lpattern(solid)) , ///
		 graphregion(fcol(white)) bgcolor(white) ylabel(0.5(0.5)1.5) xlabel(-10(2)10) xline(0) plotregion(style(none)) /// title("`t', `Atit'", margin(medium) size(medium)) ///
		 legend(off) ytitle("Outcome", size(large))  xtitle("Distance to cutoff", size(large))  ///
		saving(rdplot_all, replace)	

		stop

	* Bin means by group:

		bysort bin G: egen bin_mean_Y = mean(Y) 

		*** MEDIA PONDERADA
	
	* Bin means:

 bysort bin G: egen bin_mean_Yw = wtmean(Y), weight(weight)	
		
	* Draw plot:
	sort runvar

		graph twoway  		///(scatter bin_mean_Y  bin if  G,  msymbol(O) mcolor(gray) msize(medium)) ///
		(qfit Y runvar if  G & runvar>=-10 & runvar<=0, range(-10 0) lcolor(green) lpattern(solid)) /// 
		(qfit Y runvar if G & runvar>=0 & runvar<=10, range(0 10) lcolor(green) lpattern(solid))  ///
		///(scatter bin_mean_Y  bin if  G,  msymbol(O) mcolor(gray) msize(medium)) ///
		(qfit Y runvar if  G & runvar>=-10 & runvar<=0, range(-10 0) lcolor(green) lpattern(solid)) /// 
		(qfit Y runvar if G & runvar>=0 & runvar<=10, range(0 10) lcolor(green) lpattern(solid))  ///
		///(scatter bin_mean_Y  bin if  G==0,  msymbol(O) mcolor(gray) msize(medium)) ///
		(qfit Y runvar if  G==0 & runvar>=-10 & runvar<=0, range(-10 0) lcolor(blue) lpattern(solid)) /// 
		(qfit Y runvar if G==0 & runvar>=0 & runvar<=10, range(0 10) lcolor(blue) lpattern(solid)) ///
		///(scatter bin_mean_Yw  bin if G, msymbol(O) mcolor(gray) msize(medium)) ///
		(qfit Y  runvar [pweight = weight] if G & runvar>=-10 & runvar<=0, range(-10 0) lcolor(black) lpattern(solid)) /// 
		(qfit Y  runvar [pweight = weight] if G & runvar>=0 & runvar<=10, range(0 10) lcolor(black) lpattern(solid)) ///
		///(scatter bin_mean_Yw  bin if G==0, msymbol(O) mcolor(gray) msize(medium)) ///
		(qfit Y  runvar [pweight = weight] if G==0 & runvar>=-10 & runvar<=0, range(-10 0) lcolor(gray) lpattern(solid)) /// 
		(qfit Y  runvar [pweight = weight] if G==0 & runvar>=0 & runvar<=10, range(0 10) lcolor(gray) lpattern(solid)), ///
		 graphregion(fcol(white)) bgcolor(white) ylabel(0.5(0.5)1.5) xlabel(-10(2)10) xline(0) plotregion(style(none)) /// title("`t', `Atit'", margin(medium) size(medium)) ///
		 legend(off) ytitle("Outcome", size(large))  xtitle("Distance to cutoff", size(large))  ///
		saving(rdplot_g1_noweigth, replace)	

stop
