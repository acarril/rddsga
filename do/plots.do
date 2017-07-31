// Set root directory
if "`c(os)'" == "MacOSX" cd "/Users/alvaro/Library/Application Support/Stata/ado/personal/rddsga/"
if "`c(os)'" == "Windows" cd "C:\ado\personal\rddsga"

use rddsga_synth, clear

// Color scheme
*ssc install scheme-burd, replace
set scheme burd

gen bin=floor(Z/2)*2+1


keep if abs(Z)<=10

*** Bin mean (all)

rd Y Z, bw(10)

    bysort bin: egen bin_Y = mean(Y) 

    graph twoway  ///
    (qfitci Y Z if  Z>=-10 & Z<=0, range(-10 0) lcolor(black) lpattern(solid)) /// 
    (qfitci Y Z if  Z>=0 & Z<=10, range(0 10) lcolor(black) lpattern(solid)), ///
     xlabel(,grid) xline(0) /// title("`t', `Atit'", margin(medium) size(medium)) ///
     scheme(burd) legend(off) ylabel(-1(1)2) ytitle("Outcome")  xtitle("Distance to cutoff")  ///
    saving(figs/rdplot_unw, replace)
    graph export figs/rdplot_unw.png, replace width(960)
  * Bin means by group:


  
    bysort bin G: egen bin_mean_Y = mean(Y) 

    *** MEDIA PONDERADA
  
  * Bin means:

 bysort bin G: egen bin_mean_Yw = wtmean(Y), weight(weight) 
    
  * Draw plot:
sort Z

local c1 eltblue
local c2 erose
local c3 edkblue
local c4 maroon

local scatter_opts msymbol(Oh)

graph twoway ///
  (scatter bin_mean_Y   bin if  G, mcolor(`c1') `scatter_opts') ///
  (scatter bin_mean_Y   bin if !G, mcolor(`c2') `scatter_opts') ///
  (scatter bin_mean_Yw  bin if  G, mcolor(`c3') `scatter_opts') ///
  (scatter bin_mean_Yw  bin if !G, mcolor(`c4') `scatter_opts') ///
  (qfit Y Z if  G & Z>=-10 & Z<= 0, range(-10 0) lwidth(medthick) lcolor(`c1')) /// 
  (qfit Y Z if  G & Z>=  0 & Z<=10, range( 0 10) lwidth(medthick) lcolor(`c1'))  ///
  (qfit Y Z if !G & Z>=-10 & Z<= 0, range(-10 0) lwidth(medthick) lcolor(`c2')) /// 
  (qfit Y Z if !G & Z>=  0 & Z<=10, range( 0 10) lwidth(medthick) lcolor(`c2')) ///
  (qfit Y Z [pweight = weight] if  G & Z>=-10 & Z<= 0, range(-10 0) lwidth(medthick) lcolor(`c3')) /// 
  (qfit Y Z [pweight = weight] if  G & Z>=  0 & Z<=10, range( 0 10) lwidth(medthick) lcolor(`c3')) ///
  (qfit Y Z [pweight = weight] if !G & Z>=-10 & Z<= 0, range(-10 0) lwidth(medthick) lcolor(`c4')) /// 
  (qfit Y Z [pweight = weight] if !G & Z>=  0 & Z<=10, range( 0 10) lwidth(medthick) lcolor(`c4')), ///
    scheme(burd) ///
 xlabel(,grid) xline(0) /// title("`t', `Atit'", margin(medium) size(medium)) ///
    ytitle("Outcome") xtitle("Distance to cutoff")  ///
    legend(on order(5 "Group 1 (unweighted)" 7 "Group 0 (unweighted)" 9 "Group 1 (PSW)" 11 "Group 0 (PSW)" ) position(6)) ///
    saving(figs/rdplot_allw_bygroup, replace)
  graph export figs/rdplot_allw_bygroup.png, replace width(960)
