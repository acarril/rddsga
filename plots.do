use rddsga_synth, clear

// Color scheme
*ssc install scheme-burd, replace
set scheme burd

gen bin=floor(runvar/2)*2+1


keep if abs(runvar)<=10

*** Bin mean (all)

rd Y runvar, bw(10)

    bysort bin: egen bin_Y = mean(Y) 

    graph twoway  (scatter bin_Y  bin,  msymbol(O) mcolor(gray) msize(medium)) ///
    (qfitci Y runvar if  runvar>=-10 & runvar<=0, range(-10 0) lcolor(black) lpattern(solid)) /// 
    (qfitci Y runvar if  runvar>=0 & runvar<=10, range(0 10) lcolor(black) lpattern(solid)), ///
     ylabel(0.5(0.5)1.5) xlabel(-10(2)10) xline(0) /// title("`t', `Atit'", margin(medium) size(medium)) ///
     scheme(burd) legend(off) ytitle("Outcome")  xtitle("Distance to cutoff")  ///
    saving(figs/rdplot_unw, replace)
    graph export figs/rdplot_unw.png, replace width(480)
  * Bin means by group:

    bysort bin G: egen bin_mean_Y = mean(Y) 

    *** MEDIA PONDERADA
  
  * Bin means:

 bysort bin G: egen bin_mean_Yw = wtmean(Y), weight(weight) 
    
  * Draw plot:
sort runvar

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
  (qfit Y runvar if  G & runvar>=-10 & runvar<= 0, range(-10 0) lwidth(medthick) lcolor(`c1')) /// 
  (qfit Y runvar if  G & runvar>=  0 & runvar<=10, range( 0 10) lwidth(medthick) lcolor(`c1'))  ///
  (qfit Y runvar if !G & runvar>=-10 & runvar<= 0, range(-10 0) lwidth(medthick) lcolor(`c2')) /// 
  (qfit Y runvar if !G & runvar>=  0 & runvar<=10, range( 0 10) lwidth(medthick) lcolor(`c2')) ///
  (qfit Y runvar [pweight = weight] if  G & runvar>=-10 & runvar<= 0, range(-10 0) lwidth(medthick) lcolor(`c3')) /// 
  (qfit Y runvar [pweight = weight] if  G & runvar>=  0 & runvar<=10, range( 0 10) lwidth(medthick) lcolor(`c3')) ///
  (qfit Y runvar [pweight = weight] if !G & runvar>=-10 & runvar<= 0, range(-10 0) lwidth(medthick) lcolor(`c4')) /// 
  (qfit Y runvar [pweight = weight] if !G & runvar>=  0 & runvar<=10, range( 0 10) lwidth(medthick) lcolor(`c4')), ///
    scheme(burd) ///
    ylabel(0.5(0.5)1.5) xlabel(-10(2)10) xline(0) /// title("`t', `Atit'", margin(medium) size(medium)) ///
    ytitle("Outcome") xtitle("Distance to cutoff")  ///
    legend(on order(5 "Group 1 (unweighted)" 7 "Group 0 (unweighted)" 9 "Group 1 (PSW)" 11 "Group 0 (PSW)" ) position(6)) ///
    saving(figs/rdplot_allw_bygroup, replace)
  graph export figs/rdplot_allw_bygroup.png, replace width(480)
