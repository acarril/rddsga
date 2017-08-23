 capture program drop myreg
 program myreg, eclass
     version 13.0
     tempname bb
     quietly regress mpg turn weight price
     matrix `bb'=e(b)
		 mat list `bb'
		 matrix `bb'=`bb'[1..2,1]
		 mat list `bb'
     ereturn post `bb'
     ereturn local cmd="bootstrap"
 end

 clear
 sysuse auto
 set seed 12345
 bootstrap _b, reps(50) noisily: myreg 
