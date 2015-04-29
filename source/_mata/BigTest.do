clear all
cls
set more off

qui adopath + "D:\Github\reghdfe\source\_hdfe"
qui adopath + "D:\Github\reghdfe\source\_common"
include reghdfe.mata

// -------------------------------------------------------------------------------------------------
capture program drop Foobar
program define Foobar

	ParseAbsvars $absvars // , clusterva rs(`clustervars')
	mata: S = map_init("`weightvar'")
	mata: map_init_transform(S, "s") // k s c
	mata: map_init_acceleration(S, "cg") // no sd a cg
	mata: map_init_verbose(S, 2)
	mata: map_init_tolerance(S, 1e-6)
	mata: map_init_maxiterations(S, 1e4)
	// good matches appear to be s-cg k-ait c-cg?

	timer on 11
	mata: map_precompute(S, tokens("date_inscripcion foreign gear_ratio displacement turn trunk ciiu1 ciiu2 contribuyente raw_ubigeo") )
	timer off 11
	timer on 18
	preserve
	timer off 18

	su $varlist
	timer on 12
	mata: map_solve(S, "$varlist")
	timer off 12
	timer on 13
	regress $varlist , nocons
	matrix list e(b), format(%20.16g)
	timer off 13
	timer on 19
	restore
	timer off 19
	// areg $varlist .. , absorb(...)

end
// -------------------------------------------------------------------------------------------------
// TODO:
// STOPPING CRITERIA FOR CG
// EMBED INTO REGHDFE, ADD OPTIONS, ETC
// PROFILE
// CHECK FOR CORRECTNESS, MAYBE FORCE QUADSUM/QUADROWSUM 


	use "D:\tmp\main.dta"  /*if persona!=1 & estado==1 & ciiu3<. & date_inicio<=td(31dec2006)*/ , clear
	cou
	drop if missing(contrib+ciiu1+ciiu2+raw_ubigeo)
	gen byte www = 2
	timer clear
	sort contrib

	global absvars ciiu1 contrib#c.date_inscripcion // ciiu2 raw_ubigeo 
	global varlist ruc date_inicio deudor
	local weightvar www
	cls

	drop if date_inscripcion==.
	egen mv = rowmiss($varlist)
	drop if mv
	drop mv

	timer on 1
	Foobar
	timer off 1
	timer on 2
	reghdfe $varlist, absorb($absvars) fast dof(none) tol(1e-6)
	matrix list e(b), format(%20.16g)
	timer off 2
	timer list

exit
