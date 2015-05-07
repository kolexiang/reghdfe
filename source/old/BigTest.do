clear all
cls
set more off

qui adopath + "D:\Github\reghdfe\source\_hdfe"
qui adopath + "D:\Github\reghdfe\source\_common"
include map.mata

// -------------------------------------------------------------------------------------------------
capture program drop Foobar
program define Foobar

	ParseAbsvars $absvars // , clusterva rs(`clustervars')
	mata: S = map_init()
	mata: map_init_weightvar(S, "$weightvar")
	// mata: map_init_save_ids(S, 1)
	mata: map_init_transform(S, "k") // k s c
	mata: map_init_acceleration(S, "a") // no sd a cg
	mata: map_init_verbose(S, 2)
	mata: map_init_tolerance(S, 1e-8)
	mata: map_init_maxiterations(S, 1e4)
	mata: map_init_groupsize(S, 10)
	mata: map_init_keepvars(S, "date_inscripcion ciiu1 ciiu2 contribuyente raw_ubigeo estado")
	// good matches appear to be s-cg k-ait c-cg?

	timer on 11
	mata: map_init_save_ids(S, 1) // We need them to expand the alphas into the dataset dim.
	mata: map_precompute(S)
	timer off 11
	timer on 18

	tempvar uid
	gen double `uid' = _n
	preserve
	timer off 18

	tempfile yxb
	save "`yxb'"

	// su $varlist
	timer on 12
	mata: map_solve(S, "$varlist", "", "$partial")
	timer off 12
	timer on 13
	regress $varlist , nocons
	matrix list e(b), format(%20.16g)
	timer off 13
	timer on 19

	predict double resid, resid
	keep `uid' resid
	tempfile resid
	save "`resid'"

	use "`yxb'"
	erase "`yxb'"
	merge 1:1 `uid' using "`resid'", assert(match) nolabel nonotes noreport
	predict double resid_d, resid

	// REMOVE MEAN.. NEED TO REMEMBER TO USE WEIGHTS
	su resid_d, mean
	replace resid_d = resid_d - r(mean)

	gen double d = resid_d - resid
	// keep __ID*__ d
	su resid_d d resid
	clonevar dd = d
	mata: map_init_tolerance(S, 1e-8)
	mata: map_solve(S, "d", "", "$partial", 1) // Save FE
	su __* d // Do we need to store d?
	regress dd __hdfe*__ $partial

	
	// save sample+resid
	// after restore, merge by uid
	// then 

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

	use "D:\tmp\main.dta" if date_inicio<=td(31dec2002) , clear
	cou
	drop if missing(contrib+ciiu1+ciiu2+raw_ubigeo)
	gen byte www = 1
	timer clear
	gen long t = _n
	sort contrib

	global absvars ciiu1 contrib estado, savefe // ##c.date_inscripcion, savefe // ciiu2 raw_ubigeo 
	global varlist ruc date_inicio deudor
	global partial // t // partial no es compatible con SAVEFE xq deja d+zb en vez de d
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
	reghdfe $varlist $partial, absorb($absvars) fast dof(none) tol(1e-8)
	su __*
	matrix list e(b), format(%20.16g)
	timer off 2
	timer list

	predict dd, d
	regress dd __hdfe*__

exit
