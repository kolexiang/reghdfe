clear all
cls
set more off

qui adopath + "D:\Github\reghdfe\source\_hdfe"
qui adopath + "D:\Github\reghdfe\source\_common"
include reghdfe.mata

// -------------------------------------------------------------------------------------------------
capture program drop Foobar
program define Foobar
	local absvars ciiu1 contrib#ciiu2 raw_ubigeo
	local varlist ruc date_inicio deudor
	local weightvar www
	cls

	ParseAbsvars `absvars' // , clusterva rs(`clustervars')
	mata: S = map_init("`weightvar'")
	mata: map_init_transform(S, "sym")
	mata: map_init_acceleration(S, "sd")
	mata: map_init_verbose(S, 2)
	mata: map_init_tolerance(S, 1e-8)
	mata: map_init_maxiterations(S, 1e4)

	egen mv = rowmiss(`varlist')
	drop if mv
	drop mv
	mata: map_precompute(S, tokens("foreign gear_ratio displacement turn trunk") )
	
	//preserve

	su `varlist'
	mata: map_solve(S, "`varlist'")
	regress `varlist' , nocons
	//restore
	// areg `varlist' .. , absorb(...)

end
// -------------------------------------------------------------------------------------------------

	use "D:\tmp\main.dta" if persona!=1 & estado==1 & ciiu3<. & date_inicio<=td(31dec2005) , clear
	cou
	drop if missing(contrib+ciiu1+ciiu2+raw_ubigeo)
	gen byte www = 2
	timer clear
	sort contrib
	timer on 1
	Foobar
	timer off 1
	timer list

exit
