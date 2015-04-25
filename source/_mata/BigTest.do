clear all
cls
set more off

use "D:\tmp\main.dta", clear
drop if missing(contrib+persona+ciiu1+raw_ubigeo)

qui adopath + "D:\Github\reghdfe\source\_hdfe"
qui adopath + "D:\Github\reghdfe\source\_common"
include reghdfe.mata

* Test mapsolve init
mata:
	void function testit(struct MapProblem scalar S) {
		drop_singletons(S.fixed_effects, S.verbose)
	}
end

local absvars contrib persona#ciiu1 raw_ubigeo

timer clear

sort contrib
ParseAbsvars `absvars', clustervars(`clustervars')
mata: S = mapsolve_init(2)
//mata: mapsolve_set(S)
set rmsg on
mata: testit(S)
set rmsg off

* Check
use "D:\tmp\main.dta", clear
drop if missing(contrib+persona+ciiu1+raw_ubigeo)

global absvars1 contrib
global absvars2 persona ciiu1
global absvars3 raw_ubigeo
sort contrib
set rmsg on
NaiveDropSingletons 3
set rmsg off

exit
