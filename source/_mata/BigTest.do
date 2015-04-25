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
		timer_clear()
		drop_singletons(S.fixed_effects, S.verbose)
		timer()
	}
end

local absvars contrib // persona#ciiu1 // raw_ubigeo

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
NaiveDropSingletons 1
set rmsg off

exit

/*
timer report
  1.       12.2 /        5 =    2.4488   56%
  2.       3.06 /        5 =     .6112   14%
  3.       2.26 /        5 =     .4526   11%
  4.        1.7 /        5 =     .3402    8%
  5.       1.07 /        3 =  .3573333    5%
  6.       .089 /        3 =  .0296667    1%
  7.       1.13 /        3 =  .3756667    5%


  1.       12.6 /        5 =    2.5264
  2.       3.12 /        5 =      .624
  3.       2.32 /        5 =     .4644
  4.       1.74 /        5 =     .3472
  5.       1.05 /        3 =  .3516667
  6.       .082 /        3 =  .0273333
  7.       1.14 /        3 =  .3803333


Total ~ 21.6 secs
De eso, 21.5 esta contabilizado asi q ok


timer report
  1.       12.2 /        5 =     2.435
  2.          3 /        5 =     .6004
  3.       2.34 /        5 =     .4686
  4.       1.68 /        5 =     .3352
  5.       1.01 /        3 =  .3383333
  6.       .072 /        3 =      .024
  7.       1.16 /        3 =  .3866667
 11.          0 /        5 =         0
 12.        .31 /        5 =      .062
 13.          0 /        5 =         0
 14.          0 /        5 =         0
 15.       11.9 /        5 =     2.373
 21.       2.56 /        3 =      .854
 22.        .44 /        3 =  .1466667
 31.       1.33 /        5 =     .2666
 32.       1.01 /        5 =      .202

Basicamente el naive benchmark ocupa 87% del tiempo en los sorts
Nosotros tenemos un overhead asi que solo ocupamos el 53% (solo para el p, el collate es un poco mas tb)
