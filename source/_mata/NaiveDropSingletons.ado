// Benchmark Code
* INPUT: Number of absvars G
* It needs $absvar# to exist

capture program drop NaiveDropSingletons
program define NaiveDropSingletons
	args G
	gen byte singleton = 0
	gen byte delta = 0

	local i_last_singleton = 1
	local i = 1
	local g = 1
	
	while (`i'<`i_last_singleton'+`G') {
		if (`g'>`G') local g 1
		di as text "i=`i' g=`g'/`G'"
		local id = cond(`i'<=`G', "\${absvars`g'}", "id`g'")
		
		sort `id' // Inefficient for every i above G because we'll repeat ourselves
		qui by `id': replace delta = (_n==1)
		qui by `id': replace singleton = (_N==1) // faster: replace singleton = (delta==1 & delta[_n+1]==1)
		qui cou if singleton
		if r(N)>0 {
			di as text "  (`r(N)' singletons deleted)"
			local i_last_singleton `i'
			qui drop if singleton
		}
		else {
			di as text "  (no singletons)"
		}

		if (`i'<=`G') {
			gen long id`g' = sum(delta)
		}
		else if (`i_last_singleton'<`i') {
			qui replace id`g' = sum(delta) // only update if there is a change this is the real thing
		}

		local ++i
		local ++g
		if (`i'>100) error 123
	}
	drop singleton delta
end

exit

* Test Code
set more off
cls
clear
sysuse auto
drop if missing(rep)
global absvars1 rep foreign
global absvars2 turn trunk
global absvars3 rep
global absvars4 trunk foreign

NaiveDropSingletons 4

tab1 id*, m
cou

exit


local absvars rep#foreign   turn#rep##c.(weight gear) // B=turn (c.gear c.weight)#turn
//local absvars "foreign##c.weight"
//local absvars rep#foreign#turn

ParseAbsvars `absvars', clustervars(`clustervars')
mata: S = mapsolve_init(2)
mata: mapsolve_set(S)
mata: S.G
