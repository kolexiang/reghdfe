capture program drop ParseAbsvars
program define ParseAbsvars, rclass
syntax anything(id="absvars" name=absvars equalok everything), [SAVEfe] // [CLUSTERvars(varlist numeric fv)]
	* Logic: split absvars -> expand each into factors -> split each into parts

	di as error "TODO: parse clustervars, return them in an r(); see how it's currently done"
	
	local g 0
	local all_cvars
	local all_ivars

	while ("`absvars'"!="") {
		local ++g
		gettoken absvar absvars : absvars, bind
		local target
		if strpos("`absvar'","=") gettoken target absvar : absvar, parse("=")
		if ("`target'"!="") {
			conf new var `target'
			gettoken eqsign absvar : absvar, parse("=")
		}

		local n : word count absvar
		local hasdot = strpos("`absvar'", ".")
		local haspound = strpos("`absvar'", "#")
		if (`n'==1 & !`hasdot' & !`haspound') local absvar i.`absvar'
		
		local 0 `absvar'
		syntax varlist(numeric fv)
			//di as error "    varlist=<`varlist'>"
		
		local ivars
		local cvars
		local has_intercept 0
		foreach factor of local varlist {
			local hasdot = strpos("`factor'", ".")
			local haspound = strpos("`factor'", "#")
			local factor_has_cvars 0

			if (!`hasdot') continue
			while ("`factor'"!="") {
				gettoken part factor : factor, parse("#")
				local is_indicator = strpos("`part'", "i.")
				local is_continuous = strpos("`part'", "c.")
				local basevar = substr("`part'", 3, .)
				if (`is_indicator') local ivars `ivars' `basevar'
				if (`is_continuous') {
					local cvars `cvars' `basevar'
					local factor_has_cvars 1
				}
			}
			if (!`factor_has_cvars') local has_intercept 1

		}
		
		local ivars : list uniq ivars
		local num_slopes : word count `cvars'
		Assert "`ivars'"!="", msg("error parsing absvars: no indicator variables in absvar <`absvar'>")
		local unique_cvars : list uniq cvars
		Assert (`: list unique_cvars == cvars'), msg("error parsing absvars: factor interactions such as i.x##i.y not allowed")

		local all_cvars `all_cvars' `cvars'
		local all_ivars `all_ivars' `ivars'

		return local target`g' `target'
		return local ivars`g' `ivars'
		return local cvars`g' `cvars'
		return scalar has_intercept`g' = `has_intercept'
		return scalar num_slopes`g' = `num_slopes'
	
		local label : subinstr local ivars " " "#"
		if (`num_slopes'==1) {
			local label `label'#c.`cvars'
		}
		else if (`num_slopes'>1) {
			local label `label'#c.(`cvars')
		}
		return local varlabel`g' `label'
	
	}
	
	local all_ivars : list uniq all_ivars
	local all_cvars : list uniq all_cvars


	return scalar G = `g'
	return scalar savefe = ("`savefe'"!="")
	return local all_ivars `all_ivars'
	return local all_cvars `all_cvars'

end



/*
cls
sysuse auto, clear
set more off
set trace off
local absvars "foreign    rep##c.(weight gear) B=turn (c.gear c.length)#turn#i.foreign trunk#turn length"
ParseAbsvars `absvars'
`Integer'	order 			
`Integer'	num_slopes
`Integer'	has_intercept
`Integer'	levels			
`Varlist'	ivars			
`Varlist'	cvars			
`Varname'	varlabel		
`Varname'	estimates		
*/
