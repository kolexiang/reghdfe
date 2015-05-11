capture program drop InnerSaveCache
program define InnerSaveCache, eclass
* (note: based on Inner.ado)

* INITIAL CLEANUP
	ereturn clear // Clear previous results and drops e(sample)

* PARSE - inject opts with c_local, create Mata structure HDFE_S (use verbose>2 for details)
	Parse `0'
	assert `savecache'
	Assert !`will_save_fe', msg("savecache disallows saving FEs")

* CREATE UID - allows attaching e(sample) and the FE estimates into the restored dataset
	if (`preserve' & !`fast') {
		tempvar uid
		GenUID `uid'
	}

* COMPACT - Expand time and factor variables, and drop unused variables and obs.
	local original_depvar "`depvar'"
	Compact, basevars(`basevars') depvar(`depvar') indepvars(`indepvars') endogvars(`endogvars') instruments(`instruments') uid(`uid') timevar(`timevar') panelvar(`panelvar') weightvar(`weightvar') absorb_keepvars(`absorb_keepvars') clustervars(`clustervars') if(`if') in(`in') by(`by') verbose(`verbose') vceextra(`vceextra')
	// Injects locals: depvar indepvars endogvars instruments expandedvars

* PRECOMPUTE MATA OBJECTS (means, counts, etc.)
	mata: map_init_keepvars(HDFE_S, "`expandedvars' `uid'") 	// Non-essential vars will be deleted (e.g. interactions of a clustervar)
	mata: map_precompute(HDFE_S)
	global updated_clustervars = "`r(updated_clustervars)'"
	
* PREPARE - Compute untransformed tss, R2 of eqn w/out FEs
!!!BUGBUG TENGO QUE GUARDAR PA TODOS?!
	Prepare, weightexp(`weightexp') depvar(`depvar') stages(`stages') model(`model') expandedvars(`expandedvars') vcetype(`vcetype') endogvars(`endogvars')
	* Injects tss, tss_`endogvar' (with stages), and r2c

* COMPUTE e(stats) - Summary statistics for the all the regression variables
	if ("`stats'"!="") {
		local tabstat_weight : subinstr local weightexp "[pweight" "[aweight"
		qui tabstat `expandedvars' `tabstat_weight' , stat(`stats') col(stat) save
		tempname statsmatrix
		matrix `statsmatrix' = r(StatTotal)
	}

* MAP_SOLVE() - WITHIN TRANFORMATION (note: overwrites variables)
	qui ds `expandedvars'
	local NUM_VARS : word count `r(varlist)'
	Debug, msg("(computing residuals for `NUM_VARS' variables)")
	mata: map_solve(HDFE_S, "`expandedvars'")
end
