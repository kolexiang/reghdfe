cap pr drop EstimateDoF
program define EstimateDoF, rclass
syntax, [DOFadjustments(string) group(name) uid(varname) groupdta(string)]

	* Remember: fe2local stores the following:
	* ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels

* Starting point assumes no redundant parameters
	forv g=1/`G' {
		mata: fe2local(`g')
		local redundant`g' = 0 // will be 1 if we don't penalize at all for this absvar (i.e. if it's nested with cluster or collinear with another absvar)
		local is_slope`g' = ("`cvars'"!="") & (!`is_bivariate' | `is_mock') // two cases: i.a#c.b , i.a##c.b (which expands to <i.a i.a#c.b> and we want the second part)
		local M`g' = !`is_slope`g'' // Start with 0 with cont. interaction, 1 w/out cont interaction
		if (`g'==1) local M`g' = 0 // First FE has no redundant b/c it now includes the constant

		*For each FE, only know exactly parameters are redundant in a few cases:
		*i) nested in cluster, ii) first pure FE, iii) second pure FE if checked with connected groups
		local exact`g' 0
		local drop`g' = !(`is_bivariate' & `is_mock')
		local M`g'_nested = 0
	}

* Check if an absvar is a clustervar or is nested in a clustervar
* We *always* check if absvar is a clustervar, to prevent deleting its __FE__ variable by mistake
* But we only update the DoF if `adj_clusters' is true.

	local M_due_to_nested 0 // Redundant DoFs due to nesting within clusters
	if (`N_clustervars'>0) {

		forval i = 1/`N_clustervars' {
			mata: st_local("clustervar", clustervars[`i'])
			mata: st_local("clustervar_original", clustervars_original[`i'])
			local `clustervar'_original `clustervar_original'
		}

		forv g=1/`G' {
			mata: fe2local(`g')
			local gg = `g' - `is_mock'

			local absvar_is_clustervar : char __FE`gg'__[is_clustervar]
			local absvar_in_clustervar : char __FE`gg'__[in_clustervar]
			local nesting_clustervar : char __FE`gg'__[nesting_clustervar]

			if (`adj_clusters' & `absvar_is_clustervar') {
				Debug, level(1) msg("(categorical variable " as result "`varlabel'"as text " is also a cluster variable, so it doesn't count towards DoF)")
			}
			if (`adj_clusters' & `absvar_in_clustervar') {
				Debug, level(1) msg("(categorical variable " as result "`varlabel'" as text " is nested within cluster variable " as result "``clustervar'_original'" as text ", so it doesn't count towards DoF)")
			}

			if (`absvar_is_clustervar') local drop`g' 0

			if ( `adj_clusters' & (`absvar_is_clustervar' | `absvar_in_clustervar') ) {
				local M`g' = `levels' - (`g'==1) // First FE will always have at least one coef due to constant
				local redundant`g' 1
				local exact`g' 1
				local M_due_to_nested = `M_due_to_nested' + `levels' - 1
				local M`g'_nested = 1
			}
		} // end for over absvars
	} // end cluster adjustment

* Just indicate the first pure FE that is not nested in a cluster
	forv g=1/`G' {
		if (!`is_slope`g'' & !`redundant`g'') {
			local exact`g' 1
			continue, break
		}
	}

* Compute connected groups for the remaining FEs (except those with cont interactions)

	local dof_exact 0 // if this code never runs, it's not exact
	if (`adj_firstpairs' | `adj_pairwise') {
		Debug, level(3) msg(" - Calculating connected groups for DoF estimation")
		local dof_exact 1
		local i_comparison 0
		forv g=1/`G' {
			if (`is_slope`g'') local dof_exact 0 // We may not get all redundant vars with cont. interactions
			if (`is_slope`g'') continue
			local start_h = `g' + 1
			forv h=`start_h'/`G' {

				if (`is_slope`h'' | `redundant`h'') continue
				local ++i_comparison
				if (`i_comparison'>1) local dof_exact 0 // Only exact with one comparison
				if (`i_comparison'>1 & `adj_firstpairs') continue // -firstpairs- will only run the first comparison
				if (`i_comparison'==1) local exact`h' 1

				* ConnectedGroups does destructive operations and thus backups the dta by default
				* This is very slow with huge datasets and e.g. 4 FEs (up to 3*2*1=6 saves).
				* As a soln, use the -clear- opt and save before. Rule:
				* - Save the cache on the first comparison, OR if we are saving the connected group, on the second

				if (`i_comparison'==1 & "`group'"!="") {
					ConnectedGroups __FE`g'__ __FE`h'__ , gen(`group')
				}
				else if (`i_comparison'==1 & "`groupsp'"=="") | (`i_comparison'==2 & "`group'"!="") {
					tempfile backup
					qui save "`backup'"
					ConnectedGroups __FE`g'__ __FE`h'__ , clear
					qui use "`backup'", clear
				}
				else {
					ConnectedGroups __FE`g'__ __FE`h'__ , clear
					qui use "`backup'", clear
				}

				local candidate = r(groups)
				local M`h' = max(`M`h'', `candidate')
			}
		}
	} // end connected group comparisons

* Adjustment with cont. interactions
	if (`adj_continuous') {
		forv g=1/`G' {
			mata: fe2local(`g')
			if (!`is_slope`g'') continue
			CheckZerosByGroup, fe(`varname') cvars(`cvars') anyconstant(`is_mock')
			local M`g' = r(redundant)
		}
	}

	if (`dof_exact') {
		Debug, level(1) msg(" - DoF computation is exact")
	}
	else {
		Debug, level(1) msg(" - DoF computation not exact; DoF may be higher than reported")	
	}

	local SumM 0
	local SumK 0
	Debug, level(2) msg(" - Results of DoF adjustments:")
	forv g=1/`G' {
		mata: fe2local(`g')
		assert !missing(`M`g'') & !missing(`levels')
		local SumM = `SumM' + `M`g''
		local SumK = `SumK' + `levels'

		return scalar M`g' = `M`g''
		return scalar K`g' = `levels'
		return scalar M`g'_exact = `exact`g''
		return scalar M`g'_nested = `M`g'_nested'
		return scalar drop`g' = `drop`g''
		Debug, level(2) msg()
	}
	return scalar M = `SumM'
	local NetSumK = `SumK' - `SumM'
	Debug, level(2) msg(" - DoF loss due to FEs: Sum(Kg)=`SumK', M:Sum(Mg)=`SumM' --> KK:=SumK-SumM=`NetSumK'")
	return scalar kk = `NetSumK'

* Save mobility group if needed
	local saved_group = 0
	if ("`group'"!="") {
		conf var `group'
		tempfile backup
		qui save "`backup'"
		
		keep `uid' `group'
		sort `uid'
		la var `group' "Mobility group between `label'"
		qui save "`groupdta'" // A tempfile from the caller program
		Debug, level(2) msg(" - mobility group saved")
		qui use "`backup'", clear
		cap erase "`backup'"
		local saved_group = 1
	}
	return scalar saved_group = `saved_group'
	return scalar M_due_to_nested = `M_due_to_nested'
end

capture program drop CheckZerosByGroup
program define CheckZerosByGroup, rclass sortpreserve
syntax, fe(varname numeric) cvars(varname numeric) anyconstant(integer)
	tempvar redundant
	assert inlist(`anyconstant', 0, 1)
	if (`anyconstant') {
		qui bys `fe' (`cvars'): gen byte `redundant' = (`cvars'[1]==`cvars'[_N]) if (_n==1)
	}
	else {
		qui bys `fe' (`cvars'): gen byte `redundant' = (`cvars'[1]==0 & `cvars'[_N]==0) if (_n==1)
	}
	qui cou if `redundant'==1
	return scalar redundant = r(N)
end
