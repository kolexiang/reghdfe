// -------------------------------------------------------------------------------------------------
// Transform data and run the regression
// -------------------------------------------------------------------------------------------------
cap pr drop Estimate
program define Estimate, eclass

**** PART II - PREPARE DATA ****

* (reporting) Store dataset size
	qui de, simple
	local old_mem = string(r(width) * r(N)  / 2^20, "%6.2f")
	local RAW_N = c(N)
	local RAW_K = c(k)

* Parse arguments and create the HDFE_S Mata structure
	Parse `0' // save all arguments into locals (verbose>=3 shows them)
	mata: verbose2local(HDFE_S, "VERBOSE")

* (optional) Preserve
	if (!`clear') {
		preserve
		Debug, level(2) newline
		Debug, level(2) msg("(dataset preserved)")
	}

* (optional) Create uid so we can then attach e(sample) and/or the Zs (the FE coefs.)
	if (!`clear' & !`fast') {
		tempvar uid
		local uid_type = cond(c(N)>c(maxlong), "double", "long")
		gen `uid_type' `uid' = _n // Useful for later merges
		la var `uid' "[UID]"
	}

* Drop unused variables
	local exp "= `weightvar'"
	if ("`weightvar'"!="") {
		la var `weightvar' "[WEIGHT] `: var label `weightvar''" // so we can distinguish it with -describe-
	}
	marksample touse, novar // Uses -if- , -in- ; -weight-? and -exp- ; can't drop any var until this
	local cluster_keepvars `clustervars'
	local cluster_keepvars : subinstr local cluster_keepvars "#" " ", all
	local cluster_keepvars : subinstr local cluster_keepvars "i." "", all
	keep `uid' `touse' `basevars' `timevar' `panelvar' `weightvar' `absorb_keepvars' `cluster_keepvars' `over'

* Expand factor and time-series variables
	local expandedvars
	local sets depvar indepvars endogvars instruments // depvar MUST be first
	Debug, level(4) newline
	Debug, level(4) msg("{title:Expanding factor and time-series variables:}")
	foreach set of local sets {
		local varlist ``set''
		if ("`varlist'"=="") continue
		local original_`set' `varlist'
		* the -if- prevents creating dummies for categories that have been excluded
		ExpandFactorVariables `varlist' if `touse', setname(`set') verbose(`VERBOSE')
		local `set' "`r(varlist)'"
		local expandedvars `expandedvars' ``set''
	}

* Drop unused basevars and tsset vars (usually no longer needed)
	if ("`vceextra'"!="") local tsvars `panelvar' `timevar' // We need to keep them only with autoco-robust VCE
	keep `uid' `touse' `expandedvars' `weightvar' `absorb_keepvars' `cluster_keepvars' `tsvars' `over' 

* Drop excluded observations and observations with missing values
	markout `touse' `expandedvars' `weightvar' `absorb_keepvars' `cluster_keepvars'
	qui keep if `touse'
	drop `touse'
	if ("`weightvar'"!="") qui drop if (`weightvar'==0)
	Assert c(N)>0, rc(2000) msg("Empty sample, check for missing values or an always-false if statement")

* Precompute Mata objects
	mata: map_init_keepvars(HDFE_S, "`expandedvars' `uid'") // Non-essential vars will be deleted (e.g. interactions of a clustervar)
	mata: map_precompute(HDFE_S)

* Replace vceoption with the correct cluster names (e.g. if it's a FE or a new variable)
	if (`num_clusters'>0) {
		assert "`r(updated_clustervars)'"!=""
		local vceoption : subinstr local vceoption "<CLUSTERVARS>" "`r(updated_clustervars)'"
	}

* (reporting) memory usage demeanings
	Debug, level(2) msg("(dataset compacted: observations " as result "`RAW_N' -> `c(N)'" as text " ; variables " as result "`RAW_K' -> `c(k)'" as text ")")
	qui de, simple
	local new_mem = string(r(width) * r(N) / 2^20, "%6.2f")
	Debug, level(2) msg("(dataset compacted, c(memory): " as result "`old_mem'" as text "M -> " as result "`new_mem'" as text "M)")
	if (`VERBOSE'>3) {
		di as text "(memory usage including mata:)"
		memory
		di as text ""
	}

* Save the statistics we need before transforming the variables
	* Compute TSS of untransformed depvar
	local tmpweightexp = subinstr("`weightexp'", "[pweight=", "[aweight=", 1)
	qui su `depvar' `tmpweightexp' // BUGBUG: Is this correct?!
	local tss = r(Var)*(r(N)-1)
	assert `tss'<.
	if (`: list posof "first" in stages') {
		foreach var of varlist `endogvars' {
			qui su `var' `tmpweightexp' // BUGBUG: Is this correct?!
			local tss_`var' = r(Var)*(r(N)-1)
		}
	}

* Compute e(df_a)
	if (!`fast') mata: store_uid(HDFE_S, "`uid'")
	mata: map_estimate_dof(HDFE_S, "`dofadjustments'", "`groupvar'")

	local N_hdfe_extended = r(N_hdfe_extended)
	local M = r(M) // FEs found to be redundant
	local M_due_to_nested = r(M_due_to_nested)
	local kk = r(df_a) // FEs that were not found to be redundant (total FEs - redundant FEs)
	Assert `kk'<.
	Assert `M'>=0 & `M'<.

	forv h=1/`N_hdfe_extended' {
		local G`h' = r(G`h')
		local M`h' = r(M`h')
		local K`h' = r(K`h')
		local M`h'_exact = r(M`h'_exact)
		local M`h'_nested = r(M`h'_nested)

		assert inlist(`M`h'_exact',0,1) // 1 or 0 whether M`h' was calculated exactly or not
		assert `M`h''<. & `K`h''<.
		assert `M`h''>=0 & `K`h''>=0
	}

* Drop the fixed effect IDs (except if they are also a clustervar or we are saving their respecting alphas)
	mata: drop_ids(HDFE_S)

* (optional) If we are saving the FEs, we need to backup the untransformed variables
	if (`will_save_fe') {
		tempfile untransformed
		qui save "`untransformed'"
	}

* (optional) Compute R2/RSS to run nested Ftests on the FEs
	* a) Compute R2 of regression without FE, to build the joint FTest for all the FEs
	* b) Also, compute RSS of regressions with less FEs so we can run nested FTests on the FEs
	if ("`model'"=="ols" & inlist("`vcetype'", "unadjusted", "ols")) {
		qui _regress `expandedvars' `weightexp', noheader notable
		local r2c = e(r2)
	}

* (optional) Compute summary statistics for the all the regression variables
	if ("`stats'"!="") {
		local tabstat_weight : subinstr local weightexp "[pweight" "[aweight"
		qui tabstat `expandedvars' `tabstat_weight' , stat(`stats') col(stat) save
		tempname statsmatrix
		matrix `statsmatrix' = r(StatTotal)
	}

* Compute residuals for all variables (overwrites vars!)
	qui ds `expandedvars'
	local NUM_VARS : word count `r(varlist)'
	Debug, msg("(computing residuals for `NUM_VARS' variables)")
	mata: map_solve(HDFE_S, "`expandedvars'")


**** PART II - REGRESSION ****

* Deal with different stages
	assert "`stages'"!=""
	assert "`stages'"=="none" // TODO
	local stage none // BUGBUG
	local lhs_endogvar "<none>" // BUGBUG

* Cleanup
	ereturn clear

* Regress
	if ("`stage'"=="none") Debug, level(2) msg("(running regresion: `model'.`ivsuite')")

	local option_list ///
		depvar indepvars endogvars instruments ///
		vceoption vcetype vcesuite ///
		kk suboptions showraw vceunadjusted first weightexp ///
		estimator twicerobust // Whether to run or not two-step gmm
	
	local options
	foreach opt of local option_list {
		if ("``opt''"!="") local options `options' `opt'(``opt'')
	}

	* Five wrappers in total, two for iv (ivreg2, ivregress), three for ols (regress, avar, mwc)
	local wrapper "Wrapper_`subcmd'" // regress ivreg2 ivregress
	if ("`subcmd'"=="regress" & "`vcesuite'"=="avar") local wrapper "Wrapper_avar"
	if ("`subcmd'"=="regress" & "`vcesuite'"=="mwc") local wrapper "Wrapper_mwc"

	if (!inlist("`stage'","none", "iv")) local wrapper "Wrapper_avar" // Compatible with ivreg2
	Debug, level(3) msg(_n "call to wrapper:" _n as result "`wrapper', `options'")
	`wrapper', `options'

* Post-regression checks
	Assert e(tss)<., msg("within tss is missing (wrapper=`wrapper')")
	local subpredict = e(predict) // used to recover the FEs
	if ("`weightvar'"!="") {
		qui su `weightvar', mean
		local sumweights = r(sum)
	}

**** PART III - RECOVER FEs AND SAVE RESULTS ****

* Replace tempnames in the coefs table
	* (e.g. __00001 -> L.somevar)
	* (this needs to be AFTER predict but before deleting IDs)
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	FixVarnames `backup_colnames'
	local newnames "`r(newnames)'"
	local prettynames "`r(prettynames)'"
	matrix colnames `b' = `newnames'

* Compute corr(FE,xb) (do before rescaling by cvar or deleting)
	// Should we compute corr with the DEMEANED or the original resids!? BUGBUG
	*if ("`model'"=="ols") {
	*	tempvar xb
	*	_predict double `xb', xb // -predict- overwrites sreturn, use _predict if needed
	*	forv g=1/`N_hdfe' { 
	*		qui corr `xb' __Z`g'__
	*		local corr`g' = r(rho)
	*	}
	*	drop `xb'
	*}

* (optional) Save FEs
	if (`will_save_fe') {
		Debug, level(2) msg("(calcualting fixed effects)")
		tempvar resid
		local score = cond("`model'"=="ols", "score", "resid")
		`subpredict' double `resid', `score' // equation: y = xb + d + e, we recovered "e"
		mata: store_resid(HDFE_S, "`resid'")
		
		qui use "`untransformed'", clear
		erase "`untransformed'"
		mata: resid2dta(HDFE_S)

		tempvar resid_d
		`subpredict' double `resid_d', `score' // This is "d+e" (including constant)

		if ("`weightvar'"!="") assert "`tmpweightexp'"!=""
		su `resid_d' `tmpweightexp', mean
		qui replace `resid_d' = `resid_d' - r(mean)
		
		tempvar d
		gen double `d' = `resid_d' - `resid'
		drop `resid' `resid_d'
		//clonevar dd = `d'
		mata: map_solve(HDFE_S, "`d'", "", "", 1) // Store FEs in Mata (will fail if partial is set)
		//regress dd __hdfe*, nocons
		drop `d'
	}

	if (`fast') Debug, msg("(option {opt fast} specified; will not save e(sample) or compute correlations)")

* (optional) Restore
	if (!`clear') {
		restore
		Debug, level(2) newline
		Debug, level(2) msg("(dataset restored)")
		// TODO: Format alphas
	}

* (optional) Save mobility groups
	if ("`groupvar'"!="") mata: groupvar2dta(HDFE_S)

* (optional) Save alphas (fixed effect estimates)
	if (`will_save_fe') mata: alphas2dta(HDFE_S)

* (optional) Add e(sample)
	if (!`fast') {
		assert !`clear'
		tempvar sample
		mata: esample2dta(HDFE_S, "`sample'")
		qui replace `sample' = 0 if `sample'==.
		la var `sample' "[HDFE Sample]"
		ereturn repost , esample(`sample')

		mata: drop_uid(HDFE_S)
	}

// PART IV - ERETURN OUTPUT

	if (`c(version)'>=12) local hidden hidden // ereturn hidden requires v12+

* Ereturns common to all commands
	ereturn local cmd = "reghdfe"
	ereturn local subcmd = cond(inlist("`stage'", "none", "iv"), "`subcmd'", "regress")
	ereturn local cmdline `"`cmdline'"'
	ereturn local model = cond("`gmm2s'"=="", "`model'", "gmm2s")
	ereturn local model = cond("`cue'"=="", "`model'", "cue")
	ereturn local model = cond("`liml'"=="", "`model'", "liml")
	ereturn local dofadjustments = "`dofadjustments'"
	ereturn local title = "HDFE " + e(title)
	ereturn local title2 =  "Absorbing `N_hdfe' HDFE " + plural(`N_hdfe', "group")
	ereturn local predict = "reghdfe_p"
	ereturn local estat_cmd = "reghdfe_estat"
	ereturn local footnote = "reghdfe_footnote"
	
	ereturn local absvars = "`original_absvars'"
	ereturn `hidden' local extended_absvars = "`extended_absvars'"

	ereturn local vcesuite = "`vcesuite'"
	ereturn local maximize_options = "`maximize_options'" // In option format; tolerance(..) etc.
	if ("`stage'"!="none") ereturn local iv_depvar = "`backup_original_depvar'"
	ereturn `hidden' local diopts = "`diopts'"
	if ("`over'"!="") {
		ereturn local over = "`over'"
		if ("`over_value'"!="") ereturn local over_value = "`over_value'"
		if ("`over_label'"!="") ereturn local over_label = "`over_label'"
		local fixed_absvars = e(absvars)
		local fixed_absvars : subinstr local fixed_absvars "i.`over'#" "", all
		local fixed_absvars : subinstr local fixed_absvars "i.`over'" "", all
		local fixed_absvars `fixed_absvars' // Trim
		ereturn local absvars = "`fixed_absvars'"
	}

	if ("`e(clustvar)'"!="") {
		ereturn local clustvar `clustervars'
		//* With kiefer/dkraay we add a time clustervar (??)
		//if ("`clustvar'"!="") ereturn local clustvar "`clustervars'"
		ereturn scalar N_clustervars = `num_clusters'
	}

	* Besides each cmd's naming style (e.g. exogr, exexog, etc.) keep one common one
	foreach cat in depvar indepvars endogvars instruments {
		local vars ``cat''
		if ("`vars'"=="") continue
		ereturn local `cat' "`original_`cat''"
	}

	ereturn `hidden' local subpredict = "`subpredict'"
	ereturn `hidden' local prettynames "`prettynames'"

	* Stata uses e(vcetype) for the SE column headers
	* In the default option, leave it empty.
	* In the cluster and robust options, set it as "Robust"
	ereturn local vcetype = proper("`vcetype'") //
	if (e(vcetype)=="Cluster") ereturn local vcetype = "Robust"
	if (e(vcetype)=="Unadjusted") ereturn local vcetype
	if ("`e(vce)'"=="." | "`e(vce)'"=="") ereturn local vce = "`vcetype'" // +-+-
	Assert inlist("`e(vcetype)'", "", "Robust", "Jackknife", "Bootstrap")

	ereturn scalar N_hdfe = `N_hdfe'
	ereturn scalar N_hdfe_extended = `N_hdfe_extended'

* Absorbed-specific returns
	ereturn scalar mobility = `M'
	ereturn scalar df_a = `kk'
	forv h=1/`N_hdfe_extended' {
		*ereturn scalar G`h' = `G`h''
		ereturn scalar M`h' = `M`h''
		ereturn scalar K`h' = `K`h''
		ereturn `hidden' scalar M`h'_exact = `M`h'_exact' // 1 or 0 whether M`h' was calculated exactly or not
		ereturn `hidden' local corr`h' = "`corr`h''" //  cond("`corr`h''"=="", ., "`corr`h''")
		ereturn `hidden' local hdfe_target`h' = "`hdfe_target`h''"
		ereturn `hidden' local hdfe_cvar`h' = "`hdfe_cvar`h''"
		ereturn `hidden' scalar M`h'_nested = `M`h'_nested'
	}

	Assert e(df_r)<. , msg("e(df_r) is missing")
	ereturn `hidden' scalar tss_within = e(tss)
	ereturn scalar tss = `tss'

	ereturn scalar ll   = -0.5 * (e(N)*ln(2*_pi) + e(N)*ln(e(rss)       /e(N)) + e(N))
	ereturn scalar ll_0 = -0.5 * (e(N)*ln(2*_pi) + e(N)*ln(e(tss_within)/e(N)) + e(N))

	ereturn scalar r2 = 1 - e(rss) / e(tss)
	ereturn scalar r2_within = 1 - e(rss) / e(tss_within)
	ereturn scalar mss = e(tss) - e(rss)

	* ivreg2 uses e(r2c) and e(r2u) for centered/uncetered R2; overwrite first and discard second
	if (e(r2c)!=.) {
		ereturn scalar r2c = e(r2)
		ereturn scalar r2u = .
	}

	* Computing Adj R2 with clustered SEs is tricky because it doesn't use the adjusted inputs:
	* 1) It uses N instead of N_clust
	* 2) For the DoFs, it uses N - Parameters instead of N_clust-1
	* 3) Further, to compute the parameters, it includes those nested within clusters
	
	* Note that this adjustment is NOT PERFECT because we won't compute the mobility groups just for improving the r2a
	* (when a FE is nested within a cluster, we don't need to compute mobilty groups; but to get the same R2a as other estimators we may want to do it)
	* Instead, you can set by hand the dof() argument and remove -cluster- from the list

	if ("`model'"=="ols" & `num_clusters'>0) Assert e(unclustered_df_r)<., msg("wtf-`vcesuite'")
	local used_df_r = cond(e(unclustered_df_r)<., e(unclustered_df_r), e(df_r)) - `M_due_to_nested'
	ereturn scalar r2_a = 1 - (e(rss)/`used_df_r') / ( e(tss) / (e(N)-1) )
	ereturn scalar rmse = sqrt( e(rss) / `used_df_r' )

	ereturn scalar r2_a_within = 1 - (e(rss)/`used_df_r') / ( e(tss_within) / (`used_df_r'+e(df_m)) )

	if (e(N_clust)<.) Assert e(df_r) == e(N_clust) - 1, msg("Error, `wrapper' should have made sure that N_clust-1==df_r")
	*if (e(N_clust)<.) ereturn scalar df_r = e(N_clust) - 1

	if ("`weightvar'"!="") ereturn scalar sumweights = `sumweights'

	if ("`model'"=="ols" & inlist("`vcetype'", "unadjusted", "ols")) {
		ereturn scalar F_absorb = (e(r2)-`r2c') / (1-e(r2)) * e(df_r) / (`kk'-1) // -1 b/c we exclude constant for this
		if (`nested') {
			local rss`N_hdfe' = e(rss)
			local temp_dof = e(N) - e(df_m) // What if there are absorbed collinear with the other RHS vars?
			local j 0
			ereturn `hidden' scalar rss0 = `rss0'
			forv g=1/`N_hdfe' {
				local temp_dof = `temp_dof' - e(K`g') + e(M`g')
				*di in red "g=`g' RSS=`rss`g'' and was `rss`j''.  dof=`temp_dof'"
				ereturn `hidden' scalar rss`g' = `rss`g''
				ereturn `hidden' scalar df_a`g' = e(K`g') - e(M`g')
				local df_a_g = e(df_a`g') - (`g'==1)
				ereturn scalar F_absorb`g' = (`rss`j''-`rss`g'') / `rss`g'' * `temp_dof' / `df_a_g'
				ereturn `hidden' scalar df_r`g' = `temp_dof'
				local j `g'
			}   
		}
	}

	// There is a big assumption here, that the number of other parameters does not increase asymptotically
	// BUGBUG: We could allow the option to indicate what parameters do increase asympt.

	if ("`savefirst'"!="") ereturn `hidden' scalar savefirst = `savefirst'

	* We have to replace -unadjusted- or else subsequent calls to -suest- will fail
	Subtitle `vceoption' // will set title2, etc. Run after e(bw) and all the others are set!
	if (e(vce)=="unadjusted") ereturn local vce = "ols"

	if ("`stages'"!="none") {
		ereturn local stage = "`stage'"
		ereturn `hidden' local stages = "`stages'"
	}

* Show table and clean up
	ereturn repost b=`b', rename // why here???

	if ("`stage'"!="none") Debug, level(0) msg(_n "{title:Stage: `stage'}" _n)
	if ("`lhs_endogvar'"!="<none>") Debug, level(0) msg("{title:Endogvar: `lhs_endogvar'}")
	Replay
	Attach, notes(`notes') statsmatrix(`statsmatrix') summarize_quietly(`summarize_quietly')

*** <<<< LAST PART OF UGLY STAGE <<<<	
if (!inlist("`stage'","none", "iv")) {
	local estimate_name reghdfe_`stage'`i_endogvar'
	local stored_estimates `stored_estimates' `estimate_name'
	local cmd estimates store `estimate_name', nocopy
	Debug, level(2) msg(" - Storing estimate: `cmd'")
	`cmd'
}
else if ("`stage'"=="iv") {
	* On the last stage, save list of all stored estimates
	assert "`stored_estimates'"!=""
	ereturn `hidden' local stored_estimates = "`stored_estimates'"
}

//} // lhs_endogvar
//} // stage

* Cleanup
	mata: mata drop HDFE_S
end
