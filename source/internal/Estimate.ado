// -------------------------------------------------------------------------------------------------
// Transform data and run the regression
// -------------------------------------------------------------------------------------------------
cap pr drop Estimate
program define Estimate, eclass

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
	if ("`weightvar'"!="") la var `weightvar' "[WEIGHT] `: var label `weightvar''" // so we can distinguish it with -describe-
	marksample touse, novar // Uses -if- , -in- ; -weight-? and -exp- ; can't drop any var until this
	keep `uid' `touse' `basevars' `timevar' `panelvar' `weightvar' `absorb_keepvars' `clustervars' `over'

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
	keep `uid' `touse' `expandedvars' `weightvar' `absorb_keepvars' `clustervars' `tsvars' `over' 

* Drop excluded observations and observations with missing values
	markout `touse' `expandedvars' `weightvar' `absorb_keepvars' `clustervars'
	qui keep if `touse'
	drop `touse'
	if ("`weightvar'"!="") qui drop if (`weightvar'==0)
	Assert c(N)>0, rc(2000) msg("Empty sample, check for missing values or an always-false if statement")

* Precompute Mata objects
	mata: map_init_keepvars(HDFE_S, "`expandedvars' `uid'") // Non-essential vars will be deleted
	* Note: This is kinda redundant with the -keep- above except if a variable is somehow a cvar
	mata: map_precompute(HDFE_S)

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
	mata: map_estimate_dof(HDFE_S, "`dofadjustments'", "`groupvar'")
	TODO: SAVE GROUPVAR IN HDFE_S AND WAIT UNTIL RESTORE TO PUT IT BACK IN THE DTA
	local M = r(M) // FEs found to be redundant
	local M_due_to_nested = r(M_due_to_nested)
	local kk = r(df_a) // FEs that were not found to be redundant (total FEs - redundant FEs)
	Assert `kk'<.
	Assert `M'>=0 & `M'<.

	forv g=1/`N_hdfe' {
		local M`g' = r(M`g')
		local K`g' = r(K`g')
		local M`g'_exact = r(M`g'_exact)
		local M`g'_nested = r(M`g'_nested)

		assert inlist(`M`g'_exact',0,1) // 1 or 0 whether M`g' was calculated exactly or not
		assert `M`g''<. & `K`g''<.
		assert `M`g''>=0 & `K`g''>=0
		assert inlist(r(drop`g'), 0, 1)
	}

* Drop IDs for the absorbed FEs (except if its the clustervar)
* This is useful because the demeaning takes a lot of memory
	TODO: Make this a MATA call that checks whether fes[g].is_clustervar and also fes[g].target

* Replace vceoption with the correct cluster names (e.g. if it's a FE or a new variable)
	TODO
	//if (`num_clusters'>0) {
	//	mata: st_local("temp_clustervars", invtokens(clustervars))
	//	local vceoption : subinstr local vceoption "<CLUSTERVARS>" "`temp_clustervars'"
	//}





asdasdasdasd



* 12) Save untransformed data.
*	This allows us to:
*	i) do nested ftests for the FEs,
*	ii) recover the FEs, compute their correlations with xb, check that FE==1

	* We can avoid this if i) nested=check=0 ii) targets={} iii) fast=1
	mata: st_local("any_target_avge", strofreal(any(avge_target :!= "")) ) // saving avge?
	local any_target_hdfe 0 // saving hdfe?
	forv g=1/`N_hdfe' {
		mata: fe2local(`g')
		if (!`is_bivariate' | `is_mock') local hdfe_cvar`g' `cvars'
		// If it's the intercept part of the bivariate absorbed effect, don't add the cvar!
		local hdfe_target`g' `target'
		if ("`target'"!="") local any_target_hdfe 1
	}

	if (`fast') {
		if (`nested' | `check' | `any_target_hdfe' | `any_target_avge' | "`group'"!="") {
			Debug, msg(as text "(option {it:fast} not compatible with other options; disabled)") level(0)
			local fast 0
		}
		else {
			Debug, msg("(option {opt fast} specified; will not save e(sample) or compute correlations)")
		}
	}

	if (!`fast' | `cores'>1) {
		sort `uid'
		tempfile original_vars
		qui save "`original_vars'"
		if (`cores'>1) local parallel_opt `" filename("`original_vars'") uid(`uid') cores(`cores') "'
		Debug, msg("(untransformed dataset saved)") level(2)
	}

* 13) (optional) Compute R2/RSS to run nested Ftests on the FEs
	* a) Compute R2 of regression without FE, to build the joint FTest for all the FEs
	* b) Also, compute RSS of regressions with less FEs so we can run nested FTests on the FEs
	if ("`model'"=="ols" & !`savingcache') {
		qui _regress `vars' `weightexp', noheader notable
		local r2c = e(r2)

		if (`nested') {
			local rss0 = e(rss)
			local subZs
			forv g=1/`=`N_hdfe'-1' {
				Debug, msg("(computing nested model w/`g' FEs)")
				if (`cores'>1) {
					DemeanParallel, varlist(`vars') `maximize_options' num_fe(`g') self(reghdfe) `parallel_opt'
				}
				else {
					Demean, varlist(`vars') `maximize_options' num_fe(`g')	
				}

				qui _regress `vars' `weightexp', noheader notable
				local rss`g' = e(rss)
				qui use "`original_vars'", clear // Back to untransformed dataset
			}
		}
	}

	* Get normalized string of the absvars (i.e. turn -> i.turn)
	local original_absvars
	forv g=1/`N_hdfe' {
		mata: fe2local(`g')
		local original_absvars `original_absvars'  `varlabel'
	}

* Compute summary statistics for the all the regression variables
	if ("`stats'"!="") {
		local tabstat_weight : subinstr local weightexp "[pweight" "[aweight"
		qui tabstat `vars' `tabstat_weight' , stat(`stats') col(stat) save
		tempname statsmatrix
		matrix `statsmatrix' = r(StatTotal)
	}

* 14) Compute residuals for all variables including the AvgEs (overwrites vars!)
	qui ds `vars'
	local NUM_VARS : word count `r(varlist)'
	Debug, msg("(computing residuals for `NUM_VARS' variables)")
	Debug, msg(" - tolerance = `tolerance'")
	Debug, msg(" - max. iter = `maxiterations'")
	if ("`usecache'"=="") {
		if (`cores'>1) {
			DemeanParallel, varlist(`vars') `maximize_options' self(reghdfe) `parallel_opt'
		}
		else {
			Demean, varlist(`vars') `maximize_options'	
		}
	}
	else {
		Debug, msg("(using cache data)")
		drop `vars'
		local handshake_master : char __uid__[handshake]
		char __uid__[handshake]
		// An error in the merge most likely means different # of obs due to missing values in a group but not in other
		// try with if !missing(__uid__) // TODO: Auto-add this by default?
		// TODO: Make this fool-proof when using -over-
		if ("`over'"!="") local using using // This is dangerous
		sort __uid__ // The user may have changed the sort order of the master data
		qui merge 1:1 __uid__ using "`usecache'", keepusing(`vars') assert(match master `using') keep(master match) nolabel sorted
		qui cou if _merge!=3
		Assert r(N)==0, msg(as error "Error: the cache has `r(N)' less observations than the master data" _n ///
			as text " - This is possibly because, when created, it included variables that were missing in cases where the current ones are not.")
		qui drop if _merge!=3
		drop _merge

		local handshake_using : char __uid__[handshake]
		local tolerance_using : char __uid__[tolerance]
		local maxiterations_using : char __uid__[maxiterations]
		Assert (`handshake_master'==`handshake_using'), msg("using dataset does not have the same __uid__")
		Assert abs(`tolerance'-`tolerance_using')<epsdouble(), msg("using dataset not computed with the same tolerance (`tolerance_using')")
		Assert (`maxiterations'==`maxiterations_using'), msg("using dataset not computed with the same maxiterations (`maxiterations_using')")

		local absvar_master `original_absvars'
		local absvar_using : char __uid__[absvars_key]
		Assert ("`absvar_master'"=="`absvar_using'"), msg("using dataset not created with the same absvars")
		char __uid__[absvars_key]
	}

if (`savingcache') {
	Debug, msg("(saving cache and exiting)")
	char __uid__[absvars_key] `original_absvars'
	sort __uid__
	save "`savecache'", replace
	return clear
	ereturn clear
	ereturn local cmdline `"`cmdline'"'
	if ("`over_levels'"!="") ereturn local over_levels = "`over_levels'"
	exit
}

// PART II - REGRESSION

**** <<< START OF UGLY -stages- CODE
assert "`stages'"!=""
if ("`stages'"!="none") {
	Debug, level(2) msg(_n " {title:Stages to run}: " as result "`stages'" _n)
	local backup_fast `fast'
	local num_stages : word count `stages'
	local last_stage : word `num_stages' of `stages'
	assert "`last_stage'"=="iv"
	foreach vargroup in depvar indepvars endogvars instruments {
		local backup_`vargroup' ``vargroup''
		local backup_original_`vargroup' `original_`vargroup''
	}
	local backup_tss = `tss'
}

foreach stage of local stages {
local lhs_endogvars = cond("`stage'"=="first", "`backup_endogvars'", "<none>")

if ("`stage'"=="first") {
	local i_endogvar 0
}
else {
	local i_endogvar
}

foreach lhs_endogvar of local lhs_endogvars {
Assert inlist("`stage'", "none", "iv", "first", "ols", "reduced", "acid")

if ("`stage'"=="iv") {
	local tss = `backup_tss'
	local fast `backup_fast'
	local depvar `backup_depvar'
	local indepvars `backup_indepvars'
	local endogvars `backup_endogvars'
	local instruments `backup_instruments'
	local original_depvar `backup_original_depvar'
	local original_indepvars `backup_original_indepvars'
	local original_endogvars `backup_original_endogvars'
	local original_instruments `backup_original_instruments'
}
else if ("`stage'"=="ols") {
	local tss = `backup_tss'
	local fast 1
	local depvar `backup_depvar'
	local endogvars
	local indepvars `backup_indepvars' `backup_endogvars'
	local instruments
	local original_depvar `backup_original_depvar'
	local original_endogvars
	local original_indepvars `backup_original_indepvars' `backup_original_endogvars'
	local original_instruments
	local vcesuite avar
}
else if ("`stage'"=="reduced") {
	local tss = `backup_tss'
	local fast 1
	local depvar `backup_depvar'
	local indepvars `backup_indepvars' `backup_instruments'
	local endogvars
	local instruments
	local original_depvar `backup_original_depvar'
	local original_indepvars `backup_original_indepvars' `backup_original_instruments'
	local original_endogvars
	local original_instruments
	local vcesuite avar
}
else if ("`stage'"=="acid") {
	local tss = `backup_tss'
	local fast 1
	local depvar `backup_depvar'
	local indepvars `backup_indepvars' `backup_endogvars' `backup_instruments'
	local endogvars
	local instruments
	local original_depvar `backup_original_depvar'
	local original_indepvars `backup_original_indepvars' `backup_original_endogvars' `backup_original_instruments'
	local original_endogvars
	local original_instruments
	local vcesuite avar
}
else if ("`stage'"=="first") {
	local ++ i_endogvar
	local tss = `tss_`lhs_endogvar''
	local fast 1
	local depvar `lhs_endogvar'
	local indepvars `backup_indepvars' `backup_instruments'
	local endogvars
	local instruments
	local original_depvar : word `i_endogvar' of `backup_original_endogvars'
	local original_indepvars `backup_original_indepvars' `backup_original_endogvars' `backup_original_instruments'
	local original_endogvars
	local original_instruments
	local vcesuite avar
}
**** END OF UGLY -stages- CODE >>>> 

* Cleanup
	ereturn clear

* Regress
	if ("`stage'"=="none") Debug, level(2) msg("(running regresion: `model'.`ivsuite')")
	local avge = cond(`N_avge'>0, "__W*__", "")
	local options
	local option_list ///
		depvar indepvars endogvars instruments avgevars ///
		original_depvar original_indepvars original_endogvars ///
		original_instruments original_absvars avge_targets ///
		vceoption vcetype vcesuite ///
		kk suboptions showraw vceunadjusted first weightexp ///
		estimator twicerobust // Whether to run or not two-step gmm
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
	
	Assert e(tss)<., msg("within tss is missing (wrapper=`wrapper')")
	
	local subpredict = e(predict) // used to recover the FEs

	if ("`weightvar'"!="") {
		qui su `weightvar', mean
		local sumweights = r(sum)
	}

// PART III - RECOVER FEs AND SAVE RESULTS 

if (`fast') {
	* Copy pasted from below
	Debug, level(3) msg("(avoiding -use- of temporary dataset)")
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	FixVarnames `backup_colnames'
	local newnames "`r(newnames)'"
	local prettynames "`r(prettynames)'"
	matrix colnames `b' = `newnames'
}
else {
	assert inlist("`stage'", "iv", "none")
* 1) Restore untransformed dataset
	qui use "`original_vars'", clear

* 2) Recover the FEs

	* Predict will get (e+d) from the equation y=xb+d+e
	tempvar resid_d
	if e(df_m)>0 {
		local score = cond("`model'"=="ols", "score", "resid")
		`subpredict' double `resid_d', `score' // Auto-selects the program based on the estimation method		
	}
	else {
		gen double `resid_d' = `depvar'
	}

	** If the eqn doesn't have a constant, we need to save the mean of the resid in order to add it when predicting xb
	*if (!`addconstant') {
	*	su `resid_d', mean
	*	ereturn `hidden' scalar _cons = r(mean)
	*}

	Debug, level(2) msg("(loaded untransformed variables, predicted residuals)")

	* Absorb the residuals to obtain the FEs (i.e. run a regression on just the resids)
	Debug, level(2) tic(31)
	Demean, varlist(`resid_d') `maximize_options' save_fe(1)
	Debug, level(2) toc(31) msg("mata:make_residual on final model took")
	drop `resid_d'

* 3) Compute corr(FE,xb) (do before rescaling by cvar or deleting)
	if ("`model'"=="ols") {
		tempvar xb
		_predict double `xb', xb // -predict- overwrites sreturn, use _predict if needed
		forv g=1/`N_hdfe' { 
			qui corr `xb' __Z`g'__
			local corr`g' = r(rho)
		}
		drop `xb'
	}

* 4) Replace tempnames in the coefs table
	* (e.g. __00001 -> L.somevar)
	* (this needs to be AFTER predict but before deleting FEs and AvgEs)
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	FixVarnames `backup_colnames'
	local newnames "`r(newnames)'"
	local prettynames "`r(prettynames)'"
	matrix colnames `b' = `newnames'

* 5) Save FEs w/proper name, format
	Save, original_depvar(`original_depvar')
	local keepvars `r(keepvars)'
	if ("`keepvars'"!="") format `fe_format' `keepvars'
	
* 6) Save AvgEs
	forv g=1/`N_avge' {
		local var __W`g'__
		local target : char `var'[target]
		if ("`target'"!="") {
			rename `var' `target'
			local avge_target`g' `target' // Used by -predict-
			local keepvars `keepvars' `target'
		}
	}

	if ("`keepvars'"!="") format `fe_format' `keepvars' // The format of depvar, saved by -Parse-

* 7) Save dataset with FEs and e(sample)
	keep `uid' `keepvars'
	tempfile output
	qui save "`output'"
} // fast

* 8) Restore original dataset and merge
	if (inlist("`stage'","none", "iv")) restore // Restore user-provided dataset (since -iv- comes at the end, that is done at that stage!)
	if (!`fast') {
		// `saved_group' was created by EstimateDoF.ado
		if (!`saved_group')  local groupdta
		SafeMerge, uid(`uid') file("`output'") groupdta("`groupdta'")
		*cap tsset, noquery // we changed -sortby- when we merged (even if we didn't really resort)
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
	ereturn local title2 =  "Absorbing `N_hdfe' HDFE " + plural(`N_hdfe', "indicator")
	ereturn local predict = "reghdfe_p"
	ereturn local estat_cmd = "reghdfe_estat"
	ereturn local footnote = "reghdfe_footnote"
	ereturn local absvars = "`original_absvars'"
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
		mata: st_local("clustvar", invtokens(clustervars_original))
		* With kiefer/dkraay we add a time clustervar
		if ("`clustvar'"!="") ereturn local clustvar "`clustvar'"
		ereturn scalar N_clustervars = `num_clusters'
	}

	* Besides each cmd's naming style (e.g. exogr, exexog, etc.) keep one common one
	foreach cat in depvar indepvars endogvars instruments {
		local vars ``cat''
		if ("`vars'"=="") continue
		ereturn local `cat' "`original_`cat''"
	}
	ereturn local avgevars "`avge'" // bugbug?

	ereturn `hidden' local subpredict = "`subpredict'"
	ereturn `hidden' local prettynames "`prettynames'"
	forv g=1/`N_avge' {
		ereturn `hidden' local avge_target`g' "`avge_target`g''" // Used by -predict-
	}

	* Stata uses e(vcetype) for the SE column headers
	* In the default option, leave it empty.
	* In the cluster and robust options, set it as "Robust"
	ereturn local vcetype = proper("`vcetype'") //
	if (e(vcetype)=="Cluster") ereturn local vcetype = "Robust"
	if (e(vcetype)=="Unadjusted") ereturn local vcetype
	if ("`e(vce)'"=="." | "`e(vce)'"=="") ereturn local vce = "`vcetype'" // +-+-
	Assert inlist("`e(vcetype)'", "", "Robust", "Jackknife", "Bootstrap")

	ereturn scalar N_hdfe = `N_hdfe'
	if ("`N_avge'"!="") ereturn scalar N_avge = `N_avge'

* Absorbed-specific returns
	ereturn scalar mobility = `M'
	ereturn scalar df_a = `kk'
	forv g=1/`N_hdfe' {
		ereturn scalar M`g' = `M`g''
		ereturn scalar K`g' = `K`g''
		ereturn `hidden' scalar M`g'_exact = `M`g'_exact' // 1 or 0 whether M`g' was calculated exactly or not
		ereturn `hidden' local corr`g' = "`corr`g''" //  cond("`corr`g''"=="", ., "`corr`g''")
		ereturn `hidden' local hdfe_target`g' = "`hdfe_target`g''"
		ereturn `hidden' local hdfe_cvar`g' = "`hdfe_cvar`g''"
		ereturn `hidden' scalar M`g'_nested = `M`g'_nested'
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

} // lhs_endogvar
} // stage
*** >>>> LAST PART OF UGLY STAGE >>>>

	Stop

end

// -------------------------------------------------------------------------------------------------

* The idea of this program is to keep the sort order when doing the merges
cap pr drop SafeMerge
program define SafeMerge, eclass sortpreserve
syntax, uid(varname numeric) file(string) [groupdta(string)]
	* Merging gives us e(sample) and the FEs / AvgEs
	tempvar merge
	merge 1:1 `uid' using "`file'", assert(master match) nolabel nonotes noreport gen(`merge')
	
	* Add e(sample) from _merge
	tempvar sample
	gen byte `sample' = (`merge'==3)
	la var `sample' "[HDFE Sample]"
	ereturn repost , esample(`sample')
	drop `merge'

	* Add mobility group
	if ("`groupdta'"!="") merge 1:1 `uid' using "`groupdta'", assert(master match) nogen nolabel nonotes noreport sorted
end

capture program drop Subtitle
program define Subtitle, eclass
	* Fill e(title3/4/5) based on the info of the other e(..)

	if (inlist("`e(vcetype)'", "Robust", "Cluster")) local hacsubtitle1 "heteroskedasticity"
	if ("`e(kernel)'"!="" & "`e(clustvar)'"=="") local hacsubtitle3 "autocorrelation"
	if ("`e(kiefer)'"!="") local hacsubtitle3 "within-cluster autocorrelation (Kiefer)"
	if ("`hacsubtitle1'"!="" & "`hacsubtitle3'" != "") local hacsubtitle2 " and "
	local hacsubtitle "`hacsubtitle1'`hacsubtitle2'`hacsubtitle3'"
	if strlen("`hacsubtitle'")>30 {
		local hacsubtitle : subinstr local hacsubtitle "heteroskedasticity" "heterosk.", all word
		local hacsubtitle : subinstr local hacsubtitle "autocorrelation" "autocorr.", all word
	}
	if ("`hacsubtitle'"!="") {
		ereturn local title3 = "Statistics robust to `hacsubtitle'"
		
		if ("`e(kernel)'"!="") local notes " `notes' kernel=`e(kernel)'"
		if ("`e(bw)'"!="") local notes " `notes' bw=`e(bw)'"
		if ("`e(dkraay)'"!="") local notes " `notes' dkraay=`e(dkraay)'"
		local notes `notes' // remove initial space
		if ("`notes'"!="") ereturn local title4 = " (`notes')"
		if ("`notes'"!="") {
			if ("`_dta[_TSpanel]'"!="") local tsset panel=`_dta[_TSpanel]'
			if ("`_dta[_TStvar]'"!="") local tsset `tsset' time=`_dta[_TStvar]'
			local tsset `tsset'
			ereturn local title5 = " (`tsset')"
		}
	}
end
