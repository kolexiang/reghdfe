
* Parse options for Mata solver; this calls ParseAbsvars.ado and mata:map_init()
	
	ParseMAP, `options' weightvar(`weightvar') clustervars(`clustervars') // Creates Mata structure HDFE_S
	local options `s(options)'
	

* Preserve if asked to
	local newvars





* Check if we can save FEs
* BUGBUG: Just use HDFE_S.will_save_fe

	forval g = 1/`N_hdfe' {
		mata: fe2local(`g')
		local targets "`targets'`target'"
	}
	if ("`targets'"!="") {
		Assert ("`partial'"==""), msg("hdfe error: partial() not allowed when saving fixed effects")
		local numvars : word count `varlist'
		Assert `numvars'==1 , msg("hdfe error: to save the fixed effects, you need to demean only one variable")
		local opt_savefe "save_fe(1)"
	}










* Parse: absorb, clusters, and weights
	mata: map_precompute(HDFE_S)
	//mata: map_dof(HDFE_S)









	if (alguno pide grabar fe) clonevar y copy_y
	mata: map_solve(HDFE_S, varlist) // HACER EL BUNCHING INTERNAMENTE DENTRO DE SOLVE!!
	if (alguno pide grabar fe) {
		predict double resid, resid
		mata: map_recover_fe(HDFE_S, y, resid) // ver si alguno tiene target name, y en ese caso aplicar project a esos
	}


	Start, absorb(`absorb')  weight(`weighttype') weightvar(`weightvar') // BUGBUG SEE WHAT IT DOES EXACTLY AND INCORPORATE INTO INIT
	local absorb_keepvars = r(keepvars)
	local N_hdfe = r(N_hdfe) // We should know this already; maybe return it from the Parse?



* Keep relevant observations
	marksample touse, novar
	markout `touse' `varlist' `partial' `absorb_keepvars'
	qui keep if `touse'
	
* Keep relevant variables
	keep `varlist' `partial' `clustervars' `weightvar' `panelvar' `timevar' `uid' `absorb_keepvars'

* Populate Mata objects and auxiliary variables
	mata: map_precompute()
	// SAME; WE NEED ALL THESE OPTIONS INCORPORATED.. maybe split DOF into a separate ADO? or what?
	Precompute, ///
		keep(`varlist' `partial' `clustervars' `weightvar' `panelvar' `timevar' `uid') ///
		tsvars(`panelvar' `timevar') dofadjustments(pairwise clusters continuous)
	
* Compute e(df_a)
	EstimateDoF, dofadjustments(pairwise clusters continuous)
	* return list // what matters is r(kk) which will be e(df_a)
	local kk = r(kk)
	forval g = 1/`N_hdfe' {
		local df_a`g' = r(K`g') - r(M`g')
	}
	
* We don't need the FE variables (they are in mata objects now)
	*drop __FE*__

* Demean variables wrt to the fixed effects
	mata: map_solve("`varlist' `partial'")

	local opt varlist(`varlist' `partial') tol(`tolerance') maxiterations(`maxiterations') `options' `opt_savefe'
	if (`cores'>1) {
		DemeanParallel, `opt' self(hdfe) cores(`cores')
	}
	else {
		Demean, `opt'
	}

	if ("`opt_savefe'"!="") {
		Save, original_depvar(`varlist')
		local saved_fe = r(keepvars)
	}
	
	return scalar df_a = `kk'
	return scalar N_hdfe = `N_hdfe'
	forv g=1/`N_hdfe' {
		mata: fe2local(`g') // copies Mata structure into locals
		* Will inject the following with c_local:
		* ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		return local hdfe`g' = "`varlabel'"
		return scalar df_a`g' = `df_a`g'' // `levels'
	}

* Clean up Mata objects
	mata: mata drop HDFE_S

* Deal with partial() option (Alternative: do as ivreg-partial: precompute inv x'x )
BUGBUG: MAYBE ADD A PARTIAL OPTION DIRECTLY INTO MAP_SOLVE???
MAP_SOLVE(VARLIST, | PARTIAL)
NEED TO MAKE MAPSOLVE RESPECT BUNCHING

	if ("`partial'"!="") {
		tempvar resid
		_rmcoll `partial', forcedrop
		local partial = r(varlist)
		foreach var of local varlist {
			_regress `var' `partial' `weightexp' [`weighttype'`weightequal'`weightvar'], nohead notable
			_predict double `resid', resid
			qui replace `var' = `resid' // preserve labels
			drop `resid'
		}
		local numpartial : word count `partial'
		return scalar df_partial = `numpartial'
	}

	if ("`generate'"!="") {
		keep `varlist' `uid' `saved_fe'
		foreach var of local varlist {
			rename `var' `generate'`var'
		}

		tempfile output
		sort `uid'
		qui save "`output'"
		restore
		SafeMerge, uid(`uid') file("`output'") sample(`sample')
	}

}
if (_rc) {
	local rc = _rc
	cap mata: mata drop HDFE_S
	exit `rc'
}
end

* [SafeMerge: ADAPTED FROM THE ONE IN ESTIMATE.ADO]
* The idea of this program is to keep the sort order when doing the merges
cap pr drop SafeMerge
program define SafeMerge, eclass sortpreserve
syntax, uid(varname numeric) file(string) [sample(string)]
	* Merging gives us e(sample) and the FEs / AvgEs
	if ("`sample'"!="") {
		tempvar smpl
		merge 1:1 `uid' using "`file'", assert(master match) nolabel nonotes noreport gen(`smpl')
		gen byte `sample' = (`smpl'==3)
		drop `smpl' // redundant
	}
	else {
		merge 1:1 `uid' using "`file'", assert(master match) nolabel nonotes noreport nogen
	}
end



include "_common/Debug.ado"
include "_common/Version.ado"
include "_common/SortPreserve.ado"

include "_hdfe/ParseMAP.ado"
include "_hdfe/ConnectedGroups.ado"
include "_hdfe/GenerateID.ado"
include "_hdfe/AverageOthers.ado"
include "_hdfe/EstimateDoF.ado"

include "_hdfe/Start.ado"
include "_hdfe/ParseOneAbsvar.ado"
include "_hdfe/Precompute.ado"
include "_hdfe/Demean.ado"
include "_hdfe/Save.ado"
// include "_hdfe/Stop.ado"
// include "_hdfe/CheckCorrectOrder.ado"
f
