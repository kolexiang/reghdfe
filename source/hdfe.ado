capture program drop hdfe
program define hdfe, rclass

* Set Stata version
	version `=clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred

* Intercept version calls
	cap syntax, version
	if !c(rc) {
		Version
		exit
	}

* Parse
	syntax varlist(numeric) [if] [in] [fw aw pw/] , ///
	/// Main Options ///
		Absorb(string) ///
		[ ///
		PARTIAL(varlist numeric) /// Additional regressors besides those in absorb()
		SAMPLE(name) ///
		Generate(name) CLEAR /// Replace dataset, or just add new variables
		CLUSTERVARs(varlist numeric fv max=10) /// Used to estimate the DoF
	/// Optimization ///
		GROUPsize(string) /// Process variables in groups of #
		TRANSFORM(string) ///
		ACCELeration(string) ///
		Verbose(string) ///
		TOLerance(string) ///
		MAXITerations(string) ///
		] [*] // Remaining options 

* Time/panel variables
	local timevar `_dta[_TStvar]'
	local panelvar `_dta[_TSpanel]'

* Validation
	local clustervars : subinstr local clustervars "i." "", all // Remove i. prefixes
	if ("`options'"!="") di as error "unused options: `options'"
	if ("`sample'"!="") confirm new variable `sample'
	Assert ("`generate'"!="") + ("`clear'"!="") == 1 , ///
		msg("hdfe error: you need to specify one and only one of the following options: clear generate(...)")
	Assert "`: list varlist & partial'"=="", ///
		msg("variables in varlist cannot appear in partial()")
	if ("`weight'"!="") {
		local weightvar `exp'
		local weighttype `weight'
		confirm var `weightvar', exact // just allow simple weights
	}

* From now on, we will pollute the Mata workspace, so wrap this in case of error
	cap noi {

* Create Mata structure
	ParseAbsvars `absorb' // Stores results in r()
	mata: HDFE_S = map_init() // Reads results from r()

	if ("`weightvar'"!="") mata: map_init_weights(HDFE_S, "`weightvar'", "`weighttype'")
	* String options
	local optlist transform acceleration clustervars panelvar timevar
	foreach opt of local optlist {
		if ("``opt''"!="") mata: map_init_`opt'(HDFE_S, "``opt''")
	}
	* Numeric options
	local optlist groupsize verbose tolerance maxiterations
	foreach opt of local optlist {
		if ("``opt''"!="") mata: map_init_`opt'(HDFE_S, ``opt'')
	}

* (Optional) Preserve
	if ("`generate'"!="") {
		foreach var of varlist `varlist' {
			confirm new var `generate'`var'
			local newvars `newvars' `generate'`var'
		}
		tempvar uid
		local uid_type = cond(c(N)>c(maxlong), "double", "long")
		gen `uid_type' `uid' = _n // Useful for later merges
		la var `uid' "[UID]"
		
		preserve
	}

* Precompute Mata objects
	mata: map_init_keepvars(HDFE_S, "`varlist' `partial' `uid'") // Non-essential vars will be deleted
	mata: map_precompute(HDFE_S)
	de

* Compute e(df_a)
	// EstimateDoF, dofadjustments(pairwise clusters continuous)

* Within Transformation
	//mata: map_solve(HDFE_S, "`varlist'", "`newvars'", "`partial'")

* (Optional) Save FEs

* (Optional) Tempsave, restore and merge with 

* Cleanup after an error
	} // cap noi
	if c(rc) {
		local rc = c(rc)
		cap mata: mata drop HDFE_S // overwrites c(rc)
		exit `rc'
	}
end

include "common/Assert.ado"
include "common/Debug.ado"
include "common/Version.ado"
include "hdfe/ParseAbsvars.ado"
include "mata/map.mata"
include "hdfe/EstimateDoF.ado"
