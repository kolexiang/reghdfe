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
	syntax varlist [if] [in] [fw aw pw/] , ///
	/// Main Options ///
		Absorb(string) ///
		[ ///
		PARTIAL(varlist numeric) /// Additional regressors besides those in absorb()
		SAMPLE(name) ///
		Generate(name) CLEAR /// Replace dataset, or just add new variables
		CLUSTERVARs(string) /// Used to estimate the DoF
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
	if ("`options'"!="") di as error "unused options: `options'"
	if ("`sample'"!="") confirm new variable `sample'
	Assert ("`generate'"!="") + ("`clear'"!="") == 1 , ///
		msg("hdfe error: you need to specify one and only one of the following options: clear generate(...)")
	Assert "`: list varlist & partial'"=="", ///
		msg("variables in varlist cannot appear in partial()")
	if ("`weight'"!="") {
		local weightvar `exp'
		local weighttype `weight'
		local weightequal =
		confirm var `weightvar' // just allow simple weights
	}

* From now on, we will pollute the Mata workspace, so wrap this in case of error
	cap noi {

* Create Mata structure
	ParseAbsvars `absorb' // Stores results in r()
	mata: HDFE_S = map_init() // Reads results from r()
	local optlist groupsize weightvar transform acceleration verbose tolerance maxiterations
	foreach opt of local optlist {
		if ("``opt''"!="") di as error `"mata: map_init_`opt'(HDFE_S, "``opt''")"'
		if ("``opt''"!="") mata: map_init_`opt'(HDFE_S, "``opt''")
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
	mata: map_precompute(HDFE_S)

* Within Transformation
	mata: map_solve(HDFE_S, "`varlist'", "`newvars'", "`partial'")

* Cleanup after an error
	} // cap noi
	if c(rc) {
		local rc = c(rc)
		cap mata: mata drop HDFE_S // overwrites c(rc)
		exit `rc'
	}
end

include "common/Assert.ado"
include "mata/map.mata"
include "hdfe/ParseAbsvars.ado"
