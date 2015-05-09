// -------------------------------------------------------------
// Parsing and basic sanity checks for REGHDFE.ado
// -------------------------------------------------------------

cap pr drop Parse
program define Parse

* Remove extra spacing from cmdline (just for aesthetics)
	mata: st_local("cmdline", stritrim(`"reghdfe `0'"') )
	ereturn clear // Clear previous results and drops e(sample)

* Parse the broad syntax (also see map_init(), ParseAbsvars.ado, ParseVCE.ado, etc.)
	syntax anything(id="varlist" name=0 equalok) [if] [in] [fw aw pw/] , ///
	/// Main Options ///
		Absorb(string) ///
		[ ///
		VCE(string) ///
	/// Seldom Used ///
		DOFadjustments(string) ///
		GROUPVAR(name) /// Variable that will contain the first connected group between FEs
	/// Optimization /// Defaults are handled within Mata		
		GROUPsize(string) /// Process variables in batches of #
		TRAnsform(string) ///
		ACCELeration(string) ///
		Verbose(string) ///
		TOLerance(string) ///
		MAXITerations(string) ///
		KEEPSINgletons /// (UNDOCUMENTED) Will keep singletons
		CHECK /// TODO: Implement
		FAST /// TODO: Implement
	/// Regression ///
		ESTimator(string) /// GMM2s CUE LIML
		IVsuite(string) ///
		SAVEFIRST ///
		FIRST ///
		SHOWRAW ///
		VCEUNADJUSTED /// (UNDOCUMENTED) Option when running gmm2s with ivregress; will match results of ivreg2
		SMALL Hascons TSSCONS /// ignored options
		SUBOPTions(string) /// Options to be passed to the estimation command (e.g . to regress)
	/// Multiple regressions in one go ///
		OVER(varname numeric) CLEAR ///
		NESTED /// TODO: Implement
		STAGEs(string) ///
	/// Miscellanea ///
		NOTES(string) /// NOTES(key=value ..)
		] [*] // For display options ; and SUmmarize(stats)

	local allkeys cmdline if in

* Parse varlist: depvar indepvars (endogvars = iv_vars)
	ParseIV `0', estimator(`estimator') ivsuite(`ivsuite') `savefirst' `first' `showraw' `vceunadjusted' `small'
	local keys subcmd model ivsuite estimator depvar indepvars endogvars instruments fe_format ///
		savefirst first showraw vceunadjusted basevars
	foreach key of local keys {
		local `key' "`s(`key')'"
	}
	local allkeys `allkeys' `keys'

* Weights
	if ("`weight'"!="") {
		local weightvar `exp'
		local weighttype `weight'
		local weightexp [`weight'=`weightvar']
		confirm var `weightvar', exact // just allow simple weights

		* Check that weights are correct (e.g. with fweight they need to be integers)
		local num_type = cond("`weight'"=="fweight", "integers", "reals")
		local basenote "weight {res}`weightvar'{txt} can only contain strictly positive `num_type', but"
		qui cou if `weightvar'<0
		Assert (`r(N)'==0), msg("`basenote' `r(N)' negative values were found!")  rc(402)
		qui cou if `weightvar'==0
		if (`r(N)'>0) di as text "`basenote' `r(N)' zero values were found (will be dropped)"
		qui cou if `weightvar'>=.
		if (`r(N)'>0) di as text "`basenote' `r(N)' missing values were found (will be dropped)"
		if ("`weight'"=="fweight") {
			qui cou if mod(`weightvar',1) & `weightvar'<.
			Assert (`r(N)'==0), msg("`basenote' `r(N)' non-integer values were found!") rc(401)
		}
	}
	local allkeys `allkeys' weightvar weighttype weightexp

* Parse Absvars and optimization options
	ParseAbsvars `absorb' // Stores results in r()
		local absorb_keepvars `r(all_ivars)' `r(all_cvars)'
		local N_hdfe `r(G)'

	mata: HDFE_S = map_init() // Reads results from r()
		local will_save_fe = `r(will_save_fe)' // Returned from map_init()
		local original_absvars = "`r(original_absvars)'"
		local extended_absvars = "`r(extended_absvars)'"
		local equation_d = "`r(equation_d)'"
	local allkeys `allkeys' absorb_keepvars N_hdfe will_save_fe original_absvars extended_absvars equation_d

	* Tell Mata what weightvar we have
	if ("`weightvar'"!="") mata: map_init_weights(HDFE_S, "`weightvar'", "`weighttype'")

	* Time/panel variables (need to give them to Mata)
	local panelvar `_dta[_TSpanel]'
	local timevar `_dta[_TStvar]'

	* Parse optimization options (pass them to map_init_*)
	* String options
	local optlist transform acceleration panelvar timevar
	foreach opt of local optlist {
		if ("``opt''"!="") mata: map_init_`opt'(HDFE_S, "``opt''")
	}
	local allkeys `allkeys' `optlist'

	* Numeric options
	local keepsingletons = ("`keepsingletons'"!="")
	local optlist groupsize verbose tolerance maxiterations keepsingletons
	foreach opt of local optlist {
		if ("``opt''"!="") mata: map_init_`opt'(HDFE_S, ``opt'')
	}
	local allkeys `allkeys' `optlist'

	local fast = cond("`fast'"!="", 1, 0) // 1=Yes
	* Fast precludes i) saving FE, ii) running predict, iii) saving groupvar
	if (`fast' & ("`groupvar'"!="" | `will_save_fe'==1)) {
		di as error "(warning: option -fast- not allowed when saving FEs or mobility groups; disabled)"
		local fast 0
	}

	local allkeys `allkeys' fast

* Parse VCE options
	mata: st_local("hascomma", strofreal(strpos("`vce'", ","))) // is there a commma already in `vce'?
	if (!`hascomma') local vce `vce' ,
	ParseVCE `vce' weighttype(`weighttype')
	local keys vceoption vcetype vcesuite vceextra num_clusters clustervars bw kernel dkraay twicerobust
	foreach key of local keys {
		local `key' "`s(`key')'"
	}
	local allkeys `allkeys' `keys'
	
	* Update Mata
	if ("`clustervars'"!="") mata: map_init_clustervars(HDFE_S, "`clustervars'")
	if ("`vceextra'"!="") mata: map_init_vce_is_hac(HDFE_S, 1)

* DoF Adjustments
	if ("`dofadjustments'"=="") local dofadjustments all
	ParseDOF , `dofadjustments'
	local dofadjustments "`s(dofadjustments)'"
	* Mobility groups
	if ("`groupvar'"!="") conf new var `groupvar'
	local allkeys `allkeys' dofadjustments groupvar

* Parse summarize option: [summarize | summarize( stats... [,QUIetly])]
	* Note: ParseImplicit deals with "implicit" options and fills their default values
	local default_stats mean min max
	ParseImplicit, opt(SUmmarize) default(`default_stats') input(`options') syntax([namelist(name=stats)] , [QUIetly]) inject(stats quietly)
	local summarize_quietly = ("`quietly'"!="")
	if ("`stats'"=="" & "`quietly'"!="") local stats `default_stats'
	local allkeys `allkeys' stats summarize_quietly

* Parse over() option
	local clear = "`clear'"!=""
 	if ("`over'"!="") {
		unab over : `over', max(1)
		Assert (`clear'), msg("over() requires the -clear- option")
	}
	local allkeys `allkeys' over clear

* Nested
	local nested = cond("`nested'"!="", 1, 0) // 1=Yes
	if (`nested' & !("`model'"=="ols" & "`vcetype'"=="unadjusted") ) {
		Debug, level(0) msg("(option nested ignored, only works with OLS and conventional/unadjusted VCE)") color("error")
	}
	local allkeys `allkeys' nested

* Stages
	assert "`model'"!="" // just to be sure this goes after `model' is set
	local iv_stage iv
	local stages : list stages - iv_stage
	local valid_stages ols first acid reduced
	local wrong_stages : list stages - valid_stages
	Assert "`wrong_stages'"=="", msg("Error, invalid stages(): `wrong_stages'")
	if ("`stages'"!="") {
		Assert "`model'"=="iv", msg("Error, stages() only valid with an IV regression")
		local stages `stages' `iv_stage' // Put -iv- *last* (so it does the -restore-; note that we don't need it first to trim MVs b/c that's done earlier)
	}
	else {
		local stages none // So we can loop over stages
	}
	local allkeys `allkeys' stages

* Parse Coef Table Options (do this last!)
	_get_diopts diopts options, `options' // store in `diopts', and the rest back to `options'
	Assert `"`options'"'=="", msg(`"invalid options: `options'"')
	if ("`hascons'`tsscons'"!="") di in ye "(option `hascons'`tsscons' ignored)"
	local allkeys `allkeys' diopts


* Other keys:
	local allkeys `allkeys' suboptions notes
	// Missing keys: check

* Return values
	Debug, level(3) newline
	Debug, level(3) msg("{title:Parsed options:}")
	foreach key of local allkeys {
		if (`"``key''"'!="") Debug, level(3) msg("  `key' = " as result `"``key''"')
		c_local `key' `"``key''"' // Inject values into caller (reghdfe.ado)
	}
end

