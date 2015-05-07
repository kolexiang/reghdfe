// -------------------------------------------------------------------------------------------------
// Parse optimization options and save them into the new MAP Mata Structure
// -------------------------------------------------------------------------------------------------
capture program drop ParseMAP
program define ParseMAP, sclass
	// Note: All the option checks are done in the Mata code
	syntax, ///
		Absorb(string) ///
		[ ///
		WEIGHTVAR(string) ///
		GROUPsize(string) /// Process variables in groups of #
		TRANSFORM(string) ///
		ACCELeration(string) ///
		Verbose(string) ///
		TOLerance(string) ///
		MAXITerations(string) ///
		CLUSTERVARs(string) ///
		] [*]

	ParseAbsvars `absorb'
	mata: HDFE_S = map_init()

	local optlist groupsize weightvar transform acceleration verbose tolerance maxiterations
	foreach opt of local optlist {
		if ("`opt'"!="") mata: map_init_`opt'(S, "`opt'")
	}

	* Return remaining options
	sreturn local options `options'
end
