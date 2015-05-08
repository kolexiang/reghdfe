
// Mata code is first, then main hdfe.ado, then auxiliary .ado files
clear mata
include "mata/map.mata"

capture program drop reghdfe
program define reghdfe

* Set Stata version
	version `=clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred

* Intercept version calls
	cap syntax, version
	if !c(rc) {
		Version
		exit
	}

* Intercept replays
	if replay() {
		if (`"`e(cmd)'"'!="reghdfe") error 301
		Replay `0'
		exit
	}

* Finally, call Estimate
	cap noi Estimate `0'
	if (c(rc)) {
		local rc = c(rc)
		cap mata: mata drop HDFE_S // overwrites c(rc)
		exit `rc'
	}
end

// -------------------------------------------------------------------------------------------------
include "common/Assert.ado"
include "common/Debug.ado"
include "common/Version.ado"
// -------------------------------------------------------------------------------------------------
include "internal/Estimate.ado"
	include "internal/Parse.ado"
		include "internal/ParseIV.ado"
		include "internal/ParseVCE.ado"
		include "internal/ParseAbsvars.ado"
		include "internal/ParseDOF.ado"
		include "internal/ParseImplicit.ado"
	
	include "internal/ExpandFactorVariables.ado"
	
	include "internal/Wrapper_regress.ado"
	include "internal/Wrapper_avar.ado"
	include "internal/Wrapper_mwc.ado"
	include "internal/Wrapper_ivregress.ado"
	include "internal/Wrapper_ivreg2.ado"

	include "internal/FixVarnames.ado"
	include "internal/Attach.ado"
	include "internal/Subtitle.ado"
// -------------------------------------------------------------------------------------------------
include "internal/Replay.ado"
	include "internal/Header.ado"
// -------------------------------------------------------------------------------------------------

