*! reghdfe VERSION_NUMBER
*! Sergio Correia (sergio.correia@duke.edu)
* (built from multiple source files using build.py)

include _mata/reghdfe.mata

cap pr drop reghdfe
program define reghdfe
	local version `=clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred
	qui version `version'

	* Intercept version calls
	cap syntax, version
	local rc = _rc
	 if (`rc'==0) {
		Version
		exit
	}

	* Intercept multiprocessor/parallel calls
	cap syntax, instance [*]
	local rc = _rc
	 if (`rc'==0) {
		ParallelInstance, `options'
		exit
	}

	if replay() {
		if (`"`e(cmd)'"'!="reghdfe") error 301
		Replay `0'
	}
	else {
		* Estimate, and then clean up Mata in case of failure
		mata: st_global("reghdfe_pwd",pwd())
		Stop // clean leftovers for a possible [break]
		cap noi Estimate `0'
		if (_rc) {
			local rc = _rc
			Stop
			exit `rc'
		}
	}
end

* Note: Assert and Debug must go first
include "common/Assert.ado"
include "common/Debug.ado"
include "common/Version.ado"
include "common/SortPreserve.ado"

include "mata/fix_psd.mata"

include "reghdfe/Estimate.ado"
include "reghdfe/Parse.ado"
include "hdfe/DropSingletons.ado"
include "reghdfe/ExpandFactorVariables.ado"
include "reghdfe/FixVarnames.ado"
include "reghdfe/Wrapper_regress.ado"
include "reghdfe/Wrapper_mwc.ado"
include "reghdfe/Wrapper_avar.ado"
include "reghdfe/Wrapper_ivregress.ado"
include "reghdfe/Wrapper_ivreg2.ado"
include "reghdfe/Attach.ado"
include "reghdfe/Replay.ado"
include "reghdfe/Header.ado"

include "hdfe/ConnectedGroups.ado"
include "hdfe/GenerateID.ado"
include "hdfe/AverageOthers.ado"
include "hdfe/EstimateDoF.ado"
include "hdfe/Start.ado"
include "hdfe/ParseOneAbsvar.ado"
include "hdfe/Precompute.ado"
include "hdfe/Demean.ado"
include "hdfe/DemeanParallel.ado"
include "hdfe/ParallelInstance.ado"
include "hdfe/Save.ado"
include "hdfe/Stop.ado"
