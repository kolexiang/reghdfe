
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

include "common/Assert.ado"
include "common/Debug.ado"
include "common/Version.ado"
include "hdfe/ParseAbsvars.ado"


// ni siquiera necesito el UID en la base despues del restore! bajarmelo nomas
// es mas, CREAR EL UID JUSTO DESPUES DE PRESERVE!!!!!!!!!!!!!!!!!!!!!
// luego del restore anhado data asi:
st_store(uid, st_addvar("float", "x"), x)







* UPDATE INCLUDES
* Note: Assert and Debug must go first
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
