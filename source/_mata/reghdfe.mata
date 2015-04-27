// -------------------------------------------------------------------------------------------------
// Container for the Within-Transformation Code
// -------------------------------------------------------------------------------------------------

* Type Aliases
	local Boolean 		real scalar
	local Integer 		real scalar
	local Vector		real colvector
	local Matrix		real matrix
	local Series		real colvector // Should be N*1
	local Group			real matrix // Should be N*K
	local Varlist 		string rowvector // rowvector so they match tokens()
	local Varname 		string scalar
	local Problem		struct MapProblem scalar
	local FE			struct FixedEffect scalar

* Mata Includes
	include assert_msg.mata
	include FixedEffect.mata
	include MapProblem.mata
	include mapsolve_init.mata
	include mapsolve_precompute.mata
	include mapsolve_precompute_part1.mata
	include mapsolve_precompute_part2.mata
	include mapsolve_project.mata
// -------------------------------------------------------------------------------------------------
