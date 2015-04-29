// -------------------------------------------------------------------------------------------------
// Container for the Within-Transformation Code
// -------------------------------------------------------------------------------------------------

* Type Aliases
	local Boolean 		real scalar
	local Integer 		real scalar
	local Real 			real scalar
	local Vector		real colvector
	local Matrix		real matrix
	local Series		real colvector // Should be N*1
	local Group			real matrix // Should be N*K
	local String 		string scalar
	local Varname 		string scalar
	local Varlist 		string rowvector // rowvector so they match tokens()
	local Problem		struct MapProblem scalar
	local FE			struct FixedEffect scalar
	local FunctionPointer pointer(`Group' function) scalar // Used for the Accelerate & Transform fns

* Mata Includes
	include assert_msg.mata
	include FixedEffect.mata
	include MapProblem.mata
	include map_init.mata
	include map_precompute.mata
	include map_precompute_part1.mata
	include map_precompute_part2.mata
	include map_project.mata
	include map_solve.mata
// -------------------------------------------------------------------------------------------------
