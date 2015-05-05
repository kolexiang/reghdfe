mata:
mata set matastrict on

//  Parse absvars and initialize the almost empty MapProblem struct
`Problem' function map_init()
{
	`Integer'		g, G, num_slopes, has_intercept, i
	`Problem' 		S
	`Boolean'		auto_target // Automatically assign target names to all FEs
	`Varname'		basetarget
	`Varlist'		target
	pointer(`FE') 	fe

	S.weightvar = ""
	S.verbose = 0
	S.transform = "cimmino"
	S.acceleration = "conjugate_gradient"
	S.tolerance = 1e-7
	S.maxiterations = 1e4
	S.accel_start = 6
	S.save_ids = 0
	S.groupsize = 10

	// Specific to Aitken:
	S.accel_freq = 3
	S.pause_length = 20
	S.bad_loop_threshold = 1
	S.stuck_threshold = 5e-3
	S.N = .

	S.G = G = st_numscalar("r(G)")
	S.fes = FixedEffect(G)
	S.will_save_fe = 0
	auto_target = st_numscalar("r(savefe)")
	assert(auto_target==1 | auto_target==0)
	if (auto_target) stata(sprintf("cap drop __hdfe*__*"))

	for (g=1; g<=G; g++) {
		fe = &(S.fes[g])
		// recall a->b is the same as (*a).b
		num_slopes = st_numscalar(sprintf("r(num_slopes%f)",g))
		has_intercept = st_numscalar(sprintf("r(has_intercept%1.0f)",g))

		fe->order = g
		fe->num_slopes = num_slopes
		fe->has_intercept = has_intercept
		fe->varlabel = st_global(sprintf("r(varlabel%f)",g))
		fe->ivars = tokens(st_global( sprintf("r(ivars%f)",g) ))
		fe->cvars = tokens(st_global( sprintf("r(cvars%f)",g) ))
		fe->idvarname = sprintf("__ID%f__", g)
		fe->levels = .
		fe->target = J(0,0,"")
		
		basetarget = st_global(sprintf("r(target%f)",g))
		if (basetarget=="" & auto_target) basetarget = sprintf("__hdfe%f__", g)
		if (basetarget!="") {
			S.will_save_fe = 1
			S.save_ids = 1 // We need to save the IDs in order to copy the FEs into the dataset
			target = J(1, has_intercept + num_slopes, basetarget)
			if (has_intercept) stata(sprintf("confirm new variable %s", target[1]))
			for (i=1+has_intercept; i<=length(target); i++) {
				target[i] = target[i] + sprintf("_Slope%f", i-has_intercept)
				stata(sprintf("confirm new variable %s", target[i]))
			}
			fe->target = target
		}
	}
	return(S)
}

void function map_init_weightvar(`Problem' S, `Varname' weightvar) {
	if (weightvar!="") stata(sprintf("confirm numeric variable %s", weightvar))
	S.weightvar = weightvar
}

void function map_init_keepvars(`Problem' S, `Varname' keepvars) {
	if (keepvars!="") stata(sprintf("confirm numeric variable %s, exact", keepvars))
	S.keepvars = tokens(keepvars)
}


void function map_init_save_ids(`Problem' S, `Boolean' save_ids) {
	assert_msg(save_ids==0 | save_ids==1, "save_ids must be 0 or 1")
	S.save_ids = save_ids
}

void function map_init_transform(`Problem' S, `String' transform) {
	transform = strlower(transform)
	// Convert abbreviations
	if (strpos("cimmino", transform)==1) transform = "cimmino"
	if (strpos("kaczmarz", transform)==1) transform = "kaczmarz"
	if (strpos("symmetric_kaczmarz", transform)==1) transform = "symmetric_kaczmarz"
	assert_msg(anyof(("cimmino","kaczmarz","symmetric_kaczmarz"),transform), "invalid transform")
	S.transform = transform
}

void function map_init_verbose(`Problem' S, `Integer' verbose) {
	assert_msg(round(verbose)==verbose, "verbose must be an integer")
	assert_msg(0<=verbose & verbose<=5, "verbose must be between 0 and 5")
	S.verbose = verbose
}

void function map_init_groupsize(`Problem' S, `Integer' groupsize) {
	assert_msg(round(groupsize)==groupsize & groupsize>0, "groupsize must be a positive integer")
	S.groupsize = groupsize
}

void function map_init_acceleration(`Problem' S, `String' acceleration) {
	acceleration = strlower(acceleration)
	// Convert abbreviations
	if (strpos("conjugate_gradient", acceleration)==1 | acceleration=="cg") acceleration = "conjugate_gradient"
	if (strpos("steepest_descent", acceleration)==1 | acceleration=="sd") acceleration = "steepest_descent"
	if (strpos("aitken", acceleration)==1) acceleration = "aitken"
	if (acceleration=="no" | acceleration=="none" | acceleration=="off") acceleration = "none"
	assert_msg(anyof(("conjugate_gradient","steepest_descent", "aitken", "none"),acceleration), "invalid acceleration")
	S.acceleration = acceleration
}

void function map_init_tolerance(`Problem' S, `Real' tolerance) {
	assert_msg(1e-16<=tolerance & tolerance<1, "tolerance must be in the range [1e-16, 1).")
	S.tolerance = tolerance
}

void function map_init_maxiterations(`Problem' S, `Integer' maxiterations) {
	assert_msg(round(maxiterations)==maxiterations, "maxiterations must be an integer")
	assert_msg(maxiterations>0, "maxiterations must be positive")
	S.maxiterations = maxiterations
}

end
