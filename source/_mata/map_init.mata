mata:
mata set matastrict on

//  Parse absvars and initialize the almost empty MapProblem struct
`Problem' function map_init(| `Varname' weightvar)
{
	`Integer'		g, G
	`Problem' 		S
	pointer(`FE') 	fe

	S.weightvar = args()<1 ? "" : weightvar
	S.verbose = 0
	S.transform = "cimmino"
	S.acceleration = "conjugate_gradient"
	S.tolerance = 1e-7
	S.maxiterations = 1e4
	S.accel_start = 6

	// Specific to Aitken:
	S.accel_freq = 3
	S.pause_length = 20
	S.bad_loop_threshold = 1
	S.stuck_threshold = 5e-3

	S.G = G = st_numscalar("r(G)")
	S.fes = FixedEffect(G)
	for (g=1; g<=G; g++) {
		fe = &(S.fes[g])
		// recall a->b is the same as (*a).b
		fe->order = g
		fe->varlabel = st_global(sprintf("r(varlabel%f)",g))
		fe->num_slopes = st_numscalar(sprintf("r(num_slopes%f)",g))
		fe->has_intercept = st_numscalar(sprintf("r(has_intercept%1.0f)",g))
		fe->ivars = tokens(st_global( sprintf("r(ivars%f)",g) ))
		fe->cvars = tokens(st_global( sprintf("r(cvars%f)",g) ))
		fe->target = st_global(sprintf("r(target%f)",g))
		fe->idvarname = sprintf("__ID%f__", g)
		fe->levels = .
	}
	return(S)
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
