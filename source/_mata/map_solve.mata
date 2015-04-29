mata:
mata set matastrict on

void function map_solve(`Problem' S, `Varlist' vars, | `Varlist' newvars) {
	`Integer' i, Q
	`Group' y
	`FunctionPointer' transform
	real rowvector stdevs

	if (S.verbose>0) printf("{txt}mata: map_solve()\n")

	// Load data
	if (S.verbose>0) printf("{txt} - Loading variables into Mata\n")
	vars = tokens(vars)
	y = st_data(., vars)
	st_dropvar(vars) // We need the new ones on double precision
	if (args()<3 | newvars=="") newvars = vars
	assert_msg(length(vars)==length(newvars), "mapsolve error: newvars must have the same size as vars")
	assert_msg(S.N==st_nobs(), "dataset cannot change after map_precompute()")

	// Standardize the vars
	if (S.verbose>0) printf("{txt} - Standardizing variables\n")
	Q = cols(y)
	stdevs = J(1,Q,.)
	for (i=1; i<=Q; i++) {
		stdevs[i] = max((  sqrt(quadvariance(y[., i])) , sqrt(epsilon(1)) ))
	}
	y = y :/ stdevs

	// TODO: Report PEAK MEMORY that will be used
	// EG: de por si uso 2*size(y)
	// Luego bajo y asi que me quedo con y + ..
	// Contar cuantos vectores creo en los aceleradores y en los proyectores

	if (S.verbose>0) printf("{txt} - Solving problem (acceleration={res}%s{txt}, transform={res}%s{txt})\n", S.acceleration, S.transform)

	// Warnings
	if (S.transform=="kaczmarz" & S.acceleration=="conjugate_gradient") {
		printf("{err}(warning: convergence is {bf:unlikely} with transform=kaczmarz and accel=CG)\n")
	}

	// Load transform pointer
	if (S.transform=="cimmino") transform = &transform_cimmino()
	if (S.transform=="kaczmarz") transform = &transform_kaczmarz()
	if (S.transform=="symmetric_kaczmarz") transform = &transform_sym_kaczmarz()

	// Call acceleration routine
	if (S.acceleration=="none") y = accelerate_none(S, y, transform)
	if (S.acceleration=="conjugate_gradient") y = accelerate_cg(S, y, transform)
	if (S.acceleration=="steepest_descent") y = accelerate_sd(S, y, transform)
	if (S.acceleration=="aitken") y = accelerate_aitken(S, y, transform)

	if (S.verbose>1) printf("{txt} - Saving transformed variables\n")
	st_store(., st_addvar("double", newvars), y :* stdevs)

	// We could do without regress in the simple case:
	//y = y :* stdevs
	//qrsolve(y[., 2..cols(y)] , y[., 1])
}

// -------------------------------------------------------------------------------------------------

real rowvector safe_divide(real rowvector numerator, real rowvector denominator) {
	 // If the denominator goes below machine precision, the division explodes
	return( numerator :/ colmax(denominator \ J(1,cols(denominator),epsilon(1))) )
}

end
