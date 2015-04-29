mata:
mata set matastrict on

void function map_solve(`Problem' S, `Varlist' vars, | `Varlist' newvars) {
	`Integer' i, Q
	`Group' y
	`FunctionPointer' transform
	real rowvector stdevs

	// Load data
	vars = tokens(vars)
	y = st_data(., vars)
	st_dropvar(vars) // We need the new ones on double precision
	if (args()<3 | newvars=="") newvars = vars
	assert_msg(length(vars)==length(newvars), "mapsolve error: newvars must have the same size as vars")
	assert_msg(S.N==st_nobs(), "dataset cannot change after map_precompute()")

	// Standardize the vars
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

	if (S.verbose>1) printf("{txt} - Saving output\n")
	st_store(., st_addvar("double", newvars), y :* stdevs)

	// We could do without regress in the simple case:
	//y = y :* stdevs
	//qrsolve(y[., 2..cols(y)] , y[., 1])
}

// -------------------------------------------------------------------------------------------------
// Acceleration Schemes
// -------------------------------------------------------------------------------------------------

`Group' function accelerate_none(`Problem' S, `Group' y, `FunctionPointer' T) {
	`Integer'	iter
	`Group'		resid

	for (iter=1; iter<=S.maxiterations; iter++) {
		(*T)(S, y, resid)
		if (check_convergence(S, iter, resid, y)) break
		y = resid
	}
	return(resid)
}

// -------------------------------------------------------------------------------------------------
// Memory cost is approx = 4*size(y) (actually 3 since y is already there)
// But we need to add maybe 1 more due to u:*v
// And I also need to check how much does project and T use..
// Double check with a call to memory

// For discussion on the stopping criteria, see the following presentation:
// Arioli & Gratton, "Least-squares problems, normal equations, and stopping criteria for the conjugate gradient method". URL: https://www.stfc.ac.uk/SCD/resources/talks/Arioli-NAday2008.pdf
// Basically, we will use the Hestenes and Siefel rule

`Group' function accelerate_cg(`Problem' S, `Group' y, `FunctionPointer' T) {
	// BUGBUG iterate the first 6? without acceleration??
	`Integer'	iter, d, Q
	`Group'		r, u, v
	real rowvector alpha, beta, ssr, ssr_old, denominator, improvement_potential
	`Matrix' R, recent_ssr
	Q = cols(y)
	d = 1 // BUGBUG Set it to 2/3 // Number of recent SSR values to use for convergence criteria (lower=faster & riskier)
	improvement_potential = quadcolsum(y :* y)
	recent_ssr = J(d, Q, .)
	
	(*T)(S, y, r, 1)
	ssr = quadcolsum(r :* r) // cross(r,r) when cols(y)==1 // BUGBUG maybe diag(quadcross()) is faster?
	u = r

	for (iter=1; iter<=S.maxiterations; iter++) {
		(*T)(S, u, v, 1) // This is the hotest loop in the entire program
		// BUGBUG: colsum or quadcolsum??
		alpha = safe_divide( ssr , quadcolsum(u :* v) )
		recent_ssr[1 + mod(iter-1, d), .] = alpha :* ssr
		improvement_potential = improvement_potential - alpha :* ssr
		y = y - alpha :* u
		r = r - alpha :* v
		ssr_old = ssr
		ssr = quadcolsum(r :* r)
		beta = safe_divide( ssr , ssr_old) // Fletcher-Reeves formula, but it shouldn't matter in our problem
		u = r + beta :* u
		// Convergence if sum(recent_ssr) > tol^2 * improvement_potential
		if ( check_convergence(S, iter, sum(recent_ssr), improvement_potential, "hestenes") ) break
	}
	return(y)
}

// -------------------------------------------------------------------------------------------------

`Group' function accelerate_sd(`Problem' S, `Group' y, `FunctionPointer' T) {
	`Integer'	iter
	`Group' proj
	real rowvector t, denominator

	for (iter=1; iter<=S.maxiterations; iter++) {
		(*T)(S, y, proj, 1)
		t = safe_divide( quadcolsum(y :* proj) , quadcolsum(proj :* proj) )
		if (check_convergence(S, iter, y-proj, y)) break
		y = y - t :* proj
	}
	return(y-proj)
}

// -------------------------------------------------------------------------------------------------
// This is method 3 of Macleod (1986), a vector generalization of the Aitken-Steffensen method
// Also: "when numerically computing the sequence.. stop..  when rounding errors become too 
// important in the denominator, where the ^2 operation may cancel too many significant digits"
// Note: Sometimes the iteration gets "stuck"; can we unstuck it with adding randomness
// in the accelerate decision? There should be a better way.. (maybe symmetric kacz instead of standard one?)

`Group' function accelerate_aitken(`Problem' S, `Group' y, `FunctionPointer' T) {
	`Integer'	iter
	`Group'		resid, y_old, delta_sq
	`Boolean'	accelerate
	real rowvector t

	//S.pause_length = 20
	//S.bad_loop_threshold = 1
	//S.stuck_threshold = 5e-3
	// old_error = oldest_error = bad_loop = acceleration_countdown = 0

	for (iter=1; iter<=S.maxiterations; iter++) {
		
		(*T)(S, y, resid)
		accelerate = iter>=S.accel_start & !mod(iter,S.accel_freq)

		// Accelerate
		if (accelerate) {
			delta_sq = resid - 2 * y + y_old // = (resid - y) - (y - y_old) // Equivalent to D2.resid
			// t is just (d'd2) / (d2'd2)
			t = safe_divide( quadcolsum( (resid - y) :* delta_sq) ,  quadcolsum(delta_sq :* delta_sq) )
			resid = resid - t :*  (resid - y)
		}

		// Only check converge on non-accelerated iterations
		// BUGBUG: Do we need to disable the check when accelerating?
		if (check_convergence(S, iter, accelerate? resid :* .  : resid, y)) break
		
		// Experimental: Pause acceleration
		//if (accelerate) {
		//	improvement = max(( (old_error-update_error)/update_error , (oldest_error-update_error)/update_error ))
		//	bad_loop = improvement < stuck_threshold ? bad_loop+1 : 0
		//	// bad_loop, improvement, update_error, old_error, oldest_error
		//	// Tolerate two problems (i.e. 6+2=8 iters) and then try to unstuck
		//	if (bad_loop>bad_loop_threshold) {
		//		bad_loop = 0
		//		if (VERBOSE==3) printf(" Fixed point iteration seems stuck, acceleration paused\n")
		//		acceleration_countdown = pause_length
		//	}
		//	assert(bad_loop<=3)	
		//	oldest_error = old_error
		//	old_error = update_error
		//}
		//
		y_old = y // y_old is resid[iter-2]
		y = resid // y is resid[iter-1]
	}
	return(resid)
}

// -------------------------------------------------------------------------------------------------

`Boolean' check_convergence(`Problem' S, `Integer' iter, `Group' y_new, `Group' y_old,| `String' method) {
	`Boolean'	done, is_last_iter
	`Real'		update_error

	// max() ensures that the result when bunching vars is at least as good as when not bunching
	if (args()<5) method = "vectors" 

	if (method=="vectors") {
		update_error = max(mean(reldif(y_new, y_old)))
	}
	else if (method=="hestenes") {
		update_error = max(safe_divide( sqrt(y_new) , sqrt(y_old) ))
	}
	else {
		exit(error(100))
	}

	done = update_error <= S.tolerance
	is_last_iter = iter==S.maxiterations
	
	if (done) {
		if (S.verbose>1) printf("\n{txt} - Converged in %g iterations (last error =%3.1e)\n", iter, update_error)
		if (S.verbose==1) printf("{txt} converged in %g iterations last error =%3.1e)\n", iter, update_error)
	}
	else if (is_last_iter) {
		printf("\n{err}convergence not achieved in %g iterations (last error=%e); try increasing maxiter() or decreasing tol().\n", S.maxiterations, update_error)
		exit(430)
	}
	else {
		if ((S.verbose>=2 & S.verbose<=3 & mod(iter,1)==0) | (S.verbose==1 & mod(iter,100)==0)) {
			printf("{txt}.")
			displayflush()
		}
		if (S.verbose>=2 & S.verbose<=3 & mod(iter,100)==0) printf("{txt}%9.1f\n", update_error/S.tolerance)
	}
	return(done)
}

// -------------------------------------------------------------------------------------------------

real rowvector safe_divide(real rowvector numerator, real rowvector denominator) {
	 // If the denominator goes below machine precision, the division explodes
	return( numerator :/ colmax(denominator \ J(1,cols(denominator),epsilon(1))) )
}


// -------------------------------------------------------------------------------------------------
// Transformations: Compute RESIDUALS, not projections
// -------------------------------------------------------------------------------------------------

void function transform_cimmino(`Problem' S, `Group' y, `Group' ans,| `Boolean' get_proj) {
	`Integer' 	N, Q, g, G
	N = rows(y)
	Q = cols(y)
	G = S.G
	if (args()<4) get_proj = 0

	ans = map_project(S, 1, y)
	for (g=2; g<=G; g++) {
		ans = ans + map_project(S, g, y)
	}
	ans = get_proj ? ans / G : y - ans / G
}

// -------------------------------------------------------------------------------------------------

void function transform_kaczmarz(`Problem' S, `Group' y, `Group' ans,| `Boolean' get_proj) {
	`Integer' 	N, Q, g, G
	N = rows(y)
	Q = cols(y)
	G = S.G
	if (args()<4) get_proj = 0

	ans = y - map_project(S, 1, y)
	for (g=2; g<=G; g++) {
		ans = ans - map_project(S, g, ans)
	}
	if (get_proj) ans = y - ans
}

// -------------------------------------------------------------------------------------------------

 void function transform_sym_kaczmarz(`Problem' S, `Group' y, `Group' ans,| `Boolean' get_proj) {
	`Integer' 	N, Q, g, G
	// BUGBUG: Streamline and remove all those "ans - .." lines?
	N = rows(y)
	Q = cols(y)
	G = S.G
	if (args()<4) get_proj = 0

	ans = y - map_project(S, 1, y)
	for (g=2; g<=G; g++) {
		ans = ans - map_project(S, g, ans)
	}
	for (g=G-1; g>=1; g--) {
		ans = ans - map_project(S, g, ans)
	}
	if (get_proj) ans = y - ans
}


end
