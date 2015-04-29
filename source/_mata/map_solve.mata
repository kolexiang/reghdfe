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
		resid = (*T)(S, y)
		if (check_convergence(S, iter, resid, y)) break
		y = resid
	}
	return(resid)
}

// -------------------------------------------------------------------------------------------------

`Group' function accelerate_cg(`Problem' S, `Group' y, `FunctionPointer' T) {
	`Integer'	iter
	`Group'		resid, r, u, v, resid_old
	real rowvector alpha, beta, ssr, ssr_old, denominator
	`Matrix' R

timer_on(61)
	resid = y
	r = y - (*T)(S, resid)
	ssr = quadcolsum(r :* r) // cross(r,r) when cols(y)==1 // BUGBUG maybe diag(quadcross()) is faster?
	u = r
timer_off(61)
	for (iter=1; iter<=S.maxiterations; iter++) {
timer_on(62)
		v = u - (*T)(S, u)
		alpha = safe_divide( ssr , quadcolsum(u :* v) )
		resid = resid - alpha :* u
		ssr_old = ssr
		r = r - alpha :* v
		ssr = quadcolsum(r :* r)
		beta = safe_divide( ssr , ssr_old) // Fletcher-Reeves formula, but it shouldn't matter in our problem
		u = r + beta :* u
		if (check_convergence(S, iter, ssr)) break
timer_off(62)
	}
	return(resid)
}

// -------------------------------------------------------------------------------------------------

`Group' function accelerate_sd(`Problem' S, `Group' y, `FunctionPointer' T) {
	`Integer'	iter
	`Group' proj
	real rowvector t, denominator

	for (iter=1; iter<=S.maxiterations; iter++) {
		proj = y - (*T)(S, y)
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
		
		resid = (*T)(S, y)
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

`Boolean' check_convergence(`Problem' S, `Integer' iter, `Group' y_new,| `Group' y_old) {
	`Boolean'	done, is_last_iter
	`Real'		update_error

	// max() ensures that the result when bunching vars is at least as good as when not bunching
	if (args()>=4) {
		update_error = max(mean(reldif(y_new, y_old)))
	}
	else {
		// We don't have two vectors, instead one SSR
		assert(rows(y_new)==1)
		//update_error = max(sqrt(y_new))
		update_error = max(y_new)
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

`Group' function transform_cimmino(`Problem' S, `Group' y) {
	`Integer' 	N, Q, g, G
	`Group'		ans
	N = rows(y)
	Q = cols(y)
	G = S.G

	ans = map_project(S, 1, y)
	for (g=2; g<=G; g++) {
		ans = ans + map_project(S, g, y)
	}
	return(y - ans / G)
}

// -------------------------------------------------------------------------------------------------

`Group' function transform_kaczmarz(`Problem' S, `Group' y) {
	`Integer' 	N, Q, g, G
	`Group'		ans
	N = rows(y)
	Q = cols(y)
	G = S.G
	ans = y - map_project(S, 1, y)
	for (g=2; g<=G; g++) {
		ans = ans - map_project(S, g, ans)
	}
	return(ans)
}

// -------------------------------------------------------------------------------------------------

 `Group' function transform_sym_kaczmarz(`Problem' S, `Group' y) {
	`Integer' 	N, Q, g, G
	`Group'		ans
	timer_on(71)
	N = rows(y)
	Q = cols(y)
	G = S.G
	ans = y - map_project_fast(S, 1, y)
	for (g=2; g<=G; g++) {
		ans = ans - map_project_fast(S, g, ans)
	}
	for (g=G-1; g>=1; g--) {
		ans = ans - map_project_fast(S, g, ans)
	}
	timer_off(71)
	return(ans)
}


end
