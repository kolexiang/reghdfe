mata:
mata set matastrict on
`Group' function map_project(`Problem' S, `Integer' g, `Group' y) {
	`Integer' 	K, L, N, Q // Q is the number of depvars
	`Integer' 	j, i_lower, i_upper // j loops over levels, i loops over observations
	`Boolean' 	has_weights, sortedby, has_intercept
	`Series'	sorted_y, sorted_w
	`Group'		ans
	`Vector'	b, tmp_w, tmp_count
	real rowvector ymean // 1*Q
	real rowvector zero // 1*K
	`Matrix'	tmp_y, tmp_x
	pointer(`Series') scalar p_sorted_y, p_sorted_w

	// PROFILE TO SEE IF THIS HELPS OR NOT AT ALL
	//pointer(`Vector') scalar p_offset
	//p_offset = &(S.fes[g].offsets)

	timer_on(81)

	has_weights = S.weightvar !=""
	sortedby = S.fes[g].is_sortedby
	has_intercept = S.fes[g].has_intercept
	K = S.fes[g].num_slopes
	N = rows(y)
	Q = cols(y)
	L = S.fes[g].levels

	// Minimize copy+order operations on y
	if (!sortedby) {
		sorted_y = y[S.fes[g].p, .]
		if (has_weights) sorted_w = S.w[S.fes[g].p, .]
	}
	p_sorted_y = sortedby? &y : &sorted_y
	if (has_weights) p_sorted_w = sortedby? &(S.w) : &sorted_w
	if (K>0) zero = J(1,K,0)

	ans = J(N, Q, 0)
	i_lower = 1
	for (j=1; j<=L; j++) {
		i_upper = S.fes[g].offsets[j]
		tmp_count = S.fes[g].counts[j]
		
		if (has_weights) tmp_w = (*p_sorted_w)[| i_lower \ i_upper |]
		tmp_y = (*p_sorted_y)[| i_lower , 1 \ i_upper , . |]
		if (has_weights) {
			ymean = has_intercept ? (quadcolsum(tmp_y :* tmp_w) / tmp_count) : 0
		}
		else {
			ymean = has_intercept ? (quadcolsum(tmp_y) / tmp_count) : 0
		}

		if (K>0) {
			tmp_x = S.fes[g].x[| i_lower , 1 \ i_upper , . |]
			b = S.fes[g].inv_xx[| 1+(j-1)*K , 1 \ j*K , . |] * crossdev(tmp_x, zero, tmp_w, tmp_y, ymean)
			// betas[j, .] = b'
		}
		
		//ans[| i_lower , 1 \ i_upper , . |] = K>0 ? (ymean :+ tmp_x*b) : (ymean :+ ans[| i_lower , 1 \ i_upper , . |])
		ans[| i_lower , 1 \ i_upper , . |] = K>0 ? (ymean :+ tmp_x*b) : (ymean :+ J(i_upper-i_lower+1,Q,0))
		i_lower = i_upper + 1
	}

	
	ans = sortedby ? ans : ans[S.fes[g].inv_p, .]
	timer_off(81)
	//return(sortedby ? ans : ans[S.fes[g].inv_p, .])
	return(ans)
}


`Group' function map_project_fast(`Problem' S, `Integer' g, `Group' y) {
	`Integer' 	K, L, N, Q // Q is the number of depvars
	`Integer' 	j, i_lower, i_upper // j loops over levels, i loops over observations
	`Boolean' 	sortedby
	`Series'	sorted_y
	`Group'		ans
	`Vector'	tmp_count, asd
	real rowvector ymean // 1*Q
	`Matrix'	tmp_y
	pointer(`Series') scalar p_sorted_y

	timer_on(81)

	sortedby = S.fes[g].is_sortedby
	N = rows(y)
	Q = cols(y)
	L = S.fes[g].levels

	// Minimize copy+order operations on y
timer_on(92)
	ans = sortedby ? y : y[S.fes[g].p, .] // 50s / 140s
	i_lower = 1
timer_off(92)
asd = S.fes[g].counts[j]

	for (j=1; j<=L; j++) {
		timer_on(93)
		i_upper = S.fes[g].offsets[j]
		timer_off(93)
		timer_on(94)
		tmp_count = asd[j]
		timer_off(94)
		timer_on(95)
		tmp_y = ans[| i_lower , 1 \ i_upper , . |] // 9s
		timer_off(95)
		timer_on(96)
		ymean =  (	 colsum(tmp_y) / tmp_count) // 8s
		timer_off(96)
		timer_on(97)
		ans[| i_lower , 1 \ i_upper , . |] = (ymean :+ J(i_upper-i_lower+1,Q,0)) // 30s
		timer_off(97)
		timer_on(98)
		i_lower = i_upper + 1
		timer_off(98)
	}

	timer_on(99)
	ans = sortedby ? ans : ans[S.fes[g].inv_p, .] // 44s
	timer_off(99)
	timer_off(81)
	//return(sortedby ? ans : ans[S.fes[g].inv_p, .])
	return(ans)
}

end

