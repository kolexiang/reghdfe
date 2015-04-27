mata:
mata set matastrict on
`Group' function mapsolve_project(`Problem' S, `Integer' g, `Group' y) {
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
	//proj = has_intercept ? sum(y) / denom : 0 // ymean part
	//if K>0 {
	//	b = inv_xx[] * crossdev(x_demeaned, 0, w, y, proj)
	//	proj = proj + x_demeaned * b
	//}

		i_upper = S.fes[g].offsets[j]
		tmp_count = S.fes[g].counts[j]
		tmp_w = has_weights ? (*p_sorted_w)[| i_lower \ i_upper |] : 1
		tmp_y = (*p_sorted_y)[| i_lower , 1 \ i_upper , . |]
		ymean = has_intercept ? (quadcolsum(tmp_y :* tmp_w) / tmp_count) : 0

		if (K>0) {
			tmp_x = S.fes[g].x[| i_lower , 1 \ i_upper , . |]
			b = S.fes[g].inv_xx[| 1+(j-1)*K , 1 \ j*K , . |] * crossdev(tmp_x, zero, tmp_w, tmp_y, ymean)
			// betas[j, .] = b'
		}
		
		//ans[| i_lower , 1 \ i_upper , . |] = K>0 ? (ymean :+ tmp_x*b) : (ymean :+ ans[| i_lower , 1 \ i_upper , . |])
		ans[| i_lower , 1 \ i_upper , . |] = K>0 ? (ymean :+ tmp_x*b) : (ymean :+ J(i_upper-i_lower+1,Q,0))
		i_lower = i_upper + 1
	}

	//predicted = rowsum( (x , J(rows(predicted),1,1)) :* betas[group,.] )
	return(ans[S.fes[g].inv_p, .])
}
end


/*
Projection Strategy

General case (multiple slopes, weights):

	On Setup:
		Divide Xs (cvars) by their STDEV (not by groups, all in one go)
		
		Compute x_demeaned (BY GROUP) = x - sum(x:* w) / sum(w) (o count!)
		Compute inv_xx (BY GROUP), a KL*K matrix (K excludes constant)
			= invsym(cross(x_demeaned, w, x_demeaned))

		[Why do we compute x_demeaned? b/c we don't use x at all, only the demeaned so it's easy even for computing proj]

	When receiving Y:
		Standardize Y like we did the Xs

	On the fly:
		ymean = sum(y) / denom (where denom = N or = sum(w))
		b = inv_xx * crossdev(x_demeaned, 0, w, y, ymean)
		proj = ymean + x_demeaned * b

No weights:
	
	Easy-ish, just change all the -w- part (can we put 1? or drop it?)

No intercept, just slopes:

	In that case, we don't need the demeaning (which is just FW wrt to the constant)
	x_demeaned is just x; inv_xx stays the same; we don't compute ymean; b stays the same except we use cross(); proj has no ymean part.

Intercept and just one slope:

	The only potential optimization is that invsym can be replaced by 1 / cross(..)

Only intercept (MOST COMMON CASE):

	proj is just ymean
	we have already preocomputed denom in the same way as in prev. cases
	we just need to compute sum(y) , and we use the offset vector for that

Maybe we can do just ONE function that is general..

Eg:
proj = has_intercept ? sum(y) / denom : 0 // ymean part
if K>0 {
	b = inv_xx[] * crossdev(x_demeaned, 0, w, y, proj)
	proj = proj + x_demeaned * b
}
store back in proj vector

*/
