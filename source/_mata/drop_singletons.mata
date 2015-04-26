mata:
mata set matastrict on

// Note: So far, this function uses 5 + G rows, which may be quite a lot!

void function drop_singletons(struct FixedEffect vector fes, 
    	`Integer' verbose, | `Varlist' keepvars) {

	`Integer' G, i, j, n, g, h, i_last_singleton, num_singletons, sortedby
	`Series' singleton, sum_singleton, id, inv_p
	`Varlist' idvarnames
	`Varname' ivar
	string scalar vartype
	pointer(`Series') scalar pp // Just to shorten code
	transmorphic counter

	G = length(fes)

	// Initialize and fill counter
	counter = asarray_create()
	asarray_notfound(counter, 0)
	for (g=1; g<=G; g++) {
		keepvars = keepvars , fes[g].ivars
	}
	for (i=1; i<=length(keepvars); i++) {
		asarray(counter, keepvars[i], asarray(counter, keepvars[i])+1)
	}
	//for (i=1; i<=asarray_elements(counter); i++) {
	//	printf("{txt} - key=%s count=%f\n", keepvars[i], asarray(counter,keepvars[i]))
	//}

	// Main loop
	i = i_last_singleton = g = 1
	while (i<i_last_singleton+G) {
		if (g>G) g = 1
		if (verbose>0) printf("{txt}\ti=%f (g=%f/%f)\t(N=%f)\t", i, g, G, st_nobs())

		idvarnames = i<=G ? fes[g].ivars : fes[g].idvarname
		id = st_data(., idvarnames) // 2% of runtime
		if (i<=G) {
			for (j=1; j<=length(idvarnames); j++) {
				n = asarray(counter, idvarnames[j]) - 1
				asarray(counter, idvarnames[j], n)
				if (n==0) {
					st_dropvar(idvarnames[j])
				}
			}
		}
		if (i<=G) fes[g].is_sortedby = already_sorted(idvarnames)
		sortedby = fes[g].is_sortedby
		if (i<=G & !sortedby) {
			fes[g].p = order( id , 1..length(idvarnames) ) // 55% of function time
		}
		
		if (!sortedby) {
			_collate(id, fes[g].p) // sort id by p // 12% of function time
			inv_p = invorder(fes[g].p) // construct inv(p) that we'll use later
		}

		// Note that the lhs is actually the deltas, as in "bys id: gen delta = _n==1"
		id = rows_that_change(id) // 7% of function time
		singleton = select_singletons(id) // 5% of function time

		// Save IDs in dataset before dropping observations
		id = runningsum(id :* !singleton) // this is the ID now, not the deltas anymore
		fes[g].levels = id[length(id)]
		vartype = fes[g].levels<=100 ? "byte" : (fes[g].levels<=32740? "int" : "long")
		if (i<=G) {
			st_store(., st_addvar(vartype, fes[g].idvarname), sortedby? id : id[inv_p])
		}
		else {
			st_store(., fes[g].idvarname, sortedby? id : id[inv_p])
		}

		num_singletons = sum(singleton)
		if (num_singletons>0) {
			if (verbose>0) printf("{txt}(%f singletons)", num_singletons)
			i_last_singleton = i

			// Sort -singleton- as in the dataset, and use it to drop observations
			// 5% of function time
			singleton = sortedby? singleton : singleton[inv_p]
			st_dropobsif(singleton)
			if (!st_nobs()) {
				printf("{err}\nno observations left after dropping singletons\n")
				exit(error(2001))
			}

			// But now our precious sort orders (the p's) are useless! Fix them
			sum_singleton = runningsum(singleton)
			for (h=1;h<=G & h<=i; h++) { // 6% of function time
				if (fes[h].is_sortedby) continue
				pp = &(fes[h].p)
				(*pp) = select(*pp - sum_singleton[*pp] , !singleton[*pp] )
			}
		}
		if (verbose>0) printf("{txt}\n")

		i++
		g++
	}
}

// --------------------------------------------------------------------------------------
// (Auxiliary Functions)
// --------------------------------------------------------------------------------------

// ALREADY_SORTED:
`Integer' already_sorted(string vector vars) {
	`Varlist' sortedby
	sortedby = tokens(st_macroexpand("`" + ": sortedby" + "'"))
	return(length(vars) > length(sortedby) ? 0 : vars==sortedby[1..length(vars)])
}

// --------------------------------------------------------------------------------------

// ROWS_THAT_CHANGE: Return a 0/1 vector indicating what are different from the previous row
	// Idea: compromise between doing the operation in one go (uses lots of memory) vs doing loops (2x slower)
	// Other alternatives: loop row-by-row (slower), loop col-by-col and take max() (slower)
`Vector' rows_that_change(`Matrix' input) {
	`Vector' ans
	`Integer' i, j, K, N, stepsize
	
	// Size of blocks of matrices used (larger=faster smaller=less memory)
	// Benchmarks with 3 unsorted ivars showed 1e4 was fastest, followed by 1e3 and then 1e5
	stepsize = 1e4

	N = rows(input)
	K = cols(input)
	ans = J(N,1,0)
	ans[1] = 1
	for (i=2; i<=N;i=i+stepsize) {
		j = min((i+stepsize-1, N))
		ans[|i\j|] = rowmax(input[|i-1,1\j-1,K|] :!= input[|i,1\j,K|])
	}
	return(ans)
}

// --------------------------------------------------------------------------------------

// SELECT_SINGLETONS: 
`Vector' select_singletons(`Vector' input) {
	// Code modified from <rows_that_change>
	`Vector' ans
	`Integer' i, j, K, N, stepsize

	// Size of blocks of matrices used (larger=faster smaller=less memory)
	// Benchmarks with 3 unsorted ivars showed 1e4 was fastest, followed by 1e3 and then 1e5
	stepsize = 1e4

	N = rows(input)
	ans = J(N,1,0)
	for (i=1; i<=N-1;i=i+stepsize) {
		j = min((i+stepsize-1, N-1))
		// We need that ans[i]==1 and ans[i+1]==1
		// Since ans is either 0 or 1, this is equivalent to
		ans[|i\j|] = (input[|i\j|] + input[|i+1\j+1|] :== 2) 
	}
	ans[N] = (input[N]==1) // Special case, last obs is singleton if it's the first obs in the group
	return(ans)
}

end
