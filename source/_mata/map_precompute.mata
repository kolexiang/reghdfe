mata:
mata set matastrict on

void function map_precompute(`Problem' S, `Varlist' keepvars) {
	`Integer' i, g, G
	transmorphic counter
	if (S.verbose>0) printf("{txt}mapsolve_precompute()\n")

	// Count how many times each var is used, so we can drop them when the counter reaches zero
	G = length(S.fes)
	counter = asarray_create()
	asarray_notfound(counter, 0)
	for (g=1; g<=G; g++) {
		keepvars = keepvars , S.fes[g].ivars , S.fes[g].cvars
	}
	for (i=1; i<=length(keepvars); i++) {
		asarray(counter, keepvars[i], asarray(counter, keepvars[i])+1)
	}
	//for (i=1; i<=asarray_elements(counter); i++) {
	//	printf("{txt} - key=%s count=%f\n", keepvars[i], asarray(counter,keepvars[i]))
	//}

	// 1. Store permutation vectors and their invorder, generate ID variables, drop singletons
	if (S.verbose>0) printf("{txt}    Storing permutation vectors, generating ids, dropping singletons\n")
	map_precompute_part1(S, counter)

	// 2. Store group offsets, group counters; demeaned(x), inv(xx) if num_slopes>0; weightvars
	if (S.verbose>0) printf("{txt}    Storing counters and offsets; processing cvars\n")
	map_precompute_part2(S, counter)

	// Store N (todo: add other details) to ensure the dataset doesn't change from now on
	S.N = st_nobs()
}
end
