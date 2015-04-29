mata:
mata set matastrict on

void function map_precompute(`Problem' S) {
	`Integer' i, g, G
	`Varlist' keepvars
	transmorphic counter
	if (S.verbose>0) printf("{txt}mata: map_precompute()\n")

	// Count how many times each var is used, so we can drop them when the counter reaches zero
	G = length(S.fes)
	counter = asarray_create()
	asarray_notfound(counter, 0)
	keepvars = S.keepvars
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
	if (S.verbose>0) printf("{txt} - Storing permutation vectors, generating ids, dropping singletons\n")
	map_precompute_part1(S, counter)

	// 2. Store group offsets, group counters; demeaned(x), inv(xx) if num_slopes>0; weightvars
	if (S.verbose>0) printf("{txt} - Storing counters and offsets; processing cvars\n")
	map_precompute_part2(S, counter)

	// 3. Optionally drop IDs from data set and store precomputed inv(p)
	// We needed those IDs in part2
	if (S.verbose>0) printf("{txt} - Storing reverse permutation vectors, optionally dropping IDs\n")
	for (g=1;g<=G;g++) {
		if (!S.save_ids) st_dropvar(S.fes[g].idvarname)
		S.fes[g].inv_p = invorder(S.fes[g].p)
	}

	// Store N (todo: add other details) to ensure the dataset doesn't change from now on
	S.N = st_nobs()
}
end
