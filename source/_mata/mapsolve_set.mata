mata:
mata set matastrict on

// Fill the auxiliary vectors (count, etc.) and remove singletons
void function mapsolve_set(`Problem' S) {
	`Integer' G, g, N, i, j, K, levels, stepsize, num_singletons, last_num_singletons
	`Series' p, id, is_singleton, first_pass
	`Group'  ivars
	`Varlist' ivar_names
	real vector last_row
	string scalar vartype, newvarname

	G = S.G
	stepsize = 1e6 // size of blocks of matrices used (larger=faster smaller=less memory)
	last_num_singletons = 0
	first_pass = 1
	// singletons = J(N,1,1)

	// Generate ID Variables
	for (g=1;g<=G;g++) {
		ivar_names = S.fixed_effects[g].ivars
		ivars = st_data(., ivar_names) // Profiler: 1%
		N = rows(ivars)
		K = cols(ivars)

		if (already_sorted(ivar_names)) {
			"already sorted"
			p = 1::N
		}
		else {
			p = order(ivars, 1..K) // Profiler 80%
			_collate(ivars, p) // ivars = ivars[p, .] // Profiler: 8%
		}
		
		// Idea: compromise between doing op in one go (uses lots of memory) vs doing loops (2x slower)
		// Alternatives: loop row-by-row (slower), loop col-by-col and take max() (slower)

		id = J(N,1,0)
		id[1] = 1
		assert(N>=2)
		for (i=2; i<=N;i=i+stepsize) {
			j = min((i+stepsize-1, N))
			id[|i\j|] = rowmax(ivars[|i-1,1\j-1,K|] :!= ivars[|i,1\j,K|]) // Profiler: 6%
		}



		id = runningsum(id)

		levels = id[N]
		vartype = levels<=100 ? "byte" : (levels<=32740? "int" : "long")
		newvarname = "FE" + strofreal(g)
		st_store(., st_addvar(vartype, newvarname), id[invorder(p)]) // Profiler: 5%
		
		// S.fixed_effects[g].ivars
		// id, ivars
		//ivars
	}

	// Ya tengo los IDs
	// Ya los exporte

	// Que quiero hacer ahora?
	// 	Meter weights es facil pero hacerlo al final cuando lo necesite
	// Quiero quedarme con los IDs que no son singletons
	// hacerlo iterativamente z

}

real scalar already_sorted(string vector vars) {
	`Varlist' sortedby
	sortedby = tokens(st_macroexpand("`" + ": sortedby" + "'"))
	if (length(vars) <= length(sortedby)) sortedby[1..length(vars)]
	return(length(vars) > length(sortedby) ? 0 : vars==sortedby[1..length(vars)])
}

end


// Usually data is already sorted by the first FE
// In that case, mata is EXTREMELY slow b/c it doesnt know it
// What we should do is check if the sort is already by ivars (and something else, but starting with the ivars)
// in that case, then order is just 1..N and collate is not needed
