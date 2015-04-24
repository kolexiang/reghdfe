mata:
mata set matastrict on

// Fill the auxiliary vectors (count, etc.) and remove singletons
void function mapsolve_set(`Problem' S) {
	`Integer' G, g, N, i, j, K, levels, num_singletons, last_num_singletons, stepsize
	`Series' p, id, is_singleton, first_pass
	`Group'  ivars
	`Varlist' ivar_names
	real vector last_row
	string scalar vartype, newvarname

	// Create a Count() object (hash table with counter) with clustervars + panelvar +timevar + all ivars

	i = g = 1
	i_last_singleton = 0
	G = S.G

	while (i<=G | i<i_last_singleton+G) {
		if (i<=G) {
			ivars = S.fixed_effects[g].ivars
			id = st_data(., ivars)
		}
		
		// EQUIVALENT TO: by id: gen delta = (_n==1)
		// EQUIVALENT TO: gen delta = id!=id[_n-1]
		delta = rows_that_change(id)

		// EQUIVALENT TO: by id: gen singleton = _N==1
		// EQUIVALENT TO: gen singleton = delta!=delta[_n-1]
		singleton = rows_that_change(delta) // equivalent to: 

		if (i<=G) {
			// REMOVE ivars that are now redundant
			// purge_variables(count_table, ivars)
		}


		i++
		g = i>g : 1 : i
	}



	last_num_singletons = 0
	first_pass = 1
	// singletons = J(N,1,1)

	// Generate ID Variables
	for (g=1;g<=G;g++) {
		ivar_names = 
		ivars =  // Profiler: 1%
		N = rows(ivars)
		K = cols(ivars)

		if (already_sorted(ivar_names)) {
			p = 1::N
		}
		else {
			p = order(ivars, 1..K) // Profiler 80%
			_collate(ivars, p) // ivars = ivars[p, .] // Profiler: 8%
		}
		



		// no puedo hacer esto xq de repente el singleton esta a la mitad del pata
		// id = id :& singleton //

		// singleton = singleton :& 

		// me lo bajo asi!! dropit[invorder(p)]

		// esto es igual que dropit[invorder(p)]
		tag = J(N,1,0)
		tag[select(p,dropit)] = J(sum(dropit),1,1) // -tag- tiene a los q me bajo

		hasta runningsum(tag) puedo hacerlo comun para todos!
		p = select(p - runningsum(tag)[p], !dropit)



		id = runningsum(id) // meterle un mask a lo :* (!singleton) asi id=0 significa signleton

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

`Integer' already_sorted(string vector vars) {
	`Varlist' sortedby
	sortedby = tokens(st_macroexpand("`" + ": sortedby" + "'"))
	return(length(vars) > length(sortedby) ? 0 : vars==sortedby[1..length(vars)])
}

// -------------------------------------------------------------------------------------------------
// Return a 0/1 vector indicating what are different from the previous row
// -------------------------------------------------------------------------------------------------
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

end


