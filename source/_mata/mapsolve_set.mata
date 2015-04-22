mata:
mata set matastrict on

// Fill the auxiliary vectors (count, etc.) and remove singletons
void function mapsolve_set(`Problem' S) {
	`Integer' G, g, N, i, j, K, levels, stepsize
	`Series' p, id
	`Group'  ivars
	real vector last_row
	string scalar vartype, newvarname

	G = S.G
	stepsize = 1e6

	// Generate ID Variables
	for (g=1;g<=G;g++) {
		timer_clear()
		111
		timer_on(8)
		ivars = st_data(., S.fixed_effects[g].ivars)
		N = rows(ivars)
		K = cols(ivars)
		timer_off(8)
		timer_on(1)
		1112
		p = order(ivars, 1..cols(ivars))
		1113
		timer_off(1)
		timer_on(2)
		_collate(ivars, p) // ivars = ivars[p, .]
		timer_off(2)
		1114
		timer_on(3)
		timer_off(3)
		1115
		222
		timer_on(4)

		// Idea: compromise between doing op in one go (uses lots of memory) vs doing loops (2x slower)
		id = J(N,1,0)
		id[1] = 1
		assert(N>=2)
		for (i=2; i<=N;i=i+stepsize) {
			j = min((i+stepsize-1, N))
			id[|i\j|] = rowmax(ivars[|i-1,1\j-1,K|] :!= ivars[|i,1\j,K|])
		}
		id = runningsum(id)

		// id = J(N, 1, 0)
		// last_row = J(1,K,.)
		// for (i=1;i<=N;i++) {
		// 	if (last_row!=ivars[|i,1\i,K|]) {
		// 		last_row = ivars[|i,1\i,K|]
		// 		id[i] = 1
		// 	}
		// }

		//for (i=1;i<=K;i++) {
		//	id[|2\N|] = rowmax((id[|2\N|],  ivars[|1,i \ N-1,i|] :!= ivars[|2,i \ N,i|] ))
		//}
		timer_off(4)
		333
		timer_on(5)
		levels = id[N]
		timer_off(5)
		timer_on(6)
		vartype = levels<=100 ? "byte" : (levels<=32740? "int" : "long")
		newvarname = "FE" + strofreal(g)
		444
		st_store(., st_addvar(vartype, newvarname), id[invorder(p)])
		timer_off(6)
		555
		// S.fixed_effects[g].ivars
		// id, ivars
		//ivars
	}
	timer()
}

end
