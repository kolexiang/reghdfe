mata:
mata set matastrict on

//  Parse absvars and initialize the almost empty MapProblem struct
`Problem' function mapsolve_init(`Integer' verbose) {
	`Problem' 	S
	pointer(`FE') 	fe
	`Integer'		g, G
	`Integer'		G_expanded // Counts fixed slopes as independent categories

	G = st_numscalar("r(G)")
	S.G = G
	S.verbose = verbose
	
	if (S.verbose>0) printf("{txt}mapsolve_init()  initializing MapProblem structure\n")
	S.fixed_effects = FixedEffect(G)
	for (g=1; g<=G; g++) {
		fe = &(S.fixed_effects[g])
		// recall a->b is the same as (*a).b
		fe->order = g
		fe->varlabel = st_global(sprintf("r(varlabel%f)",g))
		fe->num_slopes = st_numscalar(sprintf("r(num_slopes%f)",g))
		fe->has_intercept = st_numscalar(sprintf("r(has_intercept%1.0f)",g))
		fe->ivars = tokens(st_global( sprintf("r(ivars%f)",g) ))
		fe->cvars = tokens(st_global( sprintf("r(cvars%f)",g) ))
		fe->target = st_global(sprintf("r(target%f)",g))
		fe->idvarname = sprintf("__ID%f__", g)
		fe->levels = .
	}

	return(S)
}

end