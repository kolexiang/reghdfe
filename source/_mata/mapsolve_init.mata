mata:
mata set matastrict on

//  Parse absvars and initialize the almost empty MapProblem struct
`Problem' function mapsolve_init() {
	`Problem' 	S
	// string vector				absvars
	// string scalar				token
	`FE' fe
	`Integer'					g, G
	`Integer'					G_expanded // Counts fixed slopes as independent categories

	G = st_numscalar("r(G)")
	S.fixed_effects = FixedEffect(G)
	for (g=1; g<=G; g++) {
		fe = S.fixed_effects[g]
		fe.has_intercept = st_numscalar(sprintf("r(G%1.0f)",g))
		S.fixed_effects[g] = fe
	}
	return(S)
}

end
