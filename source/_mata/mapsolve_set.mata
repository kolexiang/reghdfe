mata:
mata set matastrict on

// Initialize an almost empty structure; run this right after ParseAbsvars.ado to use their r(..)
function mapsolve_set(`Problem' S) {
	// string vector				absvars
	// string scalar				token
	`Integer'					g, G
	`Integer'					G_expanded // Counts fixed slopes as independent categories
	transmorphic				t

	// fe_parse(absvars)
	fe_init(absvars)
	// later split fe_parse
	return(S)
}

end
