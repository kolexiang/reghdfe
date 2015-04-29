mata:
mata set matastrict on
// -------------------------------------------------------------------------------------------------
// Transformations: Compute RESIDUALS, not projections
// -------------------------------------------------------------------------------------------------

void function transform_cimmino(`Problem' S, `Group' y, `Group' ans,| `Boolean' get_proj) {
	`Integer' 	N, Q, g, G
	N = rows(y)
	Q = cols(y)
	G = S.G
	if (args()<4) get_proj = 0

	ans = map_project(S, 1, y)
	for (g=2; g<=G; g++) {
		ans = ans + map_project(S, g, y)
	}
	ans = get_proj ? ans / G : y - ans / G
}

// -------------------------------------------------------------------------------------------------

void function transform_kaczmarz(`Problem' S, `Group' y, `Group' ans,| `Boolean' get_proj) {
	`Integer' 	N, Q, g, G
	N = rows(y)
	Q = cols(y)
	G = S.G
	if (args()<4) get_proj = 0

	ans = y - map_project(S, 1, y)
	for (g=2; g<=G; g++) {
		ans = ans - map_project(S, g, ans)
	}
	if (get_proj) ans = y - ans
}

// -------------------------------------------------------------------------------------------------

 void function transform_sym_kaczmarz(`Problem' S, `Group' y, `Group' ans,| `Boolean' get_proj) {
	`Integer' 	N, Q, g, G
	// BUGBUG: Streamline and remove all those "ans - .." lines?
	N = rows(y)
	Q = cols(y)
	G = S.G
	if (args()<4) get_proj = 0

	ans = y - map_project(S, 1, y)
	for (g=2; g<=G; g++) {
		ans = ans - map_project(S, g, ans)
	}
	for (g=G-1; g>=1; g--) {
		ans = ans - map_project(S, g, ans)
	}
	if (get_proj) ans = y - ans
}
end
