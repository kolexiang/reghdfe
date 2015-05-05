mata:
mata set matastrict on

`Integer' function map_connected_groups(`Problem' S, `Integer' g1, `Integer' g2, | `Varname' groupvar) {
	`Boolean' changed
	`Series' group, p
	`Integer' gg, g, j, i_lower, i_upper, num_groups, L
	real rowvector min_max
	`String' vartype

	changed = 1
	group = st_data(., S.fes[g1].idvarname)
	
	while (changed) {
		changed = 0
		for (gg=1;gg<=2;gg++) {
			g = gg==1 ? g2 : g1
			L = S.fes[g].levels
			if (!S.fes[g].is_sortedby) _collate(group, S.fes[g].p) // Sort it by g1 or g2
			i_lower = 1
			for (j=1;j<=L; j++) {
				i_upper = S.fes[g].offsets[j]
				min_max = minmax(group[| i_lower , 1 \ i_upper , 1 |])
				if (min_max[1]!=min_max[2]) changed = 1
				group[| i_lower , 1 \ i_upper , 1 |] = min_max[1] :* J(i_upper-i_lower+1,1,1)
				i_lower = i_upper + 1
			}
			if (!S.fes[g].is_sortedby) _collate(group, S.fes[g].inv_p) // Sort it back
		}
	}

	// Create compact group id
	p = order(group, 1)
	_collate(group, p)
	group = runningsum(rows_that_change(group))
	num_groups = group[rows(group)]
	_collate(group, invorder(p))
	
	// (optional) save group variable
	if (groupvar!="") {
		vartype = num_groups<=100 ? "byte" : (num_groups<=32740? "int" : "long")
		st_store(., st_addvar(vartype, groupvar), group)
		if (S.verbose>2) printf("{txt}    - Saving identifier for the first mobility group: {res}%s\n", groupvar)
		st_varlabel(groupvar, sprintf("Mobility Group: %s <--> %s", invtokens(S.fes[g1].ivars,"#") , invtokens(S.fes[g2].ivars,"#")))
	}
	return(num_groups)
}
end
