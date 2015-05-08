mata:
mata set matastrict on

real rowvector safe_divide(real rowvector numerator, real rowvector denominator, | real scalar epsi) {
	 // If the denominator goes below machine precision, the division explodes
	 if (args()<3) epsi = epsilon(1)
	return( numerator :/ colmax(denominator \ J(1,cols(denominator),epsi)) )
}

// -------------------------------------------------------------------------------------------------

void verbose2local(`Problem' S, string scalar loc) {
	st_local(loc, strofreal(S.verbose))
}

void function store_uid(`Problem' S, `Varname' varname) {
	S.uid = st_data(., varname)
	assert_msg(rows(S.uid)==S.N, "assertion failed: rows(S.uid)==S.N")
}

void function groupvar2dta(`Problem' S) {
	if (S.groupvar!="") {
		if (S.verbose>2) printf("{txt}    - Saving identifier for the first mobility group: {res}%s\n", S.groupvar)
		// WRONG: NEED TO INDEX BY UID
		st_store(S.uid, st_addvar(S.grouptype, S.groupvar), S.groupseries)

		S.groupseries = J(0,0,0)
		st_varlabel(S.groupvar, S.grouplabel)
	}
	else {
		if (S.verbose>=0) printf("{txt}Note: mobility group not saved as there was no need to compute it (e.g. there are less than two absvars not nested within a cluster)\n")
	}
}

void function drop_ids(`Problem' S) {
	`Integer' g
	for (g=1;g<=S.G;g++) {
		if (!S.fes[g].is_clustervar) st_dropvar(S.fes[g].idvarname)
	}
}

end
