mata:
mata set matastrict on

real rowvector safe_divide(real rowvector numerator, real rowvector denominator) {
	 // If the denominator goes below machine precision, the division explodes
	return( numerator :/ colmax(denominator \ J(1,cols(denominator),epsilon(1))) )
}

// -------------------------------------------------------------------------------------------------

void verbose2local(`Problem' S, string scalar loc) {
	st_local(loc, strofreal(S.verbose))
}

end
