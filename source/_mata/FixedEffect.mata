mata:
mata set matastrict on
	struct FixedEffect {
		`Integer'	order 			// "g", the position in varlist
		`Varname'	varlabel		// Original label of this absvar
		`Integer'	num_slopes
		`Integer'	has_intercept
		`Integer'	levels			// Number of categories spanned by the ivars
		`Varlist'	ivars			// number of i.var elements
		`Varlist'	cvars			// number of c.var elements or slopes
		`Series'	target			// Name of the variable that will hold the estimates for the FE
		`Boolean'	is_sortedby		// 1 if the dataset is sorted by the ivars of this FE
		`Series'	p 				// Permutation vector
		`Varname'	idvarname		// (optional) Name of variable with the absvar categories
	}
end
