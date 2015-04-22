mata:
mata set matastrict on

struct MapProblem {
	struct FixedEffect vector fixed_effects
	`Integer'			G 			// Number of FEs when bunching slopes
	`Integer'			G_expanded 	// Number of FEs incl. slopes
	`Varname'			weightvar 	// Variable contaning the fw/pw/aw
}
end
