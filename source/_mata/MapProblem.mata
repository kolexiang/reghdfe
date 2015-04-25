mata:
mata set matastrict on

struct MapProblem {
	struct FixedEffect vector fixed_effects
	`Integer'			G 			// Number of FEs when bunching slopes
	`Integer'			G_expanded 	// Number of FEs incl. slopes
	`Varname'			weightvar 	// Variable contaning the fw/pw/aw
	`Varlist'			clustervars	// Base vars of the clustervars (e.g. x#y z -> x y z) +- just generate the ids
	`Integer'			verbose		// Number of debug messages to show (0=None, 1=A little, 4=A *lot*)			
}
end