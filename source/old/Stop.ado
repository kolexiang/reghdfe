capture program drop Stop
program define Stop
	cap mata: mata drop HDFE_S // Structure for the HDFE problem solved with map_solve()
	cap mata: mata drop varlist_cache // Hash table with the names of the precomputed residuals
	cap mata: mata drop avge_* // Drop AvgE structures
	// BUGBUG: Maybe we can remove all but HDFE_S, and then we can remove Stop and just do cap mata: mata drop HDFE_S ?
	}
end
