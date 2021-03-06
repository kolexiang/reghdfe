capture program drop Stop
program define Stop
	cap mata: mata drop prev_numstep // Created at step 1
	cap mata: mata drop VERBOSE // Created before step 1
	cap mata: mata drop G // Num of absorbed FEs
	cap mata: mata drop FEs // Main Mata structure
	cap mata: mata drop betas // Temporary matrices used to store bi/multivariate regr coefs
	cap mata: mata drop varlist_cache // Hash table with the names of the precomputed residuals
	cap mata: mata drop avge_* // Drop AvgE structures
	cap mata: mata drop weightexp weightvar

	cap mata: mata drop clustervars
	cap mata: mata drop clustervars_ivars
	cap mata: mata drop clustervars_original

	if ("${hdfe_pwd}"!="") {
		qui cd "${hdfe_pwd}"
		global hdfe_pwd
	}

	* PARALLEL SPECIFIC CLEANUP
	cap mata: st_local("path", parallel_path)
	if ("`path'"!="") {
		mata: st_local("cores", strofreal(parallel_cores))
		assert "`cores'"!=""
		local path "`path'"
		cap erase `"`path'hdfe_mata.mo"'
		forv core=1/`cores' {
			cap erase `"`path'`core'_done.txt"'
			cap erase `"`path'`core'_ok.txt"'
			cap erase `"`path'`core'_error.txt"'
			cap erase `"`path'`core'_output.dta"'
			cap erase `"`path'`core'_log.log"'
		}
		cap rmdir `"`path'"'
		cap mata: mata drop parallel_cores
		cap mata: mata drop parallel_dta
		cap mata: mata drop parallel_vars
		cap mata: mata drop parallel_opt
		cap mata: mata drop parallel_path
	}
end
