*! hdfe 3.0.809 07may2015
*! Sergio Correia (sergio.correia@duke.edu)


// Mata code is first, then main hdfe.ado, then auxiliary .ado files
// -------------------------------------------------------------------------------------------------
// Mata Code: Method of Alternating Projections with Acceleration
// -------------------------------------------------------------------------------------------------
	//discard
	//pr drop _all
	//clear all

* Type Aliases
	local Boolean 		real scalar
	local Integer 		real scalar
	local Real 			real scalar
	local Vector		real colvector
	local Matrix		real matrix
	local Series		real colvector // Should be N*1
	local Group			real matrix // Should be N*K
	local String 		string scalar
	local Varname 		string scalar
	local Varlist 		string rowvector // rowvector so they match tokens()
	local Problem		struct MapProblem scalar
	local FE			struct FixedEffect scalar
	local FunctionPointer pointer(`Group' function) scalar // Used for the Accelerate & Transform fns

// -------------------------------------------------------------------------------------------------
	
mata:
mata set matastrict on
	void assert_msg(real scalar t, | string scalar msg)
	{
		if (args()<2 | msg=="") msg = "assertion is false"
	        if (t==0) _error(msg)
	}
	
// -------------------------------------------------------------------------------------------------
// Structure of a Fixed Effect
// -------------------------------------------------------------------------------------------------

	struct FixedEffect {
		`Integer'	order 			// "g", the position in varlist
		`Varname'	varlabel		// Original label of this absvar
		`Integer'	num_slopes
		`Integer'	has_intercept
		`Integer'	levels			// Number of categories spanned by the ivars
		`Varlist'	ivars			// number of i.var elements
		`Varlist'	cvars			// number of c.var elements or slopes
		`Boolean'	is_sortedby		// 1 if the dataset is sorted by the ivars of this FE
		`Varname'	idvarname		// (optional) Name of variable with the absvar categories
		`Varlist'	target			// Name of the variable that will hold the estimates for the FE
		
		`Series'	p 				// Permutation vector
		`Series'	inv_p 			// Precompute invorder(p) because we'll use it a lot

		`Vector'	offsets			// Pos. of last obs. with that category (when sorted)
		`Vector'	counts			// number of obs. (weighted) with that value
		`Group'		x				// Vector/Matrix of (optionally demeaned) cvars
		`Matrix'	inv_xx			// Blocks of the inv(x'x) matrix; size KL*K (k=num_slopes, L=levels)

		`Boolean'	is_clustervar, in_clustervar
		`Integer'	nesting_clustervar // Clustervar that nests this FE, if any

		// Temporary matrices for the stored FEs
		`Matrix'	alphas
		`Matrix'	tmp_alphas
	}
	
// -------------------------------------------------------------------------------------------------
// Structure of the MAP Problem
// -------------------------------------------------------------------------------------------------

struct MapProblem {
	struct FixedEffect vector fes	// The G*1 vector of FE structures
	`Integer'		G 				// Number of FEs when bunching slopes
	`Varname'		weightvar 		// Name variable contaning the fw/pw/aw
	`String'		weighttype 		// fw/pw/aw
	`String'		weights 		// "[weighttype=weightvar]"
	`Series'		w 				// Contents of variable contaning the fw/pw/aw
	`Integer'		verbose			// Number of debug messages to show (0=None, 1=A little, 4=A *lot*)			
	`Integer'		N				// Number of obs; after map_precompute() the dataset CANNOT CHANGE!
	`Varlist'		keepvars		// By default we drop cvars and ivars ASAP; this prevents it (useful for clustervars and for timevar+panelvar under HAC errors)

	`Boolean'		will_save_fe	// True if at least one FE will be saved
	`Boolean'		keepsingletons	// If set to 1, it will not drop singletons (do not touch this!)

	`Integer'		C				// Number of cluster variables
	`Varlist'		clustervars
	`Varlist'		clustervars_original // Note: need to apply tokens()
	`Varname'		panelvar
	`Varname'		timevar
	`Boolean'		vce_is_hac
	
	// Optimization parameters	
	`Integer'		groupsize 		// Group variables when demeaning (more is faster but uses more memory)
	`Real'			tolerance
	`Integer'		maxiterations
	`String'		transform		// Kaczmarz Cimmino Symmetric_kaczmarz (k c s)
	`String'		acceleration	// Acceleration method. None/No/Empty is none\
	`Integer'		accel_start		// Iteration where we start to accelerate // set it at 6? 2?3?
	
	// Specific to Aitken's acceleration
	`Integer'		accel_freq		
	`Integer'		stuck_threshold	// Call the improvement "slow" when it's less than e.g. 1%
	`Integer'		bad_loop_threshold	// If acceleration seems stuck X times in a row, pause it
	`Integer'		pause_length	// This is in terms of candidate accelerations, not iterations (i.e. x3)?

	// Temporary
	`Boolean'		storing_betas
}
	
real rowvector safe_divide(real rowvector numerator, real rowvector denominator, | real scalar epsi) {
	 // If the denominator goes below machine precision, the division explodes
	 if (args()<3) epsi = epsilon(1)
	return( numerator :/ colmax(denominator \ J(1,cols(denominator),epsi)) )
}

// -------------------------------------------------------------------------------------------------

void verbose2local(`Problem' S, string scalar loc) {
	st_local(loc, strofreal(S.verbose))
}
	
//  Parse absvars and initialize the almost empty MapProblem struct
`Problem' function map_init()
{
	`Integer'		g, G, num_slopes, has_intercept, i
	`Problem' 		S
	`Boolean'		auto_target // Automatically assign target names to all FEs
	`Varname'		basetarget
	`Varlist'		target
	pointer(`FE') 	fe

	S.weightvar = S.weighttype = S.weights = ""
	S.verbose = 0
	S.transform = "cimmino"
	S.acceleration = "conjugate_gradient"
	S.tolerance = 1e-7
	S.maxiterations = 1e4
	S.accel_start = 6
	S.groupsize = 10

	// If clustering by timevar or panelvar and VCE is HAC, we CANNOT touch the clustervars to create compact ids!
	S.timevar = ""
	S.panelvar = ""
	S.vce_is_hac = 0

	// Specific to Aitken:
	S.accel_freq = 3
	S.pause_length = 20
	S.bad_loop_threshold = 1
	S.stuck_threshold = 5e-3
	S.N = .

	S.keepsingletons = 0
	S.G = G = st_numscalar("r(G)")
	S.C = 0
	S.clustervars = S.clustervars_original = J(0,0,"")
	S.fes = FixedEffect(G)
	S.will_save_fe = 0
	auto_target = st_numscalar("r(savefe)")
	assert(auto_target==1 | auto_target==0)
	if (auto_target) stata(sprintf("cap drop __hdfe*__*"))

	for (g=1; g<=G; g++) {
		fe = &(S.fes[g])
		// recall a->b is the same as (*a).b
		num_slopes = st_numscalar(sprintf("r(num_slopes%f)",g))
		has_intercept = st_numscalar(sprintf("r(has_intercept%1.0f)",g))

		fe->order = g
		fe->num_slopes = num_slopes
		fe->has_intercept = has_intercept
		fe->varlabel = st_global(sprintf("r(varlabel%f)",g))
		fe->ivars = tokens(st_global( sprintf("r(ivars%f)",g) ))
		fe->cvars = tokens(st_global( sprintf("r(cvars%f)",g) ))
		fe->idvarname = sprintf("__ID%f__", g)
		fe->levels = .
		fe->target = J(0,0,"")
		fe->is_clustervar = 0
		fe->in_clustervar = 0
		fe->nesting_clustervar = .
		
		basetarget = st_global(sprintf("r(target%f)",g))
		if (basetarget=="" & auto_target) basetarget = sprintf("__hdfe%f__", g)
		if (basetarget!="") {
			S.will_save_fe = 1
			target = J(1, has_intercept + num_slopes, basetarget)
			if (has_intercept) stata(sprintf("confirm new variable %s", target[1]))
			for (i=1+has_intercept; i<=length(target); i++) {
				target[i] = target[i] + sprintf("_Slope%f", i-has_intercept)
				stata(sprintf("confirm new variable %s", target[i]))
			}
			fe->target = target
		}
	}
	
	st_numscalar("r(save_fe)", S.will_save_fe)
	return(S)
}

void function map_init_clustervars(`Problem' S, `String' clustervars) {
	S.clustervars = S.clustervars_original = tokens(clustervars)
	S.C = length(S.clustervars)
}

void function map_init_weights(`Problem' S, `Varname' weightvar, `String' weighttype) {
	assert_msg(weightvar!="" & weighttype!="", "map_init_weights() requires weight var and type")
	stata(sprintf("confirm numeric variable %s, exact", weightvar))
	assert_msg(anyof(("fweight", "pweight", "aweight"), weighttype), "wrong weight type")
	S.weightvar = weightvar
	S.weighttype = weighttype
	S.weights = sprintf("[%s=%s]", weighttype, weightvar)
}

void function map_init_keepvars(`Problem' S, `Varname' keepvars) {
	if (keepvars!="") stata(sprintf("confirm numeric variable %s, exact", keepvars))
	S.keepvars = tokens(keepvars)
}

void function map_init_transform(`Problem' S, `String' transform) {
	transform = strlower(transform)
	// Convert abbreviations
	if (strpos("cimmino", transform)==1) transform = "cimmino"
	if (strpos("kaczmarz", transform)==1) transform = "kaczmarz"
	if (strpos("symmetric_kaczmarz", transform)==1) transform = "symmetric_kaczmarz"
	assert_msg(anyof(("cimmino","kaczmarz","symmetric_kaczmarz"),transform), "invalid transform")
	S.transform = transform
}

void function map_init_verbose(`Problem' S, `Integer' verbose) {
	assert_msg(round(verbose)==verbose, "verbose must be an integer")
	assert_msg(0<=verbose & verbose<=5, "verbose must be between 0 and 5")
	S.verbose = verbose
}

void function map_init_groupsize(`Problem' S, `Integer' groupsize) {
	assert_msg(round(groupsize)==groupsize & groupsize>0, "groupsize must be a positive integer")
	S.groupsize = groupsize
}

void function map_init_acceleration(`Problem' S, `String' acceleration) {
	acceleration = strlower(acceleration)
	// Convert abbreviations
	if (strpos("conjugate_gradient", acceleration)==1 | acceleration=="cg") acceleration = "conjugate_gradient"
	if (strpos("steepest_descent", acceleration)==1 | acceleration=="sd") acceleration = "steepest_descent"
	if (strpos("aitken", acceleration)==1) acceleration = "aitken"
	if (acceleration=="no" | acceleration=="none" | acceleration=="off") acceleration = "none"
	assert_msg(anyof(("conjugate_gradient","steepest_descent", "aitken", "none"),acceleration), "invalid acceleration")
	S.acceleration = acceleration
}

void function map_init_tolerance(`Problem' S, `Real' tolerance) {
	assert_msg(1e-16<=tolerance & tolerance<1, "tolerance must be in the range [1e-16, 1).")
	S.tolerance = tolerance
}

void function map_init_maxiterations(`Problem' S, `Integer' maxiterations) {
	assert_msg(round(maxiterations)==maxiterations, "maxiterations must be an integer")
	assert_msg(maxiterations>0, "maxiterations must be positive")
	S.maxiterations = maxiterations
}

void function map_init_keepsingletons(`Problem' S, `Boolean' keepsingletons) {
	assert_msg(keepsingletons==0 | keepsingletons==1, "keepsingletons must be 0 or 1")
	S.keepsingletons = keepsingletons
}

void function map_init_vce_is_hac(`Problem' S, `Boolean' vce_is_hac) {
	assert_msg(vce_is_hac==0 | vce_is_hac==1, "vce_is_hac must be 0 or 1")
	S.vce_is_hac = vce_is_hac
}
	
void function map_precompute(`Problem' S) {
	`Integer' i, g, h, value
	`Varlist' keepvars, cl_ivars
	transmorphic counter, loc
	`Varname' key
	if (S.verbose>0) printf("\n{txt}{bf:mata: map_precompute()}\n")

	// Count how many times each var is used, so we can drop them when the counter reaches zero
	counter = asarray_create()
	asarray_notfound(counter, 0)
	// Variables passed through map_init_keepvars()
	keepvars = S.keepvars
	// ivars and cvars of the FEs
	for (g=1; g<=S.G; g++) {
		keepvars = keepvars, S.fes[g].ivars , S.fes[g].cvars
	}
	// Weight var
	if (S.weightvar!="") keepvars = keepvars, S.weightvar
	// Cluster vars
	for (h=1; h<=S.C; h++) {
		cl_ivars = tokens(S.clustervars[h], "#")
		cl_ivars = select(cl_ivars, cl_ivars:!="#")
		keepvars = keepvars, cl_ivars
	}
	// Time and panel vars
	if (S.vce_is_hac & S.timevar!="") keepvars = keepvars, S.timevar
	if (S.vce_is_hac & S.panelvar!="") keepvars = keepvars, S.panelvar

	// Fill
	for (i=1; i<=length(keepvars); i++) {
		asarray(counter, keepvars[i], asarray(counter, keepvars[i])+1)
	}

	// Report
	if (S.verbose>3) printf("{txt}{bf: 0. Usage count of each variable (dropped if it reaches zero)}\n")
	for (i=1; i<=asarray_elements(counter); i++) {
		if (S.verbose>3) printf("{txt}    - key=%s {col 30}count=%f\n", keepvars[i], asarray(counter,keepvars[i]))
	}

	// 1. Store permutation vectors and their invorder, generate ID variables, drop singletons
	if (S.verbose>0) printf("{txt}{bf: 1. Storing permutation vectors, generating ids, dropping singletons}\n")
	map_precompute_part1(S, counter)

	// 2. Store group offsets, group counters; demeaned(x), inv(xx) if num_slopes>0; weightvars
	if (S.verbose>0) printf("{txt}{bf: 2. Storing counters and offsets; processing cvars}\n")
	map_precompute_part2(S, counter)

	// 3. Create cluster IDs, report whether is/in/nested wrt cluster; store precomputed inv(p)
	if (S.verbose>0) printf("{txt}{bf: 3. Storing reverse permutation vectors, creating cluster IDs}\n")
	map_precompute_part3(S, counter)

	// 4. Keep only the essential variables
	keepvars  = J(1,0,"")
	for (loc=asarray_first(counter); loc!=NULL; loc=asarray_next(counter, loc)) {
		key = asarray_key(counter, loc)
		value = asarray_contents(counter, loc)
		if (value>0) keepvars = keepvars, key
	}
	if (S.verbose>3) printf("{txt}{bf: 4. Keeping the following variables}")
	st_keepvar(keepvars)
	// if (S.verbose>3) printf("\n{txt}    %s\n", invtokens(keepvars))
	if (S.verbose>3) stata("describe _all, numbers")
	if (S.verbose>3) printf("\n")

	// Store N (todo: add other details) to ensure the dataset doesn't change from now on
	S.N = st_nobs()
}
	
void map_precompute_part1(`Problem' S, transmorphic counter) {

	`Integer' G, i, j, n, g, h, i_last_singleton, num_singletons, initial_N
	`Boolean' sortedby
	`Group' id
	`Series' singleton, sum_singleton, inv_p
	`Varlist' idvarnames
	string scalar vartype
	pointer(`Series') scalar pp // Just to shorten code

	G = length(S.fes)
	i = i_last_singleton = g = 1

	// Give huge warning if keeping singletons
	if (S.keepsingletons) printf("{err}[WARNING] Singletons are not dropped; statistical significance will be biased\n")

	initial_N = st_nobs()

	// Loop until we stop discovering singletons (and then a bit more to be sure; G-1 to be exact)
	while (i<i_last_singleton+G) {
		if (g>G) g = 1
		if (S.verbose>1) printf("{txt}    - i=%f (g=%f/%f)\t(N=%f)\t", i, g, G, st_nobs())

		idvarnames = i<=G ? S.fes[g].ivars : S.fes[g].idvarname
		id = st_data(., idvarnames) // 2% of runtime
		if (i<=G) {
			for (j=1; j<=length(idvarnames); j++) {
				n = asarray(counter, idvarnames[j]) - 1
				asarray(counter, idvarnames[j], n)
				if (n==0) {
					st_dropvar(idvarnames[j])
				}
			}
		}
		if (i<=G) S.fes[g].is_sortedby = already_sorted(idvarnames)
		sortedby = S.fes[g].is_sortedby
		if (i<=G & !sortedby) {
			S.fes[g].p = order( id , 1..length(idvarnames) ) // 55% of function time
		}

		if (!sortedby) {
			_collate(id, S.fes[g].p) // sort id by p // 12% of function time
			inv_p = invorder(S.fes[g].p) // construct inv(p) that we'll use later
		}

		// Note that the lhs is actually the deltas, as in "bys id: gen delta = _n==1"
		id = rows_that_change(id) // 7% of function time
		singleton = S.keepsingletons ? J(rows(id), 1, 0) : select_singletons(id) // 5% of function time

		// Save IDs in dataset before dropping observations
		id = runningsum(id :* !singleton) // this is the ID now, not the deltas anymore
		S.fes[g].levels = id[length(id)]
		vartype = S.fes[g].levels<=100 ? "byte" : (S.fes[g].levels<=32740? "int" : "long")
		if (i<=G) {
			st_store(., st_addvar(vartype, S.fes[g].idvarname), sortedby? id : id[inv_p])
		}
		else {
			st_store(., S.fes[g].idvarname, sortedby? id : id[inv_p])
		}

		num_singletons = sum(singleton)
		if (num_singletons>0) {
			if (S.verbose>1) printf("{txt}(%f singletons)", num_singletons)
			i_last_singleton = i

			// Sort -singleton- as in the dataset, and use it to drop observations
			// 5% of function time
			singleton = sortedby? singleton : singleton[inv_p]
			st_dropobsif(singleton)
			if (!st_nobs()) {
				printf("{err}\nno observations left after dropping singletons\n")
				exit(error(2001))
			}

			// But now our precious sort orders (the p's) are useless! Fix them
			sum_singleton = runningsum(singleton)
			for (h=1;h<=G & h<=i; h++) { // 6% of function time
				if (S.fes[h].is_sortedby) continue
				pp = &(S.fes[h].p)
				(*pp) = select(*pp - sum_singleton[*pp] , !singleton[*pp] )
			}
		}

		if (S.verbose>1) printf("{txt}\n")
		i++
		g++
	}

	if (S.verbose>0) printf("{txt}    Total singleton obs. dropped: {res}%f{txt}\n", initial_N-st_nobs())

}

// -------------------------------------------------------------
// ALREADY_SORTED:
// -------------------------------------------------------------
`Integer' already_sorted(string vector vars) {
	`Varlist' sortedby
	sortedby = tokens(st_macroexpand("`" + ": sortedby" + "'"))
	return(length(vars) > length(sortedby) ? 0 : vars==sortedby[1..length(vars)])
}

// -------------------------------------------------------------
// ROWS_THAT_CHANGE: Return a 0/1 vector indicating what are diff from the previous row
// -------------------------------------------------------------
	// Idea: compromise between doing the operation in one go (uses lots of memory) vs doing loops (2x slower)
	// Other alternatives: loop row-by-row (slower), loop col-by-col and take max() (slower)
`Vector' rows_that_change(`Matrix' input) {
	`Vector' ans
	`Integer' i, j, K, N, stepsize
	
	// Size of blocks of matrices used (larger=faster smaller=less memory)
	// Benchmarks with 3 unsorted ivars showed 1e4 was fastest, followed by 1e3 and then 1e5
	stepsize = 1e4

	N = rows(input)
	K = cols(input)
	ans = J(N,1,0)
	ans[1] = 1
	for (i=2; i<=N;i=i+stepsize) {
		j = min((i+stepsize-1, N))
		ans[|i\j|] = rowmax(input[|i-1,1\j-1,K|] :!= input[|i,1\j,K|])
	}
	return(ans)
}

// -------------------------------------------------------------
// SELECT_SINGLETONS: 
// -------------------------------------------------------------
`Vector' select_singletons(`Vector' input) {
	// Code modified from <rows_that_change>
	`Vector' ans
	`Integer' i, j, N, stepsize

	// Size of blocks of matrices used (larger= hopefully faster, but smaller=less memory)
	// Benchmarks with 3 unsorted ivars showed 1e4 was fastest, followed by 1e3 and then 1e5
	stepsize = 1e4

	N = rows(input)
	ans = J(N,1,0)
	for (i=1; i<=N-1;i=i+stepsize) {
		j = min((i+stepsize-1, N-1))
		// We need that ans[i]==1 and ans[i+1]==1
		// Since ans is either 0 or 1, this is equivalent to
		ans[|i\j|] = (input[|i\j|] + input[|i+1\j+1|] :== 2) 
	}
	ans[N] = (input[N]==1) // Special case, last obs is singleton if it's the first obs in the group
	return(ans)
}
	
void map_precompute_part2(`Problem' S, transmorphic counter) {
	`Integer' G, g, k, K, n
	real scalar stdev
	`Boolean' sortedby
	`Series' id

	G = length(S.fes)
	if (S.weightvar!="") S.w = st_data(., S.weightvar)

	for (g=1; g<=G; g++) {
		sortedby = S.fes[g].is_sortedby
		K = S.fes[g].num_slopes
		assert(K==length(S.fes[g].cvars))
		if (S.verbose>1) printf("{txt}    - g=%f/%f\t\t(K=%f)\n", g, G, K)
		
		id = st_data(., S.fes[g].idvarname)
		if (!sortedby) id = id[S.fes[g].p]

		// Store offsets, counts (optionally weighted)
		S.fes[g].counts = count_by_group(id)
		S.fes[g].offsets = runningsum(S.fes[g].counts)
		if (S.weightvar!="") S.fes[g].counts = count_by_group(id, sortedby? S.w : S.w[S.fes[g].p])

		// Store cvars and related structures
		if (K>0) {

			// Store the cvars
			S.fes[g].x = st_data(., S.fes[g].cvars)

			// Drop cvars from dataset if not needed anymore
			for (k=1; k<=K; k++) {
				n = asarray(counter, S.fes[g].cvars[k]) - 1
				asarray(counter, S.fes[g].cvars[k], n)
				if (n==0) {
					st_dropvar(S.fes[g].cvars[k])
				}
			}

			// Sort the cvars if needed
			if (!sortedby) S.fes[g].x = S.fes[g].x[S.fes[g].p,]

			// Standardize
			// BUGBUG: Check that results don't change; specially on corner cases
			// EG: weights, no intercept, intercept, only one slope, etc.
			for (k=1; k<=K; k++) {
				// BUGBUG
				stdev = sqrt(quadvariance(S.fes[g].x[., k]))
				if (stdev>1e-5) S.fes[g].x[., k] = S.fes[g].x[., k] :/ stdev
			}

			// Demean X and precompute inv(X'X) (where X excludes the constant due to demeaning, if applicable)
			// Note that the demeaning is done directly in the S. structure
			S.fes[g].inv_xx = demean_and_compute_invxx(S, g)
		} // end of slope part
	
	} // next g
}

// -------------------------------------------------------------
// COUNT_BY_GROUP: alternative to mm_freq(group) (~10% of runtime)
// -------------------------------------------------------------
// Assume -id- (groups) and -w- (the weights) are sorted by id!
`Vector' function count_by_group(`Series' id, | `Series' w)
{
	`Integer' i, j, obs, levels, count
	`Boolean' has_weights
	`Vector' ans

	obs = rows(id)
	levels = id[length(id)]
	assert_msg(obs>levels, "Error: more levels of FE than observations!")
	has_weights = args()>1

	ans = J(levels, 1, 0)
	count = 0
	
	// <i> indexes observations, <j> indexes groups
	for (i=j=1; i<=obs; i++) {
		if (j<id[i]) {
			ans[j++] = count
			count = 0
		}
		count = count + (has_weights ? w[i] : 1) // optimize?
	}
	ans[j] = count // Last group
	assert( all(ans) ) // Counts *must* be positive for all levels
	return(ans)
}

// -------------------------------------------------------------
// DEMEAN_AND_COMPUTE_INVXX
// -------------------------------------------------------------
`Matrix' function demean_and_compute_invxx(`Problem' S, `Integer' g) {

	// j iterates over LEVELS; i iterates over OBS
	`Integer'	K, L, j, i_lower, i_upper
	`Boolean'	has_weights, sortedby, has_intercept
	`Matrix'	ans, tmp_x
	`Vector'	w, tmp_w
	real scalar	tmp_count
	K = S.fes[g].num_slopes // Exclude intercept
	L = S.fes[g].levels
	ans = J(L*K, K, 0)
	has_weights = S.weightvar !=""
	sortedby = S.fes[g].is_sortedby
	has_intercept = S.fes[g].has_intercept

	if (has_weights) {
		w = sortedby ? S.w : S.w[S.fes[g].p]
		assert(rows(w)==rows(S.fes[g].x))
	}
	
	i_lower = 1
	for (j=1; j<=L; j++) {
		i_upper = S.fes[g].offsets[j]
		tmp_count = S.fes[g].counts[j]
		tmp_w = has_weights ? w[| i_lower \ i_upper |] : 1
		tmp_x = S.fes[g].x[| i_lower , 1 \ i_upper , . |]
		if (has_intercept) {
			tmp_x = tmp_x :- (quadcolsum(has_weights ? tmp_x :* tmp_w : tmp_x) / tmp_count)
			S.fes[g].x[| i_lower , 1 \ i_upper , . |] = tmp_x
		}
		ans[| 1+(j-1)*K , 1 \ j*K , . |] = invsym(quadcross(tmp_x, tmp_w, tmp_x))
		i_lower = i_upper + 1

		// BUGBUG: quadcolsum???? quadcross????
		// use crossdev(x,means,w,x,means) if we don't demean beforehand
	}
	return(ans)
}
	
void map_precompute_part3(`Problem' S, transmorphic counter) {
	`Integer' g, h, i, j, n, L, i_lower, i_upper
	`Varname' var
	`Boolean' done, is_nested, sortedby
	`Vector' need_to_create_clustervar, range
	`Varlist' sorted_fe_ivars, sorted_cl_ivars, cl_ivars
	`String' vartype
	`Group' id
	`Series' p, sorted_cl_id

	need_to_create_clustervar = J(S.C, 1, 1)

	for (g=1;g<=S.G;g++) {
		S.fes[g].inv_p = invorder(S.fes[g].p)
		var = S.fes[g].idvarname
		st_varlabel(var, sprintf("[ID] %s", S.fes[g].varlabel))
		asarray(counter, var, asarray(counter, var)+1)

		done = 0
		sorted_fe_ivars = sort(S.fes[g].ivars', 1)'

		// 1. Check if the FE has the same ivars as a cluster (is_clustervar=1)
		for (h=1; h<=S.C;h++) {
			sorted_cl_ivars = tokens(S.clustervars[h], "#")
			sorted_cl_ivars = sort(select(sorted_cl_ivars, sorted_cl_ivars:!="#")', 1)'
			if (sorted_fe_ivars==sorted_cl_ivars) {
				need_to_create_clustervar[h] = 0
				S.clustervars[h] = var
				st_varlabel(var, sprintf("[CLUSTER] %s", st_varlabel(var)))
				S.fes[g].is_clustervar = 1
				done = 1
				break
			}
		}
	}

	// Create the cluster IDs if needed
	for (h=1; h<=S.C;h++) {
		cl_ivars = tokens(S.clustervars_original[h], "#")
		cl_ivars = select(cl_ivars, cl_ivars:!="#")
		
		
		for (j=1; j<=length(cl_ivars); j++) {
			n = asarray(counter, cl_ivars[j]) - 1
			assert_msg(n>=0, sprintf("counter[%s] was negative", cl_ivars[j]))
			asarray(counter, cl_ivars[j], n)
		}
		
		if (!need_to_create_clustervar[h]) continue
		if (cl_ivars==S.panelvar & S.vce_is_hac) continue
		if (cl_ivars==S.timevar  & S.vce_is_hac) continue

		id = st_data(., cl_ivars)

		// Construct and save cluster ID
		sortedby = already_sorted(cl_ivars)
		p = order( id , 1..length(cl_ivars) )
		if (!sortedby) {
			_collate(id, p) // sort id by p // 12% of function time
		}
		id = runningsum(rows_that_change(id))
		L = id[rows(id)]
		vartype = L<=100 ? "byte" : (L<=32740? "int" : "long")
		S.clustervars[h] = sprintf("__CL%f__", h)
		asarray(counter, S.clustervars[h], asarray(counter, S.clustervars[h])+1)
		st_store(., st_addvar(vartype, S.clustervars[h]), sortedby ? id : id[invorder(p)])
		st_varlabel(S.clustervars[h], sprintf("[CLUSTER] %s", S.clustervars_original[h]))
	}

	for (g=1;g<=S.G;g++) {
		var = S.fes[g].idvarname
		if (S.fes[g].is_clustervar) continue
		done = 0
		sorted_fe_ivars = sort(S.fes[g].ivars', 1)'

		// 2. Check if the FE ivars are a superset of those of the cluster (in_clustervar=1)
		for (h=1; h<=S.C;h++) {
			sorted_cl_ivars = tokens(S.clustervars_original[h], "#")
			sorted_cl_ivars = sort(select(sorted_cl_ivars, sorted_cl_ivars:!="#"), 1)
			if (length(sorted_cl_ivars)>=length(sorted_fe_ivars)) continue
			is_nested = 1
			for (i=1;i<=length(sorted_cl_ivars);i++) {
				if (!anyof(sorted_fe_ivars, sorted_cl_ivars[i])) {
					is_nested = 0
					break
				}
			}
			if (is_nested) {
				S.fes[g].in_clustervar = 1
				S.fes[g].nesting_clustervar = h
				done = 1
				break
			}
		}
		if (done) continue

		// 3. Check if the FE is nested within a cluster (e.g. cluster=state FE=zipcode)
		L = S.fes[g].levels
		for (h=1; h<=S.C; h++) {
			sorted_cl_id = st_data(., S.clustervars[h])
			if (!S.fes[g].is_sortedby) sorted_cl_id = sorted_cl_id[S.fes[g].p]
			i_lower = 1
			is_nested = 1
			for (j=1; j<=L; j++) {
				i_upper = S.fes[g].offsets[j]
				range = minmax(sorted_cl_id[| i_lower , 1 \ i_upper , . |])
				i_lower = i_upper + 1
				if (range[1]!=range[2]) {
					is_nested = 0
					break
				}
			}
			if (is_nested) {
				S.fes[g].in_clustervar = 1
				S.fes[g].nesting_clustervar = h
				break
			}
		}
	}

}
	
`Group' function map_projection(`Problem' S, `Integer' g, `Group' y) {
	`Integer' 	K, L, Q // Q is the number of depvars
	`Integer' 	j, i_lower, i_upper // j loops over levels, i loops over observations
	`Boolean' 	has_weights, sortedby, has_intercept, storing_betas
	`Series'	sorted_w
	`Group'		ans
	`Vector'	tmp_w, tmp_count
	real rowvector b
	real rowvector ymean // 1*Q
	real rowvector zero // 1*K
	`Matrix'	tmp_y, tmp_x
	pointer(`Series') scalar p_sorted_w
	pragma unset sorted_w // If we just set the pointer, what happens to the underlying data?

	// PROFILE TO SEE IF THIS HELPS OR NOT AT ALL
	//pointer(`Vector') scalar p_offset
	//p_offset = &(S.fes[g].offsets)

	has_weights = S.weightvar !=""
	sortedby = S.fes[g].is_sortedby
	has_intercept = S.fes[g].has_intercept
	K = S.fes[g].num_slopes
	Q = cols(y)
	L = S.fes[g].levels
	tmp_w = 1 // Harmless value for when there are no weights

	// Minimize copy+order operations on y
	if (has_weights) p_sorted_w = sortedby ? &(S.w) : &(sorted_w = S.w[S.fes[g].p, .])
	if (K>0) zero = J(1,K,0)

	ans = sortedby ? y : y[S.fes[g].p, .]

	i_lower = 1
	storing_betas = S.storing_betas & length(S.fes[g].target)>0
	for (j=1; j<=L; j++) {
		i_upper = S.fes[g].offsets[j]
		tmp_count = S.fes[g].counts[j]
		
		if (has_weights) tmp_w = (*p_sorted_w)[| i_lower \ i_upper |]
		tmp_y = ans[| i_lower , 1 \ i_upper , . |]
		// BUGBUG: quadcolsum or colsum ? Depends if there are dense FEs. Maybe condition it on L??
		if (has_weights) {
			ymean = has_intercept ? (colsum(tmp_y :* tmp_w) / tmp_count) : 0
		}
		else {
			ymean = has_intercept ? (colsum(tmp_y) / tmp_count) : 0
		}

		if (K>0) {
			tmp_x = S.fes[g].x[| i_lower , 1 \ i_upper , . |]
			// BUGBUG crossdev/cross or their quad version?
			if (has_intercept) {
				b = S.fes[g].inv_xx[| 1+(j-1)*K , 1 \ j*K , . |] * crossdev(tmp_x, zero, tmp_w, tmp_y, ymean)
			}
			else {
				b = S.fes[g].inv_xx[| 1+(j-1)*K , 1 \ j*K , . |] * cross(tmp_x, tmp_w, tmp_y)
			}
		}
		
		if (storing_betas) {
			if (has_intercept) S.fes[g].tmp_alphas[j, 1] = ymean
			if (K>0) S.fes[g].tmp_alphas[j, (has_intercept+1)..(has_intercept+K) ] = b'
		}

		// BUGBUG if we split this ternary will it be faster?
		ans[| i_lower , 1 \ i_upper , . |] = K>0 ? (ymean :+ tmp_x*b) : (ymean :+ J(i_upper-i_lower+1,Q,0))
		i_lower = i_upper + 1
	}
		
	return(sortedby ? ans : ans[S.fes[g].inv_p, .])
}
	
void function map_solve(`Problem' S, `Varlist' vars,
		| `Varlist' newvars, `Varlist' partial, `Boolean' save_fe) {
	`Integer' i, Q, Q_partial, offset, g
	`Group' y
	`FunctionPointer' transform, accelerate
	real rowvector stdevs
	`Varlist'	target

	if (S.verbose>0) printf("{txt}mata: map_solve()\n")
	assert_msg(S.N!=., "map_solve() needs to be run after map_precompute()")
	assert_msg(S.N==st_nobs(), "dataset cannot change after map_precompute()")

	// Load data
	// BUGBUG: This will use 2x memory for a while; partition the copy+drop based on S.groupsize?
	if (S.verbose>0) printf("{txt} - Loading variables into Mata\n")
	vars = tokens(vars)
	y = st_data(., vars)
	Q = cols(y)
	st_dropvar(vars) // We need the new ones on double precision

	if (args()<3 | newvars=="") newvars = vars
	assert_msg(length(vars)==length(newvars), "map_solve error: newvars must have the same size as vars")

	// Load additional partialled-out regressors
	Q_partial = 0
	if (args()==4 & partial!="") {
		vars = tokens(partial)
		Q_partial = cols(vars)
		y = y , st_data(., vars)
		st_dropvar(vars) // We need the new ones on double precision		
	}

	// Storing FEs and returning them requires 6 changes
	// 1) Extend the S and FE structures (add S.storing_betas, FE.alphas FE.tmp_alphas)
	// 2) Allocate them here
	// 3) Return results at the end of this function
	// 4) Within the accelerators, modify the SD to update the alphas
	// 5) Within map_projection, add a conditional to update tmp_alphas if needed
	S.storing_betas = 0
	if (args()<5) save_fe = 0
	assert_msg(save_fe==0 | save_fe==1, "map_solve error: save_fe must be either 0 or 1")
	if (save_fe) {
		assert_msg(partial=="", "map_solve error: partial must be empty if save_fe==1")
		assert_msg(length(vars)==1, "map_solve error: only one variable allowed if save_fe==1")
		if (S.verbose>0) printf("{txt} - Allocating objects to save the fixed effect estimates\n")
		S.storing_betas = 1
		for (g=1; g<=S.G; g++) {
			if (length(S.fes[g].target)>0) {
				S.fes[g].alphas = S.fes[g].tmp_alphas = 
					J(S.fes[g].levels, S.fes[g].has_intercept + S.fes[g].num_slopes, 0)
			}
		}
	}

	// Standardize all variables
	if (S.verbose>0) printf("{txt} - Standardizing variables\n")
	stdevs = J(1,cols(y),.)
	for (i=1; i<=cols(y); i++) {
		stdevs[i] = max((  sqrt(quadvariance(y[., i])) , sqrt(epsilon(1)) ))
	}
	y = y :/ stdevs

	// TODO: Report PEAK MEMORY that will be used
	// EG: de por si uso 2*size(y)
	// Luego bajo y asi que me quedo con y + ..
	// Contar cuantos vectores creo en los aceleradores y en los proyectores

	if (S.verbose>0) printf("{txt} - Solving problem (acceleration={res}%s{txt}, transform={res}%s{txt})\n", save_fe ? "steepest_descent" : S.acceleration, save_fe ? "kaczmarz" : S.transform)

	// Warnings
	if (S.transform=="kaczmarz" & S.acceleration=="conjugate_gradient") {
		printf("{err}(warning: convergence is {bf:unlikely} with transform=kaczmarz and accel=CG)\n")
	}

	// Load transform pointer
	if (S.transform=="cimmino") transform = &transform_cimmino()
	if (S.transform=="kaczmarz") transform = &transform_kaczmarz()
	if (S.transform=="symmetric_kaczmarz") transform = &transform_sym_kaczmarz()

	// Pointer to acceleration routine
	if (S.acceleration=="none") accelerate = &accelerate_none()
	if (S.acceleration=="conjugate_gradient") accelerate = &accelerate_cg()
	if (S.acceleration=="steepest_descent") accelerate = &accelerate_sd()
	if (S.acceleration=="aitken") accelerate = &accelerate_aitken()

	// Call acceleration routine
	if (save_fe) {
		y = accelerate_sd(S, y, &transform_kaczmarz()) // Only these were modified to save FEs
	}
	else if (S.groupsize>=cols(y)) {
		y = (*accelerate)(S, y, transform)
	}
	else {
		for (i=1;i<=cols(y);i=i+S.groupsize) {
			offset = min((i + S.groupsize - 1, cols(y)))
			y[., i..offset] = (*accelerate)(S, y[., i..offset], transform)
		}
	}

	// Partial-out variables
	if (Q_partial>0) {
		if (S.verbose>1) printf("{txt} - Partialling out variables\n")
		assert(cols(y)==Q+Q_partial)
		y = y[., 1..Q] - y[., (Q+1)..cols(y)] * qrsolve(y[., (Q+1)..cols(y)] , y[., 1..Q])
		stdevs =  stdevs[1..Q]
	}

	if (S.verbose>1) printf("{txt} - Saving transformed variables\n")
	// BUGBUG: This will use 2x memory for a while; partition the copy+drop based on S.groupsize?
	st_store(., st_addvar("double", newvars), y :* stdevs)

	// Store FEs
	if (save_fe) {
		if (S.verbose>1) printf("{txt} - Saving fixed effects\n")
		for (g=1; g<=S.G; g++) {
			target = S.fes[g].target
			if (length(target)>0) {
				S.fes[g].tmp_alphas = J(0,0,.)
				st_store(., st_addvar("double", target), S.fes[g].alphas[ st_data(., S.fes[g].idvarname) , . ] :* stdevs)
			}
		}
	}

}
	
// -------------------------------------------------------------------------------------------------
// Acceleration Schemes
// -------------------------------------------------------------------------------------------------

`Group' function accelerate_none(`Problem' S, `Group' y, `FunctionPointer' T) {
	`Integer'	iter
	`Group'		resid
	pragma unset resid

	for (iter=1; iter<=S.maxiterations; iter++) {
		(*T)(S, y, resid) // Faster version of "resid = S.T(y)"
		if (check_convergence(S, iter, resid, y)) break
		y = resid
	}
	return(resid)
}

// -------------------------------------------------------------------------------------------------
// Memory cost is approx = 4*size(y) (actually 3 since y is already there)
// But we need to add maybe 1 more due to u:*v
// And I also need to check how much does project and T use..
// Double check with a call to memory

// For discussion on the stopping criteria, see the following presentation:
// Arioli & Gratton, "Least-squares problems, normal equations, and stopping criteria for the conjugate gradient method". URL: https://www.stfc.ac.uk/SCD/resources/talks/Arioli-NAday2008.pdf

// Basically, we will use the Hestenes and Siefel rule

`Group' function accelerate_cg(`Problem' S, `Group' y, `FunctionPointer' T) {
	// BUGBUG iterate the first 6? without acceleration??
	`Integer'	iter, d, Q
	`Group'		r, u, v
	real rowvector alpha, beta, ssr, ssr_old, improvement_potential
	`Matrix' recent_ssr
	pragma unset r
	pragma unset v

	Q = cols(y)
	
	d = 2 // BUGBUG Set it to 2/3 // Number of recent SSR values to use for convergence criteria (lower=faster & riskier)
	// A discussion on the stopping criteria used is described in
	// http://scicomp.stackexchange.com/questions/582/stopping-criteria-for-iterative-linear-solvers-applied-to-nearly-singular-system/585#585

	improvement_potential = quadcolsum(y :* y)
	recent_ssr = J(d, Q, .)
	
	(*T)(S, y, r, 1)
	ssr = quadcolsum(r :* r) // cross(r,r) when cols(y)==1 // BUGBUG maybe diag(quadcross()) is faster?
	u = r

	for (iter=1; iter<=S.maxiterations; iter++) {
		(*T)(S, u, v, 1) // This is the hotest loop in the entire program
		// BUGBUG: colsum or quadcolsum??
		alpha = safe_divide( ssr , quadcolsum(u :* v) )
		recent_ssr[1 + mod(iter-1, d), .] = alpha :* ssr
		improvement_potential = improvement_potential - alpha :* ssr
		y = y - alpha :* u
		r = r - alpha :* v
		ssr_old = ssr
		ssr = quadcolsum(r :* r)
		beta = safe_divide( ssr , ssr_old) // Fletcher-Reeves formula, but it shouldn't matter in our problem
		u = r + beta :* u
		// Convergence if sum(recent_ssr) > tol^2 * improvement_potential
		if ( check_convergence(S, iter, sum(recent_ssr), improvement_potential, "hestenes") ) break
	}
	return(y)
}

// -------------------------------------------------------------------------------------------------

`Group' function accelerate_sd(`Problem' S, `Group' y, `FunctionPointer' T) {
	`Integer'	iter, g
	`Group' proj
	real rowvector t
	pragma unset proj

	for (iter=1; iter<=S.maxiterations; iter++) {
		(*T)(S, y, proj, 1)
		if (check_convergence(S, iter, y-proj, y)) break
		t = safe_divide( quadcolsum(y :* proj) , quadcolsum(proj :* proj) )
		// if (uniform(1,1)<0.1) t = 1 // BUGBUG: Does this help to randomly unstuck an iteration?
		y = y - t :* proj
		
		if (S.storing_betas) {
			for (g=1; g<=S.G; g++) {
				if (length(S.fes[g].target)>0) {
					S.fes[g].alphas = S.fes[g].alphas + t :* S.fes[g].tmp_alphas
				}
			}
		}
	}
	return(y-proj)
}

// -------------------------------------------------------------------------------------------------
// This is method 3 of Macleod (1986), a vector generalization of the Aitken-Steffensen method
// Also: "when numerically computing the sequence.. stop..  when rounding errors become too 
// important in the denominator, where the ^2 operation may cancel too many significant digits"
// Note: Sometimes the iteration gets "stuck"; can we unstuck it with adding randomness
// in the accelerate decision? There should be a better way.. (maybe symmetric kacz instead of standard one?)

`Group' function accelerate_aitken(`Problem' S, `Group' y, `FunctionPointer' T) {
	`Integer'	iter
	`Group'		resid, y_old, delta_sq
	`Boolean'	accelerate
	real rowvector t
	pragma unset resid

	//S.pause_length = 20
	//S.bad_loop_threshold = 1
	//S.stuck_threshold = 5e-3
	// old_error = oldest_error = bad_loop = acceleration_countdown = 0

	y_old = J(rows(y), cols(y), .)

	for (iter=1; iter<=S.maxiterations; iter++) {
		
		(*T)(S, y, resid)
		accelerate = iter>=S.accel_start & !mod(iter,S.accel_freq)

		// Accelerate
		if (accelerate) {
			delta_sq = resid - 2 * y + y_old // = (resid - y) - (y - y_old) // Equivalent to D2.resid
			// t is just (d'd2) / (d2'd2)
			t = safe_divide( quadcolsum( (resid - y) :* delta_sq) ,  quadcolsum(delta_sq :* delta_sq) )
			resid = resid - t :*  (resid - y)
		}

		// Only check converge on non-accelerated iterations
		// BUGBUG: Do we need to disable the check when accelerating?
		// if (check_convergence(S, iter, accelerate? resid :* .  : resid, y)) break
		if (check_convergence(S, iter, resid, y)) break
		
		// Experimental: Pause acceleration
		//if (accelerate) {
		//	improvement = max(( (old_error-update_error)/update_error , (oldest_error-update_error)/update_error ))
		//	bad_loop = improvement < stuck_threshold ? bad_loop+1 : 0
		//	// bad_loop, improvement, update_error, old_error, oldest_error
		//	// Tolerate two problems (i.e. 6+2=8 iters) and then try to unstuck
		//	if (bad_loop>bad_loop_threshold) {
		//		bad_loop = 0
		//		if (VERBOSE==3) printf(" Fixed point iteration seems stuck, acceleration paused\n")
		//		acceleration_countdown = pause_length
		//	}
		//	assert(bad_loop<=3)	
		//	oldest_error = old_error
		//	old_error = update_error
		//}
		//
		y_old = y // y_old is resid[iter-2]
		y = resid // y is resid[iter-1]
	}
	return(resid)
}

// -------------------------------------------------------------------------------------------------

`Boolean' check_convergence(`Problem' S, `Integer' iter, `Group' y_new, `Group' y_old,| `String' method) {
	`Boolean'	done, is_last_iter
	`Real'		update_error

	// max() ensures that the result when bunching vars is at least as good as when not bunching
	if (args()<5) method = "vectors" 

	if (method=="vectors") {
		update_error = max(mean(reldif(y_new, y_old)))
	}
	else if (method=="hestenes") {
		// If the regressor is perfectly explained by the absvars, we can have SSR very close to zero but negative
		// (so sqrt is missing)
		update_error = max(safe_divide( sqrt(y_new) , editmissing(sqrt(y_old), sqrt(epsilon(1)) ) , sqrt(epsilon(1)) ))
	}
	else {
		exit(error(100))
	}

	done = update_error <= S.tolerance
	is_last_iter = iter==S.maxiterations
	
	if (done) {
		if (S.verbose>1) printf("\n{txt} - Converged in %g iterations (last error =%3.1e)\n", iter, update_error)
		if (S.verbose==1) printf("{txt} converged in %g iterations last error =%3.1e)\n", iter, update_error)
	}
	else if (is_last_iter) {
		printf("\n{err}convergence not achieved in %g iterations (last error=%e); try increasing maxiter() or decreasing tol().\n", S.maxiterations, update_error)
		exit(430)
	}
	else {
		if ((S.verbose>=2 & S.verbose<=3 & mod(iter,1)==0) | (S.verbose==1 & mod(iter,100)==0)) {
			printf("{txt}.")
			displayflush()
		}
		if (S.verbose>=2 & S.verbose<=3 & mod(iter,100)==0) printf("{txt}%9.1f\n", update_error/S.tolerance)
		if (S.verbose==4) printf("{txt} iter={res}%4.0f{txt}\tupdate_error={res}%-9.6e\n", iter, update_error)
	}
	return(done)
}
	
// -------------------------------------------------------------------------------------------------
// Transformations: Compute RESIDUALS, not projections
// -------------------------------------------------------------------------------------------------

void function transform_cimmino(`Problem' S, `Group' y, `Group' ans,| `Boolean' get_proj) {
	`Integer' 	g, G
	G = S.G
	if (args()<4) get_proj = 0

	ans = map_projection(S, 1, y)
	for (g=2; g<=G; g++) {
		ans = ans + map_projection(S, g, y)
	}
	ans = get_proj ? ans / G : y - ans / G
}

// -------------------------------------------------------------------------------------------------

void function transform_kaczmarz(`Problem' S, `Group' y, `Group' ans,| `Boolean' get_proj) {
	`Integer' 	g, G
	G = S.G
	if (args()<4) get_proj = 0

	ans = y - map_projection(S, 1, y)
	for (g=2; g<=G; g++) {
		ans = ans - map_projection(S, g, ans)
	}
	if (get_proj) ans = y - ans
}

// -------------------------------------------------------------------------------------------------
// This seems slower than plain kaczmarz; not used currently
void function transform_rand_kaczmarz(`Problem' S, `Group' y, `Group' ans,| `Boolean' get_proj) {
	`Integer' 	g, G
	G = S.G
	`Vector' rand
	if (args()<4) get_proj = 0
	rand = sort( ( (1::G) , uniform(G,1) ) , 2 )[.,1]

	ans = y - map_projection(S, rand[1], y)
	for (g=2; g<=G; g++) {
		ans = ans - map_projection(S, rand[g], ans)
	}
	if (get_proj) ans = y - ans
}

// -------------------------------------------------------------------------------------------------

 void function transform_sym_kaczmarz(`Problem' S, `Group' y, `Group' ans,| `Boolean' get_proj) {
	`Integer' 	g, G
	// BUGBUG: Streamline and remove all those "ans - .." lines?
	G = S.G
	if (args()<4) get_proj = 0

	ans = y - map_projection(S, 1, y)
	for (g=2; g<=G; g++) {
		ans = ans - map_projection(S, g, ans)
	}
	for (g=G-1; g>=1; g--) {
		ans = ans - map_projection(S, g, ans)
	}
	if (get_proj) ans = y - ans
}
	
void map_estimate_dof(`Problem' S, string rowvector adjustments, 
		| `Varname' groupvar) {
	`Boolean' adj_firstpairs, adj_pairwise, adj_clusters, adj_continuous, belongs, already_first_constant
	string rowvector all_adjustments
	`String' adj, label, basestring
	`Integer' i, g, SuperG, SubGs, h, M_due_to_nested, num_groups, j, m, sum_levels
	`Vector' M, M_is_exact, M_is_nested, is_slope, solved, prev_g

	// Parse list of adjustments/tricks to do
	if (S.verbose>1) printf("\n")
	if (S.verbose>0) printf("{txt}{bf:mata: map_estimate_dof()}\n")
	if (S.verbose>1) printf("{txt} - Estimating degrees-of-freedom used by the fixed effects\n")
	all_adjustments = "firstpairs", "pairwise", "clusters", "continuous"
	adjustments = tokens(adjustments)
	for (i=1; i<=length(adjustments);i++) {
		assert_msg(anyof(all_adjustments, adjustments[i]), 
			sprintf("map_estimate_dof error: adjustment %s invalid", adjustments[i]))
	}
	if (S.verbose>1) printf("{txt} - Adjustments:\n")
	for (i=1;i<=length(all_adjustments);i++) {
		adj = all_adjustments[i]
		belongs = anyof(adjustments, adj)
		if (S.verbose>1) printf("{txt}    - %s:  {col 20}{res} %s\n", adj, belongs ? "yes" : "no")
		if (adj=="firstpairs") adj_firstpairs = belongs
		if (adj=="pairwise") adj_pairwise = belongs
		if (adj=="clusters") adj_clusters = belongs
		if (adj=="continuous") adj_continuous = belongs
	}

	// Assert that the clustervars exist
	for (i=1;i<=S.C;i++) {
		stata(sprintf("confirm numeric variable %s, exact", S.clustervars[i]))
	}

	// Can only save connected group if firstpairs or pairwise are active
	if (args()<3) groupvar = ""
	if (groupvar!="") {
		assert_msg(adj_firstpairs | adj_pairwise, "map_estimate_dof error: group option requires 'pairwise' or 'firstpairs' adjustments")
	}

	// Count all fixed intercepts and slopes
	SubGs = J(S.G, 1, 0) // Intercept + # of slopes in an absvar
	for (g=1;g<=S.G;g++) {
		SubGs[g] = S.fes[g].has_intercept + S.fes[g].num_slopes
	}
	SuperG = sum(SubGs)
	if (S.verbose>1) printf("{txt} - There are %f fixed intercepts and slopes in the %f absvars\n", SuperG, S.G)

	// Initialize result vectors and scalars
	M = J(SuperG, 1, 1)
	M_is_exact = J(SuperG, 1, 0)
	M_is_nested = J(SuperG, 1, 0)
	is_slope = J(SuperG, 1, .)
	solved = J(SuperG, 1, 0)

	// Initial Fill
	h = 0
	already_first_constant = 0
	for (g=1;g<=S.G;g++) {
		for (i=1;i<=SubGs[g];i++) {
			h++
			if (is_slope[h] = i>S.fes[g].has_intercept) {
				M[h] = 0
			}
			else if (!already_first_constant) {
				already_first_constant = 1
				M[h] = 0
			}

		}
	}

	// (Intercept-Only) Look for absvars that are clustervars or are nested within a clustervar
	h = 1
	M_due_to_nested = 0
	if (adj_clusters) {
		for (g=1;g<=S.G;g++) {
			if (S.fes[g].has_intercept & (S.fes[g].is_clustervar | S.fes[g].in_clustervar)) {
				M[h] = S.fes[g].levels
				M_is_exact[h] = M_is_nested[h] = 1
				M_due_to_nested = M_due_to_nested + M[h]
				solved[h] = 1
				if (S.verbose>1 & S.fes[g].is_clustervar) printf("   {txt}(categorical variable {res}%s{txt} is also a cluster variable, so it doesn't count towards DoF)\n", invtokens(S.fes[g].ivars,"#"))
				if (S.verbose>1 & S.fes[g].in_clustervar) printf("   {txt}(categorical variable {res}%s{txt} is nested within cluster {res}%s{txt}, so it doesn't count towards DoF)\n", invtokens(S.fes[g].ivars,"#"), S.clustervars_original[S.fes[g].nesting_clustervar])
			}
			h = h + SubGs[g]
		}
	}

	// (Intercept-only) Excluding those already solved, the first absvar is exact, and the second can be with pairwise/firstpairs
	i = 0
	h = 1
	prev_g = J(S.G, 1, 0)
	for (g=1;g<=S.G;g++) {
		if (!solved[h] & S.fes[g].has_intercept) {
			i++
			if (i==1) {
				M_is_exact[h] = 1
			}
			else if (i==2 & (adj_pairwise | adj_firstpairs)) {
				M_is_exact[h] = 1
				m = map_connected_groups(S, prev_g[1], g, groupvar)
				if (S.verbose>2) printf("{txt}    - Mobility groups between fixed intercept #%f and #%f: {res}%f\n", prev_g[1], g, m)
				M[h] = m
			}
			else if (i>2 & adj_pairwise) {
				// Call connected in a LOOP (but I need to save the list of those that I needed to check)
				for (j=1; j<i; j++) {
					m = map_connected_groups(S, prev_g[j], g)
					if (S.verbose>2) printf("{txt}    - Mobility groups between fixed intercept #%f and #%f: {res}%f\n", prev_g[j], g, m)
					M[h] = max((M[h], m))
				}
				if (S.verbose>2) printf("{txt}    - Maximum of mobility groups wrt fixed intercept #%f: {res}%f\n", g, M[h])

			}
			prev_g[i] = g
		}
		h = h + SubGs[g]
	}

	// See if cvars are zero (w/out intercept) or just constant (w/intercept)
	if (adj_continuous) {
		h = 0
		for (g=1;g<=S.G;g++) {
			for (i=1;i<=SubGs[g];i++) {
				h++
				// If model has intercept, redundant cvars are those that are CONSTANT
				// Without intercept, a cvar has to be zero within a FE for it to be redundant
				// Since S.fes[g].x are already demeaned IF they have intercept, we don't have to worry about the two cases
				if (is_slope[h]) {
					M[h] = count_redundant_cvars(S, g, i)
				}
			}
		}
	}

	// Report and return final results
	st_rclear()
	st_numscalar("r(M)", sum(M))
	st_numscalar("r(M_due_to_nested)", M_due_to_nested)
	sum_levels = 0
	for (g=1;g<=S.G;g++) sum_levels = sum_levels + S.fes[g].levels * (S.fes[g].has_intercept + S.fes[g].num_slopes)
	st_numscalar("r(df_a)", sum_levels - sum(M))

	h = 0
	if (S.verbose>=2) printf("{txt} - Degrees-of-freedom used by each fixed effect (K=total levels; M=redundant levels)\n")
	for (g=1;g<=S.G;g++) {
		for (i=1;i<=SubGs[g];i++) {
			h++
			st_numscalar(sprintf("r(M%f)",h), M[h])
			st_numscalar(sprintf("r(M%f_exact)",h), M_is_exact[h])
			st_numscalar(sprintf("r(M%f_nested)",h), M_is_nested[h])
			st_numscalar(sprintf("r(K%f)",h), S.fes[g].levels)
			if (S.verbose>=2) {
				label = invtokens(S.fes[g].ivars, "#")
				if (i>S.fes[g].has_intercept) label = label + "#c." + S.fes[g].cvars[i-S.fes[g].has_intercept]
				basestring = "{txt}   - FE%f ({res}%s{txt}): {col 40}K=%f {col 50}M=%f {col 60}is_exact=%f\n"
				printf(basestring, g, label, S.fes[g].levels, M[h], M_is_exact[h])
			}
		}
	}

	if (S.verbose>0) printf(" - Results: N=%f ; K=%f ; M=%f ; (K-M)==df_a=%f\n", S.N, sum_levels, sum(M), sum_levels-sum(M))
}

// -------------------------------------------------------------------------------------------------

`Integer' function count_redundant_cvars(`Problem' S, `Integer' g, `Integer' i) {
	`Integer' j, i_lower, i_upper, ans, L, ii
	real rowvector min_max
	`Series' x

	ii = i-S.fes[g].has_intercept
	ans = 0
	L = S.fes[g].levels
	x = S.fes[g].x[., ii]
	
	i_lower = 1
	for (j=1;j<=L; j++) {
		i_upper = S.fes[g].offsets[j]
		min_max = minmax(x[| i_lower \ i_upper |])
		if (sqrt(epsilon(1))>abs(min_max[1]) & sqrt(epsilon(1))>abs(min_max[2])) ans++
		i_lower = i_upper + 1
	}
	if (S.verbose>=2) printf("{txt}    - Fixed slope {res}%s#c.%s {txt}has {res}%f/%f{txt} redundant coefs.\n", invtokens(S.fes[g].ivars,"#"), S.fes[g].cvars[ii], ans, L)
	return(ans)
}

/*
	In general, we can't know the exact number of DoF lost because we don't know when multiple FEs are collinear
	When we have two pure FEs, we can use an existing algorithm, but besides that we'll just use an upper (conservative) bound

	Features:
	 - Save the first mobility group if asked
	 - Within the pure FEs, we can use the existing algorithm pairwise (FE1 vs FE2, FE3, .., FE2 vs FE3, ..)
	 - If there are n pure FEs, that means the algo gets called n! times, which may be kinda slow
	 - With FEs interacted with continuous variables, we can't do this, but can do two things:
		a) With i.a#c.b , whenever b==0 for all values of a group (of -a-), add one redundant
		b) With i.a##c.b, do the same but whenever b==CONSTANT (so not just zero)
     - With clusters, it gets trickier but in summary you don't need to penalize DoF for params that only exist within a cluster. This happens:
		a) if absvar==clustervar
		b) if absvar is nested within a clustervar. EG: if we do vce(cluster state), and -absorb(district)- or -absorb(state#year)
		c) With cont. interactions, e.g. absorb(i.state##c.year) vce(cluster year), then i) state FE is redundant, but ii) also state#c.year
		   The reason is that at the param for each "fixed slope" is shared only within a state

	Procedure:
	 - Go through all FEs and see if i) they share the same ivars as any clusters, and if not, ii) if they are nested within clusters
	 - For each pure FE in the list, run the algorithm pairwise, BUT DO NOT RUN IT BEETWEEN TWO PAIRS OF redundant
	   (since the redundants are on the left, we just need to check the rightmost FE for whether it was tagged)
	 - For the ones with cont interactions, do either of the two tests depending on the case

	Misc:
	 - There are two places where DoFs enter in the results:
		a) When computing e(V), we do a small sample adjustment (seen in Stata documentation as the -q-)
		   Instead of doing V*q with q = N/(N-k), we use q = N / (N-k-kk), so THE PURPOSE OF THIS PROGRAM IS TO COMPUTE "kk"
		   This kk will be used to adjust V and also stored in e(df_a)
		   With clusters, q = (N-1) / (N-k-kk) * M / (M-1)
		   With multiway clustering, we use the smallest N_clust as our M
	    b) In the DoF of the F and t tests (not when doing chi/normal)
	       When there are clusters, note that e(df_r) is M-1 instead of N-1-k
	       Again, here we want to use the smallest M with multiway clustering

	Inputs: +-+- if we just use -fe2local- we can avoid passing stuff around when building subroutines
	 - We need the current name of the absvars and clustervars (remember a#b is replaced by something different)
	 - Do a conf var at this point to be SURE that we didn't mess up before
	 - We need the ivars and cvars in a list
	 - For the c. interactions, we need to know if they are bivariate or univariate
	 - SOLN -> mata: fe2local(`g')  ; from mata: ivars_clustervar`i' (needed???) , and G
	 - Thus, do we really needed the syntax part??
	 - fe2local saves: ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels // Z group_k weightvar

	DOF Syntax:
	 DOFadjustments(none | all | CLUSTERs | PAIRwise | FIRSTpair | CONTinuous)
	 dof() = dof(all) = dof(cluster pairwise continuous)
	 dof(none) -> do nothing; all Ms = 0 
	 dof(first) dof(first cluster) dof(cluster) dof(continuous)

	For this to work, the program MUST be modular
*/
	
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
// -------------------------------------------------------------------------------------------------

program define reghdfe

* Set Stata version
	version `=clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred

* Intercept version calls
	cap syntax, version
	if !c(rc) {
		Version
		exit
	}

* Intercept replays
	if replay() {
		if (`"`e(cmd)'"'!="reghdfe") error 301
		Replay `0'
		exit
	}

* Finally, call Estimate
	cap noi Estimate `0'
	if (c(rc)) {
		local rc = c(rc)
		cap mata: mata drop HDFE_S // overwrites c(rc)
		exit `rc'
	}
end

// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------
// Simple assertions
// -------------------------------------------------------------

program define Assert
    syntax anything(everything equalok) [, MSG(string asis) RC(integer 198)]
    if !(`anything') {
        di as error `msg'
        exit `rc'
    }
end


// -------------------------------------------------------------
// Simple debugging
// -------------------------------------------------------------

program define Debug

	syntax, [MSG(string asis) Level(integer 1) NEWline COLOR(string)] [tic(integer 0) toc(integer 0)]
	
	mata: verbose2local(HDFE_S, "VERBOSE")
	assert "`VERBOSE'"!=""
	assert inrange(`VERBOSE',0, 5)
	
	assert inrange(`level',0, 5)
	assert (`tic'>0) + (`toc'>0)<=1

	if ("`color'"=="") local color text
	assert inlist("`color'", "text", "res", "result", "error", "input")

	if (`VERBOSE'>=`level') {

		if (`tic'>0) {
			timer clear `tic'
			timer on `tic'
		}
		if (`toc'>0) {
			timer off `toc'
			qui timer list `toc'
			local time = r(t`toc')
			if (`time'<10) local time = string(`time'*1000, "%tcss.ss!s")
			else if (`time'<60) local time = string(`time'*1000, "%tcss!s")
			else if (`time'<3600) local time = string(`time'*1000, "%tc+mm!m! SS!s")
			else if (`time'<24*3600) local time = string(`time'*1000, "%tc+hH!h! mm!m! SS!s")
			timer clear `toc'
			local time `" as result " `time'""'
		}

		if (`"`msg'"'!="") di as `color' `msg'`time'
		if ("`newline'"!="") di
	}
end


// -------------------------------------------------------------
// Report HDFE/REGHDFE version
// -------------------------------------------------------------

program define Version, eclass
    local version "3.0.809 07may2015"
    ereturn clear
    di as text "`version'"
    ereturn local version "`version'"

    di as text _n "Dependencies installed?"
    local dependencies ivreg2 avar tuples parallel
    foreach dependency of local dependencies {
    	cap findfile `dependency'.ado
    	if (_rc) {
    		di as text "{lalign 20:- `dependency'}" as error "not"
    	}
    	else {
    		di as text "{lalign 20:- `dependency'}" as result "yes"
    	}
    }

end

// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// Transform data and run the regression
// -------------------------------------------------------------------------------------------------

program define Estimate, eclass

* (reporting) Store dataset size
	qui de, simple
	local old_mem = string(r(width) * r(N)  / 2^20, "%6.2f")
	local RAW_N = c(N)
	local RAW_K = c(k)

* Parse arguments and create the HDFE_S Mata structure
	Parse `0' // save all arguments into locals (verbose>=3 shows them)
	mata: verbose2local(HDFE_S, "VERBOSE")

* (optional) Preserve
	if (!`clear') {
		preserve
		Debug, level(2) newline
		Debug, level(2) msg("(dataset preserved)")
	}

* (optional) Create uid so we can then attach e(sample) and/or the Zs (the FE coefs.)
	if (!`clear' & !`fast') {
		tempvar uid
		local uid_type = cond(c(N)>c(maxlong), "double", "long")
		gen `uid_type' `uid' = _n // Useful for later merges
		la var `uid' "[UID]"
	}

* Drop unused variables
	local exp "= `weightvar'"
	if ("`weightvar'"!="") la var `weightvar' "[WEIGHT] `: var label `weightvar''" // so we can distinguish it with -describe-
	marksample touse, novar // Uses -if- , -in- ; -weight-? and -exp- ; can't drop any var until this
	keep `uid' `touse' `basevars' `timevar' `panelvar' `weightvar' `absorb_keepvars' `clustervars' `over'

* Expand factor and time-series variables
	local expandedvars
	local sets depvar indepvars endogvars instruments // depvar MUST be first
	Debug, level(4) newline
	Debug, level(4) msg("{title:Expanding factor and time-series variables:}")
	foreach set of local sets {
		local varlist ``set''
		if ("`varlist'"=="") continue
		local original_`set' `varlist'
		* the -if- prevents creating dummies for categories that have been excluded
		ExpandFactorVariables `varlist' if `touse', setname(`set') verbose(`VERBOSE')
		local `set' "`r(varlist)'"
		local expandedvars `expandedvars' ``set''
	}

* Drop unused basevars and tsset vars (usually no longer needed)
	if ("`vceextra'"!="") local tsvars `panelvar' `timevar' // We need to keep them only with autoco-robust VCE
	keep `uid' `touse' `expandedvars' `weightvar' `absorb_keepvars' `clustervars' `tsvars' `over' 

* Drop excluded observations and observations with missing values
	markout `touse' `expandedvars' `weightvar' `absorb_keepvars' `clustervars'
	qui keep if `touse'
	drop `touse'
	if ("`weightvar'"!="") qui drop if (`weightvar'==0)
	Assert c(N)>0, rc(2000) msg("Empty sample, check for missing values or an always-false if statement")

* Precompute Mata objects
	mata: map_init_keepvars(HDFE_S, "`expandedvars' `uid'") // Non-essential vars will be deleted
	* Note: This is kinda redundant with the -keep- above except if a variable is somehow a cvar
	mata: map_precompute(HDFE_S)

* (reporting) memory usage demeanings
	Debug, level(2) msg("(dataset compacted: observations " as result "`RAW_N' -> `c(N)'" as text " ; variables " as result "`RAW_K' -> `c(k)'" as text ")")
	qui de, simple
	local new_mem = string(r(width) * r(N) / 2^20, "%6.2f")
	Debug, level(2) msg("(dataset compacted, c(memory): " as result "`old_mem'" as text "M -> " as result "`new_mem'" as text "M)")
	if (`VERBOSE'>3) {
		di as text "(memory usage including mata:)"
		memory
		di as text ""
	}

* Save the statistics we need before transforming the variables
	* Compute TSS of untransformed depvar
	local tmpweightexp = subinstr("`weightexp'", "[pweight=", "[aweight=", 1)
	qui su `depvar' `tmpweightexp' // BUGBUG: Is this correct?!
	local tss = r(Var)*(r(N)-1)
	assert `tss'<.
	if (`: list posof "first" in stages') {
		foreach var of varlist `endogvars' {
			qui su `var' `tmpweightexp' // BUGBUG: Is this correct?!
			local tss_`var' = r(Var)*(r(N)-1)
		}
	}

* Compute e(df_a)
	mata: map_estimate_dof(HDFE_S, "`dofadjustments'", "`groupvar'")
	TODO: SAVE GROUPVAR IN HDFE_S AND WAIT UNTIL RESTORE TO PUT IT BACK IN THE DTA
	local M = r(M) // FEs found to be redundant
	local M_due_to_nested = r(M_due_to_nested)
	local kk = r(df_a) // FEs that were not found to be redundant (total FEs - redundant FEs)
	Assert `kk'<.
	Assert `M'>=0 & `M'<.

	forv g=1/`N_hdfe' {
		local M`g' = r(M`g')
		local K`g' = r(K`g')
		local M`g'_exact = r(M`g'_exact)
		local M`g'_nested = r(M`g'_nested)

		assert inlist(`M`g'_exact',0,1) // 1 or 0 whether M`g' was calculated exactly or not
		assert `M`g''<. & `K`g''<.
		assert `M`g''>=0 & `K`g''>=0
		assert inlist(r(drop`g'), 0, 1)
	}

* Drop IDs for the absorbed FEs (except if its the clustervar)
* This is useful because the demeaning takes a lot of memory
	TODO: Make this a MATA call that checks whether fes[g].is_clustervar and also fes[g].target

* Replace vceoption with the correct cluster names (e.g. if it's a FE or a new variable)
	TODO
	//if (`num_clusters'>0) {
	//	mata: st_local("temp_clustervars", invtokens(clustervars))
	//	local vceoption : subinstr local vceoption "<CLUSTERVARS>" "`temp_clustervars'"
	//}





asdasdasdasd



* 12) Save untransformed data.
*	This allows us to:
*	i) do nested ftests for the FEs,
*	ii) recover the FEs, compute their correlations with xb, check that FE==1

	* We can avoid this if i) nested=check=0 ii) targets={} iii) fast=1
	mata: st_local("any_target_avge", strofreal(any(avge_target :!= "")) ) // saving avge?
	local any_target_hdfe 0 // saving hdfe?
	forv g=1/`N_hdfe' {
		mata: fe2local(`g')
		if (!`is_bivariate' | `is_mock') local hdfe_cvar`g' `cvars'
		// If it's the intercept part of the bivariate absorbed effect, don't add the cvar!
		local hdfe_target`g' `target'
		if ("`target'"!="") local any_target_hdfe 1
	}

	if (`fast') {
		if (`nested' | `check' | `any_target_hdfe' | `any_target_avge' | "`group'"!="") {
			Debug, msg(as text "(option {it:fast} not compatible with other options; disabled)") level(0)
			local fast 0
		}
		else {
			Debug, msg("(option {opt fast} specified; will not save e(sample) or compute correlations)")
		}
	}

	if (!`fast' | `cores'>1) {
		sort `uid'
		tempfile original_vars
		qui save "`original_vars'"
		if (`cores'>1) local parallel_opt `" filename("`original_vars'") uid(`uid') cores(`cores') "'
		Debug, msg("(untransformed dataset saved)") level(2)
	}

* 13) (optional) Compute R2/RSS to run nested Ftests on the FEs
	* a) Compute R2 of regression without FE, to build the joint FTest for all the FEs
	* b) Also, compute RSS of regressions with less FEs so we can run nested FTests on the FEs
	if ("`model'"=="ols" & !`savingcache') {
		qui _regress `vars' `weightexp', noheader notable
		local r2c = e(r2)

		if (`nested') {
			local rss0 = e(rss)
			local subZs
			forv g=1/`=`N_hdfe'-1' {
				Debug, msg("(computing nested model w/`g' FEs)")
				if (`cores'>1) {
					DemeanParallel, varlist(`vars') `maximize_options' num_fe(`g') self(reghdfe) `parallel_opt'
				}
				else {
					Demean, varlist(`vars') `maximize_options' num_fe(`g')	
				}

				qui _regress `vars' `weightexp', noheader notable
				local rss`g' = e(rss)
				qui use "`original_vars'", clear // Back to untransformed dataset
			}
		}
	}

	* Get normalized string of the absvars (i.e. turn -> i.turn)
	local original_absvars
	forv g=1/`N_hdfe' {
		mata: fe2local(`g')
		local original_absvars `original_absvars'  `varlabel'
	}

* Compute summary statistics for the all the regression variables
	if ("`stats'"!="") {
		local tabstat_weight : subinstr local weightexp "[pweight" "[aweight"
		qui tabstat `vars' `tabstat_weight' , stat(`stats') col(stat) save
		tempname statsmatrix
		matrix `statsmatrix' = r(StatTotal)
	}

* 14) Compute residuals for all variables including the AvgEs (overwrites vars!)
	qui ds `vars'
	local NUM_VARS : word count `r(varlist)'
	Debug, msg("(computing residuals for `NUM_VARS' variables)")
	Debug, msg(" - tolerance = `tolerance'")
	Debug, msg(" - max. iter = `maxiterations'")
	if ("`usecache'"=="") {
		if (`cores'>1) {
			DemeanParallel, varlist(`vars') `maximize_options' self(reghdfe) `parallel_opt'
		}
		else {
			Demean, varlist(`vars') `maximize_options'	
		}
	}
	else {
		Debug, msg("(using cache data)")
		drop `vars'
		local handshake_master : char __uid__[handshake]
		char __uid__[handshake]
		// An error in the merge most likely means different # of obs due to missing values in a group but not in other
		// try with if !missing(__uid__) // TODO: Auto-add this by default?
		// TODO: Make this fool-proof when using -over-
		if ("`over'"!="") local using using // This is dangerous
		sort __uid__ // The user may have changed the sort order of the master data
		qui merge 1:1 __uid__ using "`usecache'", keepusing(`vars') assert(match master `using') keep(master match) nolabel sorted
		qui cou if _merge!=3
		Assert r(N)==0, msg(as error "Error: the cache has `r(N)' less observations than the master data" _n ///
			as text " - This is possibly because, when created, it included variables that were missing in cases where the current ones are not.")
		qui drop if _merge!=3
		drop _merge

		local handshake_using : char __uid__[handshake]
		local tolerance_using : char __uid__[tolerance]
		local maxiterations_using : char __uid__[maxiterations]
		Assert (`handshake_master'==`handshake_using'), msg("using dataset does not have the same __uid__")
		Assert abs(`tolerance'-`tolerance_using')<epsdouble(), msg("using dataset not computed with the same tolerance (`tolerance_using')")
		Assert (`maxiterations'==`maxiterations_using'), msg("using dataset not computed with the same maxiterations (`maxiterations_using')")

		local absvar_master `original_absvars'
		local absvar_using : char __uid__[absvars_key]
		Assert ("`absvar_master'"=="`absvar_using'"), msg("using dataset not created with the same absvars")
		char __uid__[absvars_key]
	}

if (`savingcache') {
	Debug, msg("(saving cache and exiting)")
	char __uid__[absvars_key] `original_absvars'
	sort __uid__
	save "`savecache'", replace
	return clear
	ereturn clear
	ereturn local cmdline `"`cmdline'"'
	if ("`over_levels'"!="") ereturn local over_levels = "`over_levels'"
	exit
}

// PART II - REGRESSION

**** <<< START OF UGLY -stages- CODE
assert "`stages'"!=""
if ("`stages'"!="none") {
	Debug, level(2) msg(_n " {title:Stages to run}: " as result "`stages'" _n)
	local backup_fast `fast'
	local num_stages : word count `stages'
	local last_stage : word `num_stages' of `stages'
	assert "`last_stage'"=="iv"
	foreach vargroup in depvar indepvars endogvars instruments {
		local backup_`vargroup' ``vargroup''
		local backup_original_`vargroup' `original_`vargroup''
	}
	local backup_tss = `tss'
}

foreach stage of local stages {
local lhs_endogvars = cond("`stage'"=="first", "`backup_endogvars'", "<none>")

if ("`stage'"=="first") {
	local i_endogvar 0
}
else {
	local i_endogvar
}

foreach lhs_endogvar of local lhs_endogvars {
Assert inlist("`stage'", "none", "iv", "first", "ols", "reduced", "acid")

if ("`stage'"=="iv") {
	local tss = `backup_tss'
	local fast `backup_fast'
	local depvar `backup_depvar'
	local indepvars `backup_indepvars'
	local endogvars `backup_endogvars'
	local instruments `backup_instruments'
	local original_depvar `backup_original_depvar'
	local original_indepvars `backup_original_indepvars'
	local original_endogvars `backup_original_endogvars'
	local original_instruments `backup_original_instruments'
}
else if ("`stage'"=="ols") {
	local tss = `backup_tss'
	local fast 1
	local depvar `backup_depvar'
	local endogvars
	local indepvars `backup_indepvars' `backup_endogvars'
	local instruments
	local original_depvar `backup_original_depvar'
	local original_endogvars
	local original_indepvars `backup_original_indepvars' `backup_original_endogvars'
	local original_instruments
	local vcesuite avar
}
else if ("`stage'"=="reduced") {
	local tss = `backup_tss'
	local fast 1
	local depvar `backup_depvar'
	local indepvars `backup_indepvars' `backup_instruments'
	local endogvars
	local instruments
	local original_depvar `backup_original_depvar'
	local original_indepvars `backup_original_indepvars' `backup_original_instruments'
	local original_endogvars
	local original_instruments
	local vcesuite avar
}
else if ("`stage'"=="acid") {
	local tss = `backup_tss'
	local fast 1
	local depvar `backup_depvar'
	local indepvars `backup_indepvars' `backup_endogvars' `backup_instruments'
	local endogvars
	local instruments
	local original_depvar `backup_original_depvar'
	local original_indepvars `backup_original_indepvars' `backup_original_endogvars' `backup_original_instruments'
	local original_endogvars
	local original_instruments
	local vcesuite avar
}
else if ("`stage'"=="first") {
	local ++ i_endogvar
	local tss = `tss_`lhs_endogvar''
	local fast 1
	local depvar `lhs_endogvar'
	local indepvars `backup_indepvars' `backup_instruments'
	local endogvars
	local instruments
	local original_depvar : word `i_endogvar' of `backup_original_endogvars'
	local original_indepvars `backup_original_indepvars' `backup_original_endogvars' `backup_original_instruments'
	local original_endogvars
	local original_instruments
	local vcesuite avar
}
**** END OF UGLY -stages- CODE >>>> 

* Cleanup
	ereturn clear

* Regress
	if ("`stage'"=="none") Debug, level(2) msg("(running regresion: `model'.`ivsuite')")
	local avge = cond(`N_avge'>0, "__W*__", "")
	local options
	local option_list ///
		depvar indepvars endogvars instruments avgevars ///
		original_depvar original_indepvars original_endogvars ///
		original_instruments original_absvars avge_targets ///
		vceoption vcetype vcesuite ///
		kk suboptions showraw vceunadjusted first weightexp ///
		estimator twicerobust // Whether to run or not two-step gmm
	foreach opt of local option_list {
		if ("``opt''"!="") local options `options' `opt'(``opt'')
	}

	* Five wrappers in total, two for iv (ivreg2, ivregress), three for ols (regress, avar, mwc)
	local wrapper "Wrapper_`subcmd'" // regress ivreg2 ivregress
	if ("`subcmd'"=="regress" & "`vcesuite'"=="avar") local wrapper "Wrapper_avar"
	if ("`subcmd'"=="regress" & "`vcesuite'"=="mwc") local wrapper "Wrapper_mwc"

	if (!inlist("`stage'","none", "iv")) local wrapper "Wrapper_avar" // Compatible with ivreg2
	Debug, level(3) msg(_n "call to wrapper:" _n as result "`wrapper', `options'")
	`wrapper', `options'
	
	Assert e(tss)<., msg("within tss is missing (wrapper=`wrapper')")
	
	local subpredict = e(predict) // used to recover the FEs

	if ("`weightvar'"!="") {
		qui su `weightvar', mean
		local sumweights = r(sum)
	}

// PART III - RECOVER FEs AND SAVE RESULTS 

if (`fast') {
	* Copy pasted from below
	Debug, level(3) msg("(avoiding -use- of temporary dataset)")
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	FixVarnames `backup_colnames'
	local newnames "`r(newnames)'"
	local prettynames "`r(prettynames)'"
	matrix colnames `b' = `newnames'
}
else {
	assert inlist("`stage'", "iv", "none")
* 1) Restore untransformed dataset
	qui use "`original_vars'", clear

* 2) Recover the FEs

	* Predict will get (e+d) from the equation y=xb+d+e
	tempvar resid_d
	if e(df_m)>0 {
		local score = cond("`model'"=="ols", "score", "resid")
		`subpredict' double `resid_d', `score' // Auto-selects the program based on the estimation method		
	}
	else {
		gen double `resid_d' = `depvar'
	}

	** If the eqn doesn't have a constant, we need to save the mean of the resid in order to add it when predicting xb
	*if (!`addconstant') {
	*	su `resid_d', mean
	*	ereturn `hidden' scalar _cons = r(mean)
	*}

	Debug, level(2) msg("(loaded untransformed variables, predicted residuals)")

	* Absorb the residuals to obtain the FEs (i.e. run a regression on just the resids)
	Debug, level(2) tic(31)
	Demean, varlist(`resid_d') `maximize_options' save_fe(1)
	Debug, level(2) toc(31) msg("mata:make_residual on final model took")
	drop `resid_d'

* 3) Compute corr(FE,xb) (do before rescaling by cvar or deleting)
	if ("`model'"=="ols") {
		tempvar xb
		_predict double `xb', xb // -predict- overwrites sreturn, use _predict if needed
		forv g=1/`N_hdfe' { 
			qui corr `xb' __Z`g'__
			local corr`g' = r(rho)
		}
		drop `xb'
	}

* 4) Replace tempnames in the coefs table
	* (e.g. __00001 -> L.somevar)
	* (this needs to be AFTER predict but before deleting FEs and AvgEs)
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	FixVarnames `backup_colnames'
	local newnames "`r(newnames)'"
	local prettynames "`r(prettynames)'"
	matrix colnames `b' = `newnames'

* 5) Save FEs w/proper name, format
	Save, original_depvar(`original_depvar')
	local keepvars `r(keepvars)'
	if ("`keepvars'"!="") format `fe_format' `keepvars'
	
* 6) Save AvgEs
	forv g=1/`N_avge' {
		local var __W`g'__
		local target : char `var'[target]
		if ("`target'"!="") {
			rename `var' `target'
			local avge_target`g' `target' // Used by -predict-
			local keepvars `keepvars' `target'
		}
	}

	if ("`keepvars'"!="") format `fe_format' `keepvars' // The format of depvar, saved by -Parse-

* 7) Save dataset with FEs and e(sample)
	keep `uid' `keepvars'
	tempfile output
	qui save "`output'"
} // fast

* 8) Restore original dataset and merge
	if (inlist("`stage'","none", "iv")) restore // Restore user-provided dataset (since -iv- comes at the end, that is done at that stage!)
	if (!`fast') {
		// `saved_group' was created by EstimateDoF.ado
		if (!`saved_group')  local groupdta
		SafeMerge, uid(`uid') file("`output'") groupdta("`groupdta'")
		*cap tsset, noquery // we changed -sortby- when we merged (even if we didn't really resort)
	}

// PART IV - ERETURN OUTPUT

	if (`c(version)'>=12) local hidden hidden // ereturn hidden requires v12+

* Ereturns common to all commands
	ereturn local cmd = "reghdfe"
	ereturn local subcmd = cond(inlist("`stage'", "none", "iv"), "`subcmd'", "regress")
	ereturn local cmdline `"`cmdline'"'
	ereturn local model = cond("`gmm2s'"=="", "`model'", "gmm2s")
	ereturn local model = cond("`cue'"=="", "`model'", "cue")
	ereturn local model = cond("`liml'"=="", "`model'", "liml")
	ereturn local dofadjustments = "`dofadjustments'"
	ereturn local title = "HDFE " + e(title)
	ereturn local title2 =  "Absorbing `N_hdfe' HDFE " + plural(`N_hdfe', "indicator")
	ereturn local predict = "reghdfe_p"
	ereturn local estat_cmd = "reghdfe_estat"
	ereturn local footnote = "reghdfe_footnote"
	ereturn local absvars = "`original_absvars'"
	ereturn local vcesuite = "`vcesuite'"
	ereturn local maximize_options = "`maximize_options'" // In option format; tolerance(..) etc.
	if ("`stage'"!="none") ereturn local iv_depvar = "`backup_original_depvar'"
	ereturn `hidden' local diopts = "`diopts'"
	if ("`over'"!="") {
		ereturn local over = "`over'"
		if ("`over_value'"!="") ereturn local over_value = "`over_value'"
		if ("`over_label'"!="") ereturn local over_label = "`over_label'"
		local fixed_absvars = e(absvars)
		local fixed_absvars : subinstr local fixed_absvars "i.`over'#" "", all
		local fixed_absvars : subinstr local fixed_absvars "i.`over'" "", all
		local fixed_absvars `fixed_absvars' // Trim
		ereturn local absvars = "`fixed_absvars'"
	}

	if ("`e(clustvar)'"!="") {
		mata: st_local("clustvar", invtokens(clustervars_original))
		* With kiefer/dkraay we add a time clustervar
		if ("`clustvar'"!="") ereturn local clustvar "`clustvar'"
		ereturn scalar N_clustervars = `num_clusters'
	}

	* Besides each cmd's naming style (e.g. exogr, exexog, etc.) keep one common one
	foreach cat in depvar indepvars endogvars instruments {
		local vars ``cat''
		if ("`vars'"=="") continue
		ereturn local `cat' "`original_`cat''"
	}
	ereturn local avgevars "`avge'" // bugbug?

	ereturn `hidden' local subpredict = "`subpredict'"
	ereturn `hidden' local prettynames "`prettynames'"
	forv g=1/`N_avge' {
		ereturn `hidden' local avge_target`g' "`avge_target`g''" // Used by -predict-
	}

	* Stata uses e(vcetype) for the SE column headers
	* In the default option, leave it empty.
	* In the cluster and robust options, set it as "Robust"
	ereturn local vcetype = proper("`vcetype'") //
	if (e(vcetype)=="Cluster") ereturn local vcetype = "Robust"
	if (e(vcetype)=="Unadjusted") ereturn local vcetype
	if ("`e(vce)'"=="." | "`e(vce)'"=="") ereturn local vce = "`vcetype'" // +-+-
	Assert inlist("`e(vcetype)'", "", "Robust", "Jackknife", "Bootstrap")

	ereturn scalar N_hdfe = `N_hdfe'
	if ("`N_avge'"!="") ereturn scalar N_avge = `N_avge'

* Absorbed-specific returns
	ereturn scalar mobility = `M'
	ereturn scalar df_a = `kk'
	forv g=1/`N_hdfe' {
		ereturn scalar M`g' = `M`g''
		ereturn scalar K`g' = `K`g''
		ereturn `hidden' scalar M`g'_exact = `M`g'_exact' // 1 or 0 whether M`g' was calculated exactly or not
		ereturn `hidden' local corr`g' = "`corr`g''" //  cond("`corr`g''"=="", ., "`corr`g''")
		ereturn `hidden' local hdfe_target`g' = "`hdfe_target`g''"
		ereturn `hidden' local hdfe_cvar`g' = "`hdfe_cvar`g''"
		ereturn `hidden' scalar M`g'_nested = `M`g'_nested'
	}

	Assert e(df_r)<. , msg("e(df_r) is missing")
	ereturn `hidden' scalar tss_within = e(tss)
	ereturn scalar tss = `tss'

	ereturn scalar ll   = -0.5 * (e(N)*ln(2*_pi) + e(N)*ln(e(rss)       /e(N)) + e(N))
	ereturn scalar ll_0 = -0.5 * (e(N)*ln(2*_pi) + e(N)*ln(e(tss_within)/e(N)) + e(N))

	ereturn scalar r2 = 1 - e(rss) / e(tss)
	ereturn scalar r2_within = 1 - e(rss) / e(tss_within)
	ereturn scalar mss = e(tss) - e(rss)

	* ivreg2 uses e(r2c) and e(r2u) for centered/uncetered R2; overwrite first and discard second
	if (e(r2c)!=.) {
		ereturn scalar r2c = e(r2)
		ereturn scalar r2u = .
	}

	* Computing Adj R2 with clustered SEs is tricky because it doesn't use the adjusted inputs:
	* 1) It uses N instead of N_clust
	* 2) For the DoFs, it uses N - Parameters instead of N_clust-1
	* 3) Further, to compute the parameters, it includes those nested within clusters
	
	* Note that this adjustment is NOT PERFECT because we won't compute the mobility groups just for improving the r2a
	* (when a FE is nested within a cluster, we don't need to compute mobilty groups; but to get the same R2a as other estimators we may want to do it)
	* Instead, you can set by hand the dof() argument and remove -cluster- from the list

	if ("`model'"=="ols" & `num_clusters'>0) Assert e(unclustered_df_r)<., msg("wtf-`vcesuite'")
	local used_df_r = cond(e(unclustered_df_r)<., e(unclustered_df_r), e(df_r)) - `M_due_to_nested'
	ereturn scalar r2_a = 1 - (e(rss)/`used_df_r') / ( e(tss) / (e(N)-1) )
	ereturn scalar rmse = sqrt( e(rss) / `used_df_r' )

	ereturn scalar r2_a_within = 1 - (e(rss)/`used_df_r') / ( e(tss_within) / (`used_df_r'+e(df_m)) )

	if (e(N_clust)<.) Assert e(df_r) == e(N_clust) - 1, msg("Error, `wrapper' should have made sure that N_clust-1==df_r")
	*if (e(N_clust)<.) ereturn scalar df_r = e(N_clust) - 1

	if ("`weightvar'"!="") ereturn scalar sumweights = `sumweights'

	if ("`model'"=="ols" & inlist("`vcetype'", "unadjusted", "ols")) {
		ereturn scalar F_absorb = (e(r2)-`r2c') / (1-e(r2)) * e(df_r) / (`kk'-1) // -1 b/c we exclude constant for this
		if (`nested') {
			local rss`N_hdfe' = e(rss)
			local temp_dof = e(N) - e(df_m) // What if there are absorbed collinear with the other RHS vars?
			local j 0
			ereturn `hidden' scalar rss0 = `rss0'
			forv g=1/`N_hdfe' {
				local temp_dof = `temp_dof' - e(K`g') + e(M`g')
				*di in red "g=`g' RSS=`rss`g'' and was `rss`j''.  dof=`temp_dof'"
				ereturn `hidden' scalar rss`g' = `rss`g''
				ereturn `hidden' scalar df_a`g' = e(K`g') - e(M`g')
				local df_a_g = e(df_a`g') - (`g'==1)
				ereturn scalar F_absorb`g' = (`rss`j''-`rss`g'') / `rss`g'' * `temp_dof' / `df_a_g'
				ereturn `hidden' scalar df_r`g' = `temp_dof'
				local j `g'
			}   
		}
	}

	// There is a big assumption here, that the number of other parameters does not increase asymptotically
	// BUGBUG: We could allow the option to indicate what parameters do increase asympt.

	if ("`savefirst'"!="") ereturn `hidden' scalar savefirst = `savefirst'

	* We have to replace -unadjusted- or else subsequent calls to -suest- will fail
	Subtitle `vceoption' // will set title2, etc. Run after e(bw) and all the others are set!
	if (e(vce)=="unadjusted") ereturn local vce = "ols"

	if ("`stages'"!="none") {
		ereturn local stage = "`stage'"
		ereturn `hidden' local stages = "`stages'"
	}

* Show table and clean up
	ereturn repost b=`b', rename // why here???

	if ("`stage'"!="none") Debug, level(0) msg(_n "{title:Stage: `stage'}" _n)
	if ("`lhs_endogvar'"!="<none>") Debug, level(0) msg("{title:Endogvar: `lhs_endogvar'}")
	Replay
	Attach, notes(`notes') statsmatrix(`statsmatrix') summarize_quietly(`summarize_quietly')

*** <<<< LAST PART OF UGLY STAGE <<<<	
if (!inlist("`stage'","none", "iv")) {
	local estimate_name reghdfe_`stage'`i_endogvar'
	local stored_estimates `stored_estimates' `estimate_name'
	local cmd estimates store `estimate_name', nocopy
	Debug, level(2) msg(" - Storing estimate: `cmd'")
	`cmd'
}
else if ("`stage'"=="iv") {
	* On the last stage, save list of all stored estimates
	assert "`stored_estimates'"!=""
	ereturn `hidden' local stored_estimates = "`stored_estimates'"
}

} // lhs_endogvar
} // stage
*** >>>> LAST PART OF UGLY STAGE >>>>

	Stop

end

// -------------------------------------------------------------------------------------------------

* The idea of this program is to keep the sort order when doing the merges

program define SafeMerge, eclass sortpreserve
syntax, uid(varname numeric) file(string) [groupdta(string)]
	* Merging gives us e(sample) and the FEs / AvgEs
	tempvar merge
	merge 1:1 `uid' using "`file'", assert(master match) nolabel nonotes noreport gen(`merge')
	
	* Add e(sample) from _merge
	tempvar sample
	gen byte `sample' = (`merge'==3)
	la var `sample' "[HDFE Sample]"
	ereturn repost , esample(`sample')
	drop `merge'

	* Add mobility group
	if ("`groupdta'"!="") merge 1:1 `uid' using "`groupdta'", assert(master match) nogen nolabel nonotes noreport sorted
end

program define Subtitle, eclass
	* Fill e(title3/4/5) based on the info of the other e(..)

	if (inlist("`e(vcetype)'", "Robust", "Cluster")) local hacsubtitle1 "heteroskedasticity"
	if ("`e(kernel)'"!="" & "`e(clustvar)'"=="") local hacsubtitle3 "autocorrelation"
	if ("`e(kiefer)'"!="") local hacsubtitle3 "within-cluster autocorrelation (Kiefer)"
	if ("`hacsubtitle1'"!="" & "`hacsubtitle3'" != "") local hacsubtitle2 " and "
	local hacsubtitle "`hacsubtitle1'`hacsubtitle2'`hacsubtitle3'"
	if strlen("`hacsubtitle'")>30 {
		local hacsubtitle : subinstr local hacsubtitle "heteroskedasticity" "heterosk.", all word
		local hacsubtitle : subinstr local hacsubtitle "autocorrelation" "autocorr.", all word
	}
	if ("`hacsubtitle'"!="") {
		ereturn local title3 = "Statistics robust to `hacsubtitle'"
		
		if ("`e(kernel)'"!="") local notes " `notes' kernel=`e(kernel)'"
		if ("`e(bw)'"!="") local notes " `notes' bw=`e(bw)'"
		if ("`e(dkraay)'"!="") local notes " `notes' dkraay=`e(dkraay)'"
		local notes `notes' // remove initial space
		if ("`notes'"!="") ereturn local title4 = " (`notes')"
		if ("`notes'"!="") {
			if ("`_dta[_TSpanel]'"!="") local tsset panel=`_dta[_TSpanel]'
			if ("`_dta[_TStvar]'"!="") local tsset `tsset' time=`_dta[_TStvar]'
			local tsset `tsset'
			ereturn local title5 = " (`tsset')"
		}
	}
end

	
// -------------------------------------------------------------
// Parsing and basic sanity checks for REGHDFE.ado
// -------------------------------------------------------------

program define Parse

* Remove extra spacing from cmdline (just for aesthetics)
	mata: st_local("cmdline", stritrim(`"reghdfe `0'"') )
	ereturn clear // Clear previous results and drops e(sample)

* Parse the broad syntax (also see map_init(), ParseAbsvars.ado, ParseVCE.ado, etc.)
	syntax anything(id="varlist" name=0 equalok) [if] [in] [fw aw pw/] , ///
	/// Main Options ///
		Absorb(string) ///
		[ ///
		VCE(string) ///
	/// Seldom Used ///
		DOFadjustments(string) ///
		GROUPVAR(name) /// Variable that will contain the first connected group between FEs
	/// Optimization /// Defaults are handled within Mata		
		GROUPsize(string) /// Process variables in batches of #
		TRANSFORM(string) ///
		ACCELeration(string) ///
		Verbose(string) ///
		TOLerance(string) ///
		MAXITerations(string) ///
		KEEPSINGLETONS(string) /// (UNDOCUMENTED) Will keep singletons
		CHECK /// TODO: Implement
		FAST /// TODO: Implement
	/// Regression ///
		ESTimator(string) /// GMM2s CUE LIML
		IVsuite(string) ///
		SAVEFIRST ///
		FIRST ///
		SHOWRAW ///
		VCEUNADJUSTED /// (UNDOCUMENTED) Option when running gmm2s with ivregress; will match results of ivreg2
		SMALL Hascons TSSCONS /// ignored options
		SUBOPTions(string) /// Options to be passed to the estimation command (e.g . to regress)
	/// Multiple regressions in one go ///
		OVER(varname numeric) CLEAR ///
		NESTED /// TODO: Implement
		STAGEs(string) ///
	/// Miscellanea ///
		NOTES(string) /// NOTES(key=value ..)
		] [*] // For display options ; and SUmmarize(stats)

	local allkeys cmdline if in

* Parse varlist: depvar indepvars (endogvars = iv_vars)
	ParseIV `0', estimator(`estimator') ivsuite(`ivsuite') `savefirst' `first' `showraw' `vceunadjusted' `small'
	local keys subcmd model ivsuite estimator depvar indepvars endogvars instruments fe_format ///
		savefirst first showraw vceunadjusted basevars
	foreach key of local keys {
		local `key' "`s(`key')'"
	}
	local allkeys `allkeys' `keys'

* Weights
	if ("`weight'"!="") {
		local weightvar `exp'
		local weighttype `weight'
		local weightexp [`weight'=`weightvar']
		confirm var `weightvar', exact // just allow simple weights

		* Check that weights are correct (e.g. with fweight they need to be integers)
		local num_type = cond("`weight'"=="fweight", "integers", "reals")
		local basenote "weight {res}`weightvar'{txt} can only contain strictly positive `num_type', but"
		qui cou if `weightvar'<0
		Assert (`r(N)'==0), msg("`basenote' `r(N)' negative values were found!")
		qui cou if `weightvar'==0
		if (`r(N)'>0) di as text "`basenote' `r(N)' zero values were found (will be dropped)"
		qui cou if `weightvar'>=.
		if (`r(N)'>0) di as text "`basenote' `r(N)' missing values were found (will be dropped)"
		if ("`weight'"=="fweight") {
			qui cou if mod(`weightvar',1) & `weightvar'<.
			Assert (`r(N)'==0), msg("`basenote' `r(N)' non-integer values were found!")
		}
	}
	local allkeys `allkeys' weightvar weighttype weightexp

* Parse VCE options: (Needs to be BEFORE ParseAbsvars, because of -clustervars-)
	mata: st_local("hascomma", strofreal(strpos("`vce'", ","))) // is there a commma already in `vce'?
	if (!`hascomma') local vce `vce' ,
	ParseVCE `vce' weighttype(`weighttype') // Might call map_init_*()
	local keys vceoption vcetype vcesuite vceextra num_clusters clustervars bw kernel dkraay twicerobust
	foreach key of local keys {
		local `key' "`s(`key')'"
	}
	local allkeys `allkeys' `keys'

* Parse Absvars and optimization options
	ParseAbsvars `absorb' // Stores results in r()
	local absorb_keepvars `r(all_ivars)' `r(all_cvars)'
	local N_hdfe `r(G)'
	mata: HDFE_S = map_init() // Reads results from r()
	local allkeys `allkeys' absorb_keepvars N_hdfe

	* Tell Mata what weightvar we have
	if ("`weightvar'"!="") mata: map_init_weights(HDFE_S, "`weightvar'", "`weighttype'")

	* Time/panel variables (need to give them to Mata)
	local panelvar `_dta[_TSpanel]'
	local timevar `_dta[_TStvar]'

	* Parse optimization options (pass them to map_init_*)
	* String options
	local optlist transform acceleration clustervars panelvar timevar
	foreach opt of local optlist {
		if ("``opt''"!="") mata: map_init_`opt'(HDFE_S, "``opt''")
	}
	local allkeys `allkeys' `optlist'

	* Numeric options
	local optlist groupsize verbose tolerance maxiterations keepsingletons
	foreach opt of local optlist {
		if ("``opt''"!="") mata: map_init_`opt'(HDFE_S, ``opt'')
	}
	local allkeys `allkeys' `optlist'

	local fast = cond("`fast'"!="", 1, 0) // 1=Yes
	local allkeys `allkeys' fast

* DoF Adjustments
	if ("`dofadjustments'"=="") local dofadjustments all
	ParseDOF , `dofadjustments'
	local dofadjustments "`s(dofadjustments)'"
	* Mobility groups
	if ("`groupvar'"!="") conf new var `groupvar'
	local allkeys `allkeys' dofadjustments groupvar

* Parse summarize option: [summarize | summarize( stats... [,QUIetly])]
	* Note: ParseImplicit deals with "implicit" options and fills their default values
	local default_stats mean min max
	ParseImplicit, opt(SUmmarize) default(`default_stats') input(`options') syntax([namelist(name=stats)] , [QUIetly]) inject(stats quietly)
	local summarize_quietly = ("`quietly'"!="")
	if ("`stats'"=="" & "`quietly'"!="") local stats `default_stats'
	local allkeys `allkeys' stats summarize_quietly

* Parse over() option
	local clear = "`clear'"!=""
 	if ("`over'"!="") {
		unab over : `over', max(1)
		Assert (`clear'), msg("over() requires the -clear- option")
	}
	local allkeys `allkeys' over clear

* Nested
	local nested = cond("`nested'"!="", 1, 0) // 1=Yes
	if (`nested' & !("`model'"=="ols" & "`vcetype'"=="unadjusted") ) {
		Debug, level(0) msg("(option nested ignored, only works with OLS and conventional/unadjusted VCE)") color("error")
	}
	local allkeys `allkeys' nested

* Stages
	assert "`model'"!="" // just to be sure this goes after `model' is set
	local iv_stage iv
	local stages : list stages - iv_stage
	local valid_stages ols first acid reduced
	local wrong_stages : list stages - valid_stages
	Assert "`wrong_stages'"=="", msg("Error, invalid stages(): `wrong_stages'")
	if ("`stages'"!="") {
		Assert "`model'"=="iv", msg("Error, stages() only valid with an IV regression")
		local stages `stages' `iv_stage' // Put -iv- *last* (so it does the -restore-; note that we don't need it first to trim MVs b/c that's done earlier)
	}
	else {
		local stages none // So we can loop over stages
	}
	local allkeys `allkeys' stages

* Parse Coef Table Options (do this last!)
	_get_diopts diopts options, `options' // store in `diopts', and the rest back to `options'
	Assert `"`options'"'=="", msg(`"invalid options: `options'"')
	if ("`hascons'`tsscons'"!="") di in ye "(option `hascons'`tsscons' ignored)"
	local allkeys `allkeys' diopts


* Other keys:
	local allkeys `allkeys' suboptions notes
	// Missing keys: check

* Return values
	Debug, level(3) newline
	Debug, level(3) msg("{title:Parsed options:}")
	foreach key of local allkeys {
		if (`"``key''"'!="") Debug, level(3) msg("  `key' = " as result `"``key''"')
		c_local `key' `"``key''"' // Inject values into caller (reghdfe.ado)
	}
end

program define ParseIV, sclass
	syntax anything(id="varlist" name=0 equalok), [ ///
		estimator(string) ivsuite(string) ///
		savefirst first showraw vceunadjusted small ]

	* Parses varlist: depvar indepvars [(endogvars = instruments)]
		* depvar: dependent variable
		* indepvars: included exogenous regressors
		* endogvars: included endogenous regressors
		* instruments: excluded exogenous regressors

	* Model: OLS or IV-type?
	local model ols
	foreach _ of local 0 {
		if (substr(`"`_'"', 1, 1)=="(") {
			local model iv
			continue, break
		}
	}

	* IV Suite
	if ("`model'"=="iv") {
		if ("`ivsuite'"=="") local ivsuite ivreg2 // Set default
		Assert inlist("`ivsuite'","ivreg2","ivregress") , msg("error: wrong IV routine (`ivsuite'), valid options are -	ivreg2- and -ivregress-")
		cap findfile `ivsuite'.ado
		Assert !_rc , msg("error: -`ivsuite'- not installed, please run {stata ssc install `ivsuite'} or change the option 	-ivsuite-")

		local savefirst = ("`savefirst'"!="")
		local first = ("`first'"!="")
		if (`savefirst') Assert `first', msg("Option -savefirst- requires -first-")
		local subcmd `ivsuite'
	}
	else {
		local savefirst
		local first
		local subcmd regress
	}

	* Estimator
	if ("`estimator'"=="" & "`model'"=="iv") local estimator 2sls // Set default
	if ("`estimator'"!="") {
		Assert "`model'"=="iv", ///
			msg("reghdfe error: estimator() requires an instrumental-variable regression")
		if (substr("`estimator'", 1, 3)=="gmm") local estimator gmm2s
		Assert inlist("`estimator'", "2sls", "gmm2s", "liml", "cue"), ///
			msg("reghdfe error: invalid estimator `estimator'")
		if ("`estimator'"=="cue") Assert "`ivsuite'"=="ivreg2", ///
			msg("reghdfe error: estimator `estimator' only available with the ivreg2 command, not ivregress")
	}

	* For this, _iv_parse would have been useful, but I don't want to do factor expansions when parsing
	if ("`model'"=="iv") {

		* get part before parentheses
		local wrongparens 1
		while (`wrongparens') {
			gettoken tmp 0 : 0 ,p("(")
			local left `left'`tmp'
			* Avoid matching the parens of e.g. L(-1/2) and L.(var1 var2)
			* Using Mata to avoid regexm() and trim() space limitations
			mata: st_local("tmp1", subinstr("`0'", " ", "") ) // wrong parens if ( and then a number
			mata: st_local("tmp2", substr(strtrim("`left'"), -1) ) // wrong parens if dot
			local wrongparens = regexm("`tmp1'", "^\([0-9-]") | ("`tmp2'"==".")
			if (`wrongparens') {
				gettoken tmp 0 : 0 ,p(")")
				local left `left'`tmp'
			}
		}

		* get part in parentheses
		gettoken right 0 : 0 ,bind match(parens)
		Assert trim(`"`0'"')=="" , msg("error: remaining argument: `0'")

		* now parse part in parentheses
		gettoken endogvars instruments : right ,p("=")
		gettoken equalsign instruments : instruments ,p("=")

		Assert "`endogvars'"!="", msg("iv: endogvars required")
		local 0 `endogvars'
		syntax varlist(fv ts numeric)
		local endogvars `varlist'

		Assert "`instruments'"!="", msg("iv: instruments required")
		local 0 `instruments'
		syntax varlist(fv ts numeric)
		local instruments `varlist'
		
		local 0 `left' // So OLS part can handle it
	}

* OLS varlist
	syntax varlist(fv ts numeric)
	gettoken depvar indepvars : varlist
	_fv_check_depvar `depvar'

* Extract format of depvar so we can format FEs like this
	fvrevar `depvar', list
	local fe_format : format `r(varlist)' // The format of the FEs and AvgEs that will be saved

* Variables shouldn't be repeated
* This is not perfect (e.g. doesn't deal with "x1-x10") but still helpful
	local allvars `depvar' `indepvars' `endogvars' `instruments'
	local dupvars : list dups allvars
	Assert "`dupvars'"=="", msg("error: there are repeated variables: <`dupvars'>")

* More IV options
	if ("`small'"!="") di in ye "(note: reghdfe will always use the option -small-, no need to specify it)"

* Parse -showraw- : shows raw output of called subcommand (e.g. ivreg2)
	local showraw = ("`showraw'"!="")

* Parse -unadjusted-
	* If true, will use wmatrix(...) vce(unadjusted) instead of the default of setting vce contents equal to wmatrix
	* This basically undoes the extra adjustment that ivregress does, so it's comparable with ivreg2
	*
	* Note: Cannot match exactly the -ivregress- results without vceunadjusted (see test-gmm.do)
	* Thus, I will set this to true ALWAYS
	local vceunadjusted = 1 // ("`vceunadjusted'"!="")

* Get base variables of time and factor variables (e.g. i.foo L(1/3).bar -> foo bar)
	foreach vars in depvar indepvars endogvars instruments {
		if ("``vars''"!="") {
			fvrevar ``vars'' , list
			local basevars `basevars' `r(varlist)'
		}
	}

	local keys subcmd model ivsuite estimator depvar indepvars endogvars instruments fe_format ///
		savefirst first showraw vceunadjusted basevars
	foreach key of local keys {
		sreturn local `key' ``key''
	}
end 

program define ParseVCE, sclass
	* Note: bw=1 *usually* means just do HC instead of HAC
	* BUGBUG: It is not correct to ignore the case with "bw(1) kernel(Truncated)"
	* but it's too messy to add -if-s everywhere just for this rare case (see also Mark Schaffer's email)

	syntax 	[anything(id="VCE type")] , ///
			[bw(integer 1) KERnel(string) dkraay(integer 1) kiefer] ///
			[suite(string) TWICErobust] ///
			[weighttype(string)]

	if ("`anything'"=="") local anything unadjusted
	Assert `bw'>0, msg("VCE bandwidth must be a positive integer")
	gettoken vcetype clustervars : anything
	* Expand variable abbreviations; but this adds unwanted i. prefixes
	if ("`clustervars'"!="") {
		fvunab clustervars : `clustervars'
		local clustervars : subinstr local clustervars "i." "", all
	}

	* vcetype abbreviations:
	if (substr("`vcetype'",1,3)=="ols") local vcetype unadjusted
	if (substr("`vcetype'",1,2)=="un") local vcetype unadjusted
	if (substr("`vcetype'",1,1)=="r") local vcetype robust
	if (substr("`vcetype'",1,2)=="cl") local vcetype cluster
	if ("`vcetype'"=="conventional") local vcetype unadjusted // Conventional is the name given in e.g. xtreg
	Assert strpos("`vcetype'",",")==0, msg("Unexpected contents of VCE: <`vcetype'> has a comma")

	* Sanity checks on vcetype
	if ("`vcetype'"=="" & "`weighttype'"=="pweight") local vcetype robust
	Assert !("`vcetype'"=="unadjusted" & "`weighttype'"=="pweight"), msg("pweights do not work with unadjusted errors, use a different vce()")
	if ("`vcetype'"=="") local vcetype unadjusted
	Assert inlist("`vcetype'", "unadjusted", "robust", "cluster"), msg("VCE type not supported: `vcetype'")

	* Cluster vars
	local num_clusters : word count `clustervars'
	Assert inlist( (`num_clusters'>0) + ("`vcetype'"=="cluster") , 0 , 2), msg("Can't specify cluster without clustervars and viceversa") // XOR

	* VCE Suite
	local vcesuite `suite'
	if ("`vcesuite'"=="") local vcesuite default
	if ("`vcesuite'"=="default") {
		if (`bw'>1 | `dkraay'>1 | "`kiefer'"!="" | "`kernel'"!="") {
			local vcesuite avar
		}
		else if (`num_clusters'>1) {
			local vcesuite mwc
		}
	}

	Assert inlist("`vcesuite'", "default", "mwc", "avar"), msg("Wrong vce suite: `vcesuite'")

	if ("`vcesuite'"=="mwc") {
		cap findfile tuples.ado
		Assert !_rc , msg("error: -tuples- not installed, please run {stata ssc install tuples} to estimate multi-way clusters.")
	}
	
	if ("`vcesuite'"=="avar" | "`stages'"!="none") {
		cap findfile avar.ado // We use -avar- as default with stages (on the non-iv stages)
		Assert !_rc , msg("error: -avar- not installed, please run {stata ssc install avar} or change the option -vcesuite-")
	}

	* Some combinations are not coded
	Assert !("`ivsuite'"=="ivregress" & (`num_clusters'>1 | `bw'>1 | `dkraay'>1 | "`kiefer'"!="" | "`kernel'"!="") ), msg("option vce(`vce') incompatible with ivregress")
	Assert !("`ivsuite'"=="ivreg2" & (`num_clusters'>2) ), msg("ivreg2 doesn't allow more than two cluster variables")
	Assert !("`model'"=="ols" & "`vcesuite'"=="avar" & (`num_clusters'>2) ), msg("avar doesn't allow more than two cluster variables")
	Assert !("`model'"=="ols" & "`vcesuite'"=="default" & (`bw'>1 | `dkraay'>1 | "`kiefer'"!="" | "`kernel'"!="") ), msg("to use those vce options you need to use -avar- as the vce suite")
	if (`num_clusters'>0) local temp_clustervars " <CLUSTERVARS>"
	if (`bw'==1 & `dkraay'==1 & "`kernel'"!="") local kernel // No point in setting kernel here 
	if (`bw'>1 | "`kernel'"!="") local vceextra `vceextra' bw(`bw') 
	if (`dkraay'>1) local vceextra `vceextra' dkraay(`dkraay') 
	if ("`kiefer'"!="") local vceextra `vceextra' kiefer 
	if ("`kernel'"!="") local vceextra `vceextra' kernel(`kernel')
	if ("`vceextra'"!="") local vceextra , `vceextra'
	local vceoption "`vcetype'`temp_clustervars'`vceextra'" // this excludes "vce(", only has the contents

	if ("`vceextra'"!="") mata: map_init_vce_is_hac(HDFE_S, 1)

	local keys vceoption vcetype vcesuite vceextra num_clusters clustervars bw kernel dkraay twicerobust
	foreach key of local keys {
		sreturn local `key' ``key''
	}
end

program define ParseAbsvars, rclass
syntax anything(id="absvars" name=absvars equalok everything), [SAVEfe] // [CLUSTERvars(varlist numeric fv)]
	* Logic: split absvars -> expand each into factors -> split each into parts

	local g 0
	local all_cvars
	local all_ivars

	while ("`absvars'"!="") {
		local ++g
		gettoken absvar absvars : absvars, bind
		local target
		if strpos("`absvar'","=") gettoken target absvar : absvar, parse("=")
		if ("`target'"!="") {
			conf new var `target'
			gettoken eqsign absvar : absvar, parse("=")
		}

		local n : word count absvar
		local hasdot = strpos("`absvar'", ".")
		local haspound = strpos("`absvar'", "#")
		if (`n'==1 & !`hasdot' & !`haspound') local absvar i.`absvar'
		
		local 0 `absvar'
		syntax varlist(numeric fv)
			//di as error "    varlist=<`varlist'>"
		
		local ivars
		local cvars
		local has_intercept 0
		foreach factor of local varlist {
			local hasdot = strpos("`factor'", ".")
			local haspound = strpos("`factor'", "#")
			local factor_has_cvars 0

			if (!`hasdot') continue
			while ("`factor'"!="") {
				gettoken part factor : factor, parse("#")
				local is_indicator = strpos("`part'", "i.")
				local is_continuous = strpos("`part'", "c.")
				local basevar = substr("`part'", 3, .)
				if (`is_indicator') local ivars `ivars' `basevar'
				if (`is_continuous') {
					local cvars `cvars' `basevar'
					local factor_has_cvars 1
				}
			}
			if (!`factor_has_cvars') local has_intercept 1

		}
		
		local ivars : list uniq ivars
		local num_slopes : word count `cvars'
		Assert "`ivars'"!="", msg("error parsing absvars: no indicator variables in absvar <`absvar'>")
		local unique_cvars : list uniq cvars
		Assert (`: list unique_cvars == cvars'), msg("error parsing absvars: factor interactions such as i.x##i.y not allowed")

		local all_cvars `all_cvars' `cvars'
		local all_ivars `all_ivars' `ivars'

		return local target`g' `target'
		return local ivars`g' `ivars'
		return local cvars`g' `cvars'
		return scalar has_intercept`g' = `has_intercept'
		return scalar num_slopes`g' = `num_slopes'
	
		local label : subinstr local ivars " " "#", all
		if (`num_slopes'==1) {
			local label `label'#c.`cvars'
		}
		else if (`num_slopes'>1) {
			local label `label'#c.(`cvars')
		}
		return local varlabel`g' `label'
	
	}
	
	local all_ivars : list uniq all_ivars
	local all_cvars : list uniq all_cvars

	return scalar G = `g'
	return scalar savefe = ("`savefe'"!="")
	return local all_ivars `all_ivars'
	return local all_cvars `all_cvars'
end

program define ParseDOF, sclass
	syntax, [ALL NONE] [PAIRwise FIRSTpair] [CLusters] [CONTinuous]
	opts_exclusive "`all' `none'" dofadjustments
	opts_exclusive "`pairwise' `firstpair'" dofadjustments
	if ("`none'"!="") {
		Assert "`pairwise'`firstpair'`clusters'`continuous'"=="", msg("option {bf:dofadjustments()} invalid; {bf:none} not allowed with other alternatives")
		local dofadjustments
	}
	if ("`all'"!="") {
		Assert "`pairwise'`firstpair'`clusters'`continuous'"=="", msg("option {bf:dofadjustments()} invalid; {bf:all} not allowed with other alternatives")
		local dofadjustments pairwise clusters continuous
	}
	else {
		local dofadjustments `pairwise' `firstpair' `clusters' `continuous'
	}
	sreturn local dofadjustments "`dofadjustments'"
end

program define ParseImplicit
* Parse options in the form NAME|NAME(arguments)
	* opt()			name of the option (so if opt=spam, we can have spam or spam(...))
	* default()		default value for the implicit form (in case we don't have a parenthesis)
	* syntax()		syntax of the contents of the parenthesis
	* input()		text to parse (usually `options', the result of a previous syntax .. , .. [*] )
	* inject()		what locals to inject on the caller (depend on -syntax)
	* xor			opt is mandatory (one of the two versions must occur)
	syntax, opt(name local) default(string) syntax(string asis) [input(string asis)] inject(namelist local) [xor]

	* First see if the implicit version is possible
	local lower_opt = lower("`opt'")
	local 0 , `input'
	cap syntax, `opt' [*]
	if ("`xor'"=="") local capture capture
	local rc = _rc
	if (`rc') {
		`capture' syntax, `opt'(string asis) [*]
		if ("`capture'"!="" & _rc) exit
	}
	else {
		local `lower_opt' `default'
	}
	local 0 ``lower_opt''
	syntax `syntax'
	foreach loc of local inject {
		c_local `loc' ``loc''
	}
	c_local options `options'
end

	
// -------------------------------------------------------------------------------------------------
// Expand factor time-series variables
// -------------------------------------------------------------------------------------------------
* Steps:
* 1) Call -fvrevar-
* 2) Label newly generated temporary variables
* 3) Drop i) omitted variables, and ii) base variables (if not part of a #c.var interaction)

program define ExpandFactorVariables, rclass
syntax varlist(min=1 numeric fv ts) [if] [,setname(string)] [CACHE] verbose(integer)
	
	* If saving the data for later regressions -savecache(..)- we will need to match each expansion to its newvars
	* This mata array is used for that
	* Note: This explains why we need to wrap -fvrevar- in a loop
	if ("`cache'"!="") mata: varlist_cache = asarray_create()

	local expanded_msg `"" - variable expansion for `setname': {res}`varlist'{txt} ->""'
	while (1) {
		gettoken factorvar varlist : varlist, bind
		if ("`factorvar'"=="") continue, break

		fvrevar `factorvar' `if' // , stub(__V__) // stub doesn't work in Stata 11.2
		local contents
		foreach var of varlist `r(varlist)' {
			LabelRenameVariable `var' // Tempvars not renamed will be dropped automatically
			if !r(is_dropped) local contents `contents' `r(varname)'
			* Yellow=Already existed, White=Created, Red=NotCreated (omitted or base)
			local color = cond(r(is_dropped), "error", cond(r(is_newvar), "input", "result"))
			if (`verbose'>3) {
				local expanded_msg `"`expanded_msg' as `color' " `r(name)'" as text " (`r(varname)')""'
			}
		}
		Assert "`contents'"!="", msg("error: variable -`fvvar'- in varlist -`varlist'- in category -`setname'- is  empty after factor/time expansion")
		if ("`cache'"!="") mata: asarray(varlist_cache, "`fvvar'", "`contents'")
		local newvarlist `newvarlist' `contents'
	}

	Debug, level(4) msg(`expanded_msg')
	return clear
	return local varlist "`newvarlist'"
end

program define LabelRenameVariable, rclass
syntax varname
	local var `varlist'
	local fvchar : char `var'[fvrevar]
	local tschar : char `var'[tsrevar]
	local is_newvar = ("`fvchar'`tschar'"!="") & substr("`var'", 1, 2)=="__"
	local name "`var'"
	local will_drop 0

	if (`is_newvar') {
		local name "`fvchar'`tschar'"
		local parts : subinstr local fvchar "#" " ", all
		local has_cont_interaction = strpos("`fvchar'", "c.")>0
		local is_omitted 0
		local is_base 0
		foreach part of local parts {
			if (regexm("`part'", "b.*\.")) local is_base 1
			if (regexm("`part'", "o.*\.")) local is_omitted 1
		}

		local will_drop = (`is_omitted') | (`is_base' & !`has_cont_interaction')
		if (!`will_drop') {
			char `var'[name] `name'
			la var `var' "[TEMPVAR] `name'"
			local newvar : subinstr local name "." "__", all
			local newvar : subinstr local newvar "#" "_X_", all
			* -permname- selects newname# if newname is taken (# is the first number available)
			local newvar : permname __`newvar', length(30)
			rename `var' `newvar'
			local var `newvar'
		}
	}

	return scalar is_newvar = `is_newvar'
	return scalar is_dropped = `will_drop'
	return local varname "`var'"
	return local name "`name'"
end

// -------------------------------------------------------------------------------------------------

