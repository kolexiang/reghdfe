capture program drop Prepare
program define Prepare, sclass

syntax, depvar(string) stages(string) model(string) expandedvars(string) vcetype(string) [weightexp(string) endogvars(string)]

* Save the statistics we need before transforming the variables
	* Compute TSS of untransformed depvar
	local tmpweightexp = subinstr("`weightexp'", "[pweight=", "[aweight=", 1)
	qui su `depvar' `tmpweightexp' // BUGBUG: Is this correct?!
	c_local tss = r(Var)*(r(N)-1)
	if (`: list posof "first" in stages') {
		foreach var of varlist `endogvars' {
			qui su `var' `tmpweightexp' // BUGBUG: Is this correct?!
			c_local tss_`var' = r(Var)*(r(N)-1)
		}
	}

* (optional) Compute R2/RSS to run nested Ftests on the FEs
	* a) Compute R2 of regression without FE, to build the joint FTest for all the FEs
	* b) Also, compute RSS of regressions with less FEs so we can run nested FTests on the FEs
	if ("`model'"=="ols" & inlist("`vcetype'", "unadjusted", "ols")) {
		qui _regress `expandedvars' `weightexp', noheader notable
		c_local r2c = e(r2)
	}
end
