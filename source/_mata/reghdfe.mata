clear all

* Type Aliases
local Varlist 		string scalar
local Varname 		string scalar

local Integer 		real scalar
local VarByFE 		real colvector // Should be levels*1
local Series		real colvector // Should be N*1
local Matrix		real matrix

local SharedData 	external struct FixedEffect vector

mata:
mata set matastrict on

void assert_msg(real scalar t, | string scalar msg)
{
	if (args()<2 | msg=="") msg = "assertion is false"
        if (t==0) _error(msg)
}

struct FixedEffect {
	`Integer'	order 			// "g", the position in varlist
	`Integer'	num_slopes
	`Integer'	has_intercept
	`Integer'	levels			// Number of categories spanned by the ivars
	`Varlist'	ivars			// number of i.var elements
	`Varlist'	cvars			// number of c.var elements or slopes
	`Varname'	varlabel		// Original label of this absvar
	`Varname'	estimates		// Name of the variable that will hold the estimates for the FE

}

struct MapProblem {
	struct FixedEffect	vector fixed_effects
	`Integer'			G 			// Number of FEs when bunching slopes
	`Integer'			G_expanded 	// Number of FEs incl. slopes
	`Varname'			weightvar 	// Variable contaning the fw/pw/aw
	real scalar foobar
	real scalar spam

	* fixed_effects = J(G, 1, FixedEffect())
}

struct MapProblem scalar function mapsolve_init(string scalar absvars) {
	struct MapProblem scalar 	S
	string vector				toks
	`Integer'					g, G, Kg
	`Integer'					G_expanded // Counts fixed slopes as independent categories
	
	toks = tokens(absvars)
	G = cols(toks)

	tok

	S.foobar = 2
	S.spam = 7
	return(S)
}

function borrar() {
	struct MapProblem scalar x // or transmorphic
	x = mapsolve_init("turn i.trunk")
	x.foobar + x.spam
}



end

clear
sysuse auto
mata: borrar()


* S = mapsolve_init(transform={Kaczmarz,Cimmino,Symmetric_kaczmarz} , acceleration=...lol, tol=..., fe_varlist=...?)
* project(series, g, S)
* transform(series, S)
* .. = mapsolve(series, S)
* 



// -------------------------------------------------------------------------------------------------
// References
// -------------------------------------------------------------------------------------------------
/*

Hernández-Ramos, Luis M., René Escalante, and Marcos Raydan. "Unconstrained optimization techniques for the acceleration of alternating projection methods." Numerical functional analysis and optimization 32.10 (2011): 1041-1066.

*/



/*

como me bajo singletons recursivamente?
como genero los IDs?

presort por (fe1 fe2 ..)?
asi me aseguro q venga presorteado en una dim aunque sea

grabar los fes como 
i) pares de variables (FE1,Z1). Cuando FE# es 0, es dropped por singleton
ii) alphas en e(a#) , donde la primera columna es alpha, y las columnas 2...H son para cada ivar (no para cvars claro!) 
WAIT PUTA MARE, el max matsize es 11k
Asi que e(a) entonces te da nomas el nombre de la matriz de mata L*(1+#ivar)
*/
