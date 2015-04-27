clear all
cls
set more off

*use "D:\tmp\main.dta", clear
*rename  contrib /*ciiu1*/ rep
*rename dni /*persona */ foreign
*rename raw_ubigeo /*condicion*/ turn
*drop if missing(rep+foreign+turn)
sysuse auto
drop if missing(rep)

qui adopath + "D:\Github\reghdfe\source\_hdfe"
qui adopath + "D:\Github\reghdfe\source\_common"
include reghdfe.mata

* Test mapsolve init
mata:
	void function testit(struct MapProblem scalar S) {
		mapsolve_precompute(S, tokens("foreign gear_ratio displacement") )
		101010
		S.verbose
		S.G
		S.fes[1].levels
		S.fes[1].inv_xx
		quadcross(S.fes[1].x, S.w)
		// liststruct(S)
		S.fes[1].cvars
		1337
		// S.fes[2].x
		mapsolve_project(S, 1, st_data(., tokens("price u") ))
	}
end

//local absvars rep#foreign   turn#rep##c.(weight gear) // B=turn (c.gear c.weight)#turn
//local absvars "foreign##c.weight"
//local absvars rep#foreign#turn
local absvars rep#foreign##c.(gear displacement) turn##c.weight trunk // rep#trunk turn#trunk trunk#foreign

set rmsg off
timer clear
sort rep foreign, stable
gen www = 1 // + int(head)
gen u = uniform()
ParseAbsvars `absvars', clustervars(`clustervars')
mata: S = mapsolve_init(2, "www")
set rmsg on
egen rf = group(rep foreign)
mata: testit(S)

reg price rf##c.(gear_ratio displacement)
predict double xb, xb
br xb

exit
cou
tab1 __ID*__ [fw=www], m

exit

* Check
sysuse auto, clear
drop if missing(rep)

*global absvars1 rep foreign
*global absvars2 turn rep

global absvars1 rep foreign
global absvars2 rep trunk
global absvars3 turn trunk
global absvars4 trunk foreign
sort rep foreign
set rmsg on
NaiveDropSingletons 4
set rmsg off
tab1 id*, m



*GenerateID rep foreign turn, gen(R1)
*assert FE1==R1

exit
/*
// -------------------------------------------------------------------------------------------------


assert altid==R1

drop FE*
mata: mapsolve_set(S)

set rmsg off
// mata: liststruct(S)
// mata: testit()







exit
mata:
mata set matastrict on

																																																																																																							
function borrar() {
	`Problem' x // or transmorphic
	x = mapsolve_init()
	x.foobar + x.spam
}	

function lee_absvars(string vector cvars) {
	4321
	`Matrix' X
	`Vector' p
	// Falta meter weights..
	X = st_data(., cvars)
	p = order(X, 1..cols(X))
	// _collate(X, p)
	return(X)
}




end

clear
sysuse auto
gen byte ok = rep >2
tab ok, m
drop if rep==.
mata: borrar()

set more off
mata: z = lee_absvars(tokens("foreign rep"))
mata: w = st_data(.,"ok")
mata: mm_freq(z,w,levels=.)
mata: levels
mata: z

mata: st_view(ivars=., ., tokens("foreign rep"))
mata: p = order(ivars, 1..cols(ivars))
mata: id = J(rows(ivars), 1, 0)
// mata: id = id + 

exit

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

# Al inicio de cada variable:

input y
y = y[p1]

# En cada loop:

for j=1:levels(FE1)
	start = ..
	end =  ..
	projection[j] = sum(y(start:end))
next

projection = projection / count_fe1
projection = repeat(projection, countfe1)

Para el FE2, repetir pero hacer antes
y = y[p2] donde p2 no va de 0 a 2 sino de 1 a 2!


es decir que para cada FE necesito siempre una variable mas..


al inicio de todo que necesito?
1) calcular promedios por grupo
2) matar singletons
3) calcular niveles

para que necesitaria pasar los niveles a la data?
- al devolver los FEs (y x lo tanto para hacer predicts)





levels = _mm_uniqrows
para generar mm_freq en algun momento genere el ID

weights as first-class citizens

me basta con generar el map ivar -> ID una vez

tengo
F = ivar1,ivar2,..,ivarN,id,freq (id es tacito xq es row #)
p (permutacion)

notar que puedo hacer
gen id1 = 0
mata: id1[invorder(p)] = 1..rows(freq_table)

WAIT:
si la tabla de freq esta sorteada, lo que tengo q hacer es (freq[col=index])[invorder(p)]


supon que hago el MM-freq o algo asi con weights

cuando un grupo tenga freq==1, me lo bajo
que pasa si trato a los singletons como weight=0 ??!?



entender que puedo meter assignments dentro de evals.. 

order(X,idx)
_collate(X,p) = X[p, .]
_sort usa collate

invorder
revorder

Y = X[p, .]
X = Y[invorder(p), .]

En matrix notation, llama P a la permutation matrix (una I trucha)
Entonces, Y = PX , X = inv(P) Y

Codigo:
p2[p] = (r>1 ? 1::r : 1..c)
Lo unico que hace



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



exit




capture program drop AltID
program define AltID, sortpreserve
	bys `0': gen long altid = (_n==1)
	replace altid = sum(altid)
	compress altid
	gen long index = _n
end







//	// Generate ID Variables
//	for (g=1;g<=G;g++) {
//		ivar_names = 
//		ivars =  // Profiler: 1%
//		N = rows(ivars)
//		K = cols(ivars)
//
//		// no puedo hacer esto xq de repente el singleton esta a la mitad del pata
//		// id = id :& singleton //
//
//		// singleton = singleton :& 
//
//		// me lo bajo asi!! dropit[invorder(p)]
//
//
//
//
//		id = runningsum(id) // meterle un mask a lo :* (!singleton) asi id=0 significa signleton
//
//		
//		// fe.ivars
//		// id, ivars
//		//ivars
//	}

	// Ya tengo los IDs
	// Ya los exporte

	// Que quiero hacer ahora?
	// 	Meter weights es facil pero hacerlo al final cuando lo necesite
	// Quiero quedarme con los IDs que no son singletons
	// hacerlo iterativamente z
