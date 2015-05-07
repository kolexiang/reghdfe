


## Literature Review

### Unconstrained Optimization Techniques for the Acceleration of Alternating Projection Methods (Hernandez-Ramos et al, 2011)

General problem: Find the closest point, to a given point, in the intersection of several sets that belong to a given space H.

Note that the intersection part is key, it means we are talking about M and not P!!!!
EG: We are in R^N and you have G sets of fixed effects (each with L_g distinct groups).
We can decompose y = yhat + resid, and what we want is resid, the vector that is orthogonal to M1, ..., MG!

So the notation is exactly the opposite... P means M!

They accelerate Cimmino's version of MAP with Conj. Gradient, when applied to "large-scale Saddle Point Problems". See [2] for uses solving large+sparse linear systems.

THM 1:
MG .. T := M1 y -> M y
THM 2:
T := (MG + .. M1)/G y -> M y

MAPs are slow without acceleration


ACceleration schemes:
Dado un T,

x(k+1) = A(x(k); T)

Donde:
A(x;T) = t_x T(x) + (1-t_x) x

Es decir, un weighted avg entre T(x) y x.
Pero cual es el weight??
Sin acelerar, el weight es 1!

Medio raro.. dicen que
Si Tx=x, t=1 (aca da igual, debe haber el problema de accuracy , ver mi mu)
Sino,

t_x = [x ' (x-Tx)] / [(x-Tx)'(x-Tx)]

Notar que la formula que uso es: ans = ans - map_project(S, g, ans)

Asi que DELTA === x - Tx = map_project(..)

Pero eso seria solo cuando G=1 :( :(

old_y = y
y = T(y)
delta = y - old_y
Uso delta para convergence
Uso delta, y para t()
Uso y, old_y, t() para A()


delta = 
cross(y, delta) / cross(delta, delta)


OJO: ESTE METODO DE ACELERACION LLAMARLO "STEEPEST DESCENT"
Tambien se puede (DEBE) reescribir mas facilmente:
NOTAR QUE ESO SIGNIFICA QUE NO NECESITO el nuevo y! Solo calculo delta directo:
Accel = x + t * delta
delta = T(y) - y
t = - cross(x, delta) / cross(delta, delta)

CG:

Para hacer CG necesitamos:
e = My = y - PX; asi que el prob (5)
f(y) = 1/2 || x - Tx ||^2
 se ve como
f(y) = 1/2 || y - Py ||^2 = e'e / 2 ; e=y-Px=Mx


gradient: derivada de e'e = (I-T)*(I-T)x = Mx ???

hessian: M ?!?
nota rque esta vaina no asume que T es self adjoint (SIMETRICO) y idempotente!
creooo q lo q quiero es accel #2!

Metodos:
kacz con no accel, steep desc, CG
cimm con *3
sym kacz con *3
es decir en total hay 9 opciones

Segun ellos el mejor es CIMMINO-CG (PERO EL MODIFICADO QUE INCORPORA EL HECHO DE QUE EL OPERADOR ES OK, ALGO 2!)



De nuevo:

Notacion: min f(x) = 1/2 || (I-T) x ||2
Debe ser || e ||
e = My = y - Py

ASI QUE EN RESUMEN, T ES P!!!!!!!!!!!!!!!!

Tambien, SELF ADJOINT = SYM

grad = My
hess = M

Notar que M es no solo SYM y IDEMP sino tambien PSD.
Prueba: yMy = e'e que es >=0

Aplicar la optimizacion directamente a My = 0

ALGO 2: NOTAR QUE ACA T es P! wtf

y = y0
r = T(y) - y
rr = r'r
u = r

for k=1..
	if convergence: return y
	v = u - T(u)
	alpha = rr / u'v
	y = y + alpha u
	r0 = rr
	r = r - alpha v
	rr = r'r
	beta = rr / r0
	u = r + beta u



minus_resid = -(y - T(y))
ssr = minus_resid' minus_resid
u = minus_resid

while ()
	
	// 
	v = u - T(u)
	alpha = ssr / u'v

	y += alpha u

	ssr_old = ssr
	minus_resid = minus_resid - alpha v
	ssr = minus_resid' minus_resid

	beta = ssr / ssr_old
	u = minus_resid + beta * u


interpretacion

	EMPIEZO CON
	yhat = 0
	resid = y - T(y)
	...

	AL FINAL QUIERO TENER ALGO ASI
	y = y - alpha u

	ESO ES LO Q TENIA CON STEEPEST DESCENT, en cuyo caso tenia
	u = y - T(y)
	alpha = x'u / u'u

	PERO SD no es optimo, QUIERO
	ANGULOS DE ATAQUE ORTOGONALES a mi anterior angulo de ataque


STRANG
	1) Compute the scaling of the orthogonal search direction d
	v = Md?? = d - T(d)
	alpha = r'r / d'v
	(esto no es algo asi como r'r / d'Md ?)

	2) Update y
	y = y + alpha d

	3) If I know the change in x, I know the change in r
	r = r - alpha v
	
	4) Beta
	beta = r'r / old_r'old_r

	5) Update orthogonal search direction
	d = r + beta * d

Cost per cycle
	v = d-T(d) // "Multiplication by A"
	Inner product r'r (guardo el anterior), inner product d'v
	2/3 vector updates






Para recuperar los FEs, me basta con cojer al resid y hacerle UNA proy wrt un FE
Es mas podria usar map_project en vez de map_solve!!!!








primero pasar las funciones a otro lado para que sean faciles de usar

luego hacer la iteracion con dos FEs,
y verificar que si repito el power iteration me sale lo mismo que con hdfe

luego hacer la super eigenaproximacion

luego de que funcione y benchmarkear, ver si pasar de 2(G-1) ops por iteracion a G ops por iteracion hace que sea mas rapido o no

mantener los dos programas al mismo tiempo???

eg: hdfe.mata 



puedo sacar el 2do eigenpair.. para q lo quiero?  pa q sirve?


el segundo eigenvalue me da una indicacion de que tan "DENSO" es el sistema??



key insight:
P is a laplacian matrix

It's a symmetric diagonally dominant matrix
