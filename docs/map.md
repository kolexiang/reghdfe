# Method of Alternating Projections

As shown in e.g. Hernandez-Ramos et al (2011), the least square residuals of a variable *y* with respect to groups of regressors X1, ..., XG, denoted "M" can be achieved either with direct methods (e.g computing `inv(X'X)X'y`) or with iterated methods such as the Method of Alternating Projections (MAP).

## Projection Schemes

The authors describe three alternating projections schemes that converge into the joint projection:

1. Kaczmarz: T = M_G M_{G-1} ... M1 
2. Cimmino: T = (M_G + M_{G-1} + ... + M1) / G
3. Symmetric Kaczmarz: T = M_G M_{G-1} ... M1 M2 ... M_G

(Note that only schemes 2 and 3 result in symmetric operators.)

In all cases, lim_{k} T^k y = M y. In other words, My is the fixed point of the Ty transformation.

However, under certain conditions (related to the second eigenvalue of the underlying graph +-+-), convergence can be slow.

## Acceleration Techniques

To improve convergence, the authors discuss two acceleration schemes:

1. Steep Descent. They show that the technique proposed by Gearahrt and Koshy, and discussed by Bauschke (and similar to the algorithms used by e.g. `lfe`, `reg2hdfe`, and older versions of `reghdfe`) is closely related to the Steep Descent method.
2. Conjugate Gradient. Note that to guarantee convergence, CG requires a symmetric operator (so Kaczmarz is ruled out). There is a version of CG that doesn't require symmetric matrices, and that is considerably slower, but I haven't tried that one.

There is also a third scheme, that works slightly better than steep descent, where the scalar `t` that determines the acceleration is computed with a variant of the Aitken's delta-squared formula. However, this acceleration may get stuck or fail to converge under some conditions, which requires aditional tuning in each iteration.

## Future

Besides getting closer to the metal (rewrite in C, multicore support besides what's already given by Mata's functions), there is an interesting alternative.

There is a string of recent papers (2013-2014) that show simple ways to solve symmetric diagonally dominant (SDD) systems in nearly-linear time.

Initial papers are very elaborate, and require tools such as complicated preconditioners, ultra-sparsifiers, fast computation of low-stretch spanning trees, efficient local clustering algorithms, etc. More recent papers reduce the complexity without compromising speed.

For two sets of fixed effects, the X'X matrix is SDD, and can thus be converted into a Laplacian with a simple trick. However, for G>2, I am not entirely sure about how to proceed. It's possible that there is no easy way to convert it into a Laplacian, but we can still use some of the tools to build a good preconditioner and thus achieve a good approximate solution before switching to an exact method.

Note: See section 9 of Kelner et al (2013) for a parallel of these tools with MAP and Kaczmarz.


