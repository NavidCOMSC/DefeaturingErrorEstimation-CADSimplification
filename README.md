# DefeaturingErrorEstimation-CADSimplification
This code is to estimate the simplification error. It’s an an extension to ngsolve / netgen to estimate the effect of model simplifications on a quantity of interest aimed at electrostatics problems.

## CAD Model Simplification Error Estimation for Electrostatic Problems
Simplifying the geometry of a CAD model using defeaturing techniques enables more efficient discretization and subsequent simulation for engineering analysis problems. Understanding the effect this simplification has on the solution helps one to decide whether the simplification is suitable for a specific simulation problem. It can also help one to understand the functional effect of a geometry feature. The effect of the simplification is quantified by a user-defined quantity of interest which is assumed to be (approximately) linear in the solution. A bound on the difference between the quantity of interest of the original and simplified solutions based on the energy norm is derived. The approach is presented in the context of electrostatics problems but can be applied in general to a range of elliptic partial differential equations. Numerical results on the efficiency of the bound are provided for electrostatics problems with simplifications involving changes inside the problem domain as well as changes to the boundaries.

N. Rahimi, P. Kerfriden, F. C. Langbein, R.R. Martin. CAD Model simplification error estimation for electrostatics problems. SIAM J. Sci. Comput. 40(1):B196–B227, 2018 [pdf](https://epubs.siam.org/doi/10.1137/16M1078641)

![SimpErrES-4 25](https://user-images.githubusercontent.com/8395231/112233930-455ee600-8c33-11eb-8da9-3fa3072f3628.png)
