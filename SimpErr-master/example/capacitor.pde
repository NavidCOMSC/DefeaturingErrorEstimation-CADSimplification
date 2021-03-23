geometry = capacitor.in2d
mesh = capacitor.vol

shared = libsimplification


# geometry = capacitor_round.in2d
# mesh = capacitor_round.vol


define constant heapsize = 100000000

# for curved elements
# define constant geometryorder = 3

# for geometric hp-mesh refinement
# define constant hpref = 5

# dielectric constant
define constant eps0 = 8.854e-12


define coefficient coef_eps
(eps0),

# Dirichlet values according to bc number in geometry file
define coefficient coef_dirichlet
0, 1, -1,


define fespace v  -type=h1ho -order=3  -dirichlet=[2,3]

define gridfunction u -fespace=v

define linearform f -fespace=v


define bilinearform a -fespace=v -symmetric 
laplace coef_eps


#################################################
define coefficient rho
1, 1

define bilinearform b -fespace=v -symmetric
mass rho
#################################################


define preconditioner c -type=direct -bilinearform=a
# define preconditioner c -type=local -bilinearform=a
# define preconditioner c -type=multigrid -bilinearform=a -smoothingsteps=1


# sets Dirichlet boundary values
numproc setvalues npsv -coefficient=coef_dirichlet -gridfunction=u -boundary

# solves system of linear equations for free dofs
numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -maxsteps=1000 

# visualization of fluxes
numproc drawflux np2 -bilinearform=a -solution=u -label=flux


# visualization dialog box setting
numproc visualization npv1 -scalarfunction=u -subdivision=2 -nolineartexture


# evaluate energy:
numproc evaluate npeval  -bilinearform=a -gridfunction=u -gridfunction2=u -outputprecision=16


#################################################
# Evaluate alt energy from bilinear form b
numproc evaluate npeval2  -bilinearform=b -gridfunction=u -gridfunction2=u -outputprecision=16

# QOI
numproc quantity_of_interest qoi1 -gridfunction=u -domain=2
#################################################




# error estimator:

define fespace verr -type=l2 -order=0
define gridfunction err -fespace=verr

numproc zzerrorestimator np3 -bilinearform=a -solution=u -error=err -minlevel=1 -filename=errest.out
numproc markelements np5 -error=err -minlevel=1 -factor=0.5




# reference solutions
# rectangular electrodes  (p=20, hpref=10)
# 1.4229280603066e-10    
# round electrodes (p=13, geom=12, hpref=5)
# 1.488633523184388e-10

