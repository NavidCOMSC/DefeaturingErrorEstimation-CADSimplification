geometry = capacitor_with_pos_feature_S.in2d
mesh = capacitor_with_pos_feature_S.vol

shared = libsimplification

define constant heapsize = 100000000

# dielectric constant
#define constant eps0 = 8.854e-12
define constant eps0 = 1.00

define coefficient coef_featured
(eps0), (eps0), (5*eps0)

define coefficient coef_defeatured
(eps0), (eps0), (eps0)


#Test Coef Integral
define coefficient coef_integral
0, 0, (eps0)

# Dirichlet values according to bc number in geometry file
define coefficient coef_dirichlet_primal
0, 220, -220, 0

# Dirichlet values according to bc number in geometry file
define coefficient coef_dirichlet_dual
0, 0, 0, 0

#Test Dirichlet Integration
define coefficient coef_dirichlet_integral
0, 0, 0, 0

define fespace v  -type=h1ho -order=3  -dirichlet=[2,3]

# Primal featured problem
#define gridfunction phi -fespace=v

#define linearform l -fespace=v

#define bilinearform a -fespace=v -symmetric
#laplace coef_featured

#define preconditioner pc_a -type=direct -bilinearform=a

# sets Dirichlet boundary values
#numproc setvalues npsv -coefficient=coef_dirichlet_primal -gridfunction=phi -boundary

# solves system of linear equations for free dofs
#numproc bvp np2 -bilinearform=a -linearform=l -gridfunction=phi -preconditioner=pc_a -maxsteps=1000


# Primal defeatured problem
define gridfunction phibar -fespace=v

define linearform lbar -fespace=v

define bilinearform abar -fespace=v -symmetric
laplace coef_defeatured

define preconditioner pc_abar -type=direct -bilinearform=abar

# sets Dirichlet boundary values
numproc setvalues npsv -coefficient=coef_dirichlet_primal -gridfunction=phibar -boundary

# solves system of linear equations for free dofs
numproc bvp np2 -bilinearform=abar -linearform=lbar -gridfunction=phibar -preconditioner=pc_abar -maxsteps=1000



# Calculate linear form for dual problems
#define gridfunction l -fespace=v

#numproc rhs calc_l -primal=phibar -lin=l -pbfdf=abar -domain=2 -verbose



# Dual Featured problem
#define gridfunction psi -fespace=v

#define linearform L -fespace=v
#source l

#define bilinearform A -fespace=v -symmetric
#laplace coef_featured

#define preconditioner pc_A -type=direct -bilinearform=A

# sets Dirichlet boundary values
#numproc setvalues npsv -coefficient=coef_dirichlet_dual -gridfunction=psi -boundary

# solves system of linear equations for free dofs
#numproc bvp np3 -bilinearform=A -linearform=L -gridfunction=psi -preconditioner=pc_A -maxsteps=1000



# Dual DeFeatured problem
define gridfunction psibar -fespace=v

define gridfunction myLF_phibar -fespace=v
numproc restrict calc_myLF -in=phibar -out=myLF_phibar -domain=2
define linearform Lbar -fespace=v
myLF myLF_phibar

define bilinearform Abar -fespace=v -symmetric
laplace coef_defeatured

define preconditioner pc_Abar -type=direct -bilinearform=Abar

# sets Dirichlet boundary values
numproc setvalues npsv -coefficient=coef_dirichlet_dual -gridfunction=psibar -boundary

# solves system of linear equations for free dofs
numproc bvp np3 -bilinearform=Abar -linearform=Lbar -gridfunction=psibar -preconditioner=pc_Abar -maxsteps=1000


# Integral bilinear form for feature
define gridfunction phi_integral -fespace=v

define linearform l_integral -fespace=v

define bilinearform a_integral -fespace=v -symmetric
laplace coef_integral

define preconditioner pc_a_integral -type=direct -bilinearform=a_integral

# sets Dirichlet boundary values
numproc setvalues npsv -coefficient=coef_dirichlet_integral -gridfunction=phi_integral -boundary

# solves system of linear equations for free dofs
numproc bvp np2 -bilinearform=a_integral -linearform=l_integral -gridfunction=phi_integral -preconditioner=pc_a_integral -maxsteps=1000


# visualization of fluxes
#numproc drawflux np2 -bilinearform=abar -solution=phibar -label=flux

# QOI
#numproc quantity_of_interest qoi1 -gridfunction=u -domain=2

define gridfunction phibar_f -fespace=v
define gridfunction phibar_I -fespace=v
define gridfunction phibar_C -fespace=v
define gridfunction psibar_f -fespace=v
define gridfunction psibar_I -fespace=v
define gridfunction psibar_C -fespace=v


# Estimate the primal and dual bounded energy in the interested area
numproc nu_pos nuepos -primal=phibar -primal.rest=phibar_f -interface=phibar_I -primal.Cap=phibar_C -bfdfp=a_integral -factor=5 -domain=3 -capacitor=1 -verbose

numproc nu_psi_pos nued -dual=psibar -dual.rest=psibar_f -interface=psibar_I -dual.Cap=psibar_C -bfdfd=a_integral -factor=5 -domain=3 -capacitor=1 -verbose

# visualization dialog box setting
numproc visualization npv1 -scalarfunction=phibar -subdivision=2 -nolineartexture
numproc visualization npv2 -scalarfunction=phibar_f -subdivision=2 -nolineartexture
numproc visualization npv2 -scalarfunction=phibar_C -subdivision=2 -nolineartexture
numproc visualization npv2 -scalarfunction=phibar_I -subdivision=2 -nolineartexture
numproc visualization npv3 -scalarfunction=psibar -subdivision=2 -nolineartexture
numproc visualization npv4 -scalarfunction=psibar_f -subdivision=2 -nolineartexture
numproc visualization npv4 -scalarfunction=psibar_I -subdivision=2 -nolineartexture
numproc visualization npv4 -scalarfunction=psibar_C -subdivision=2 -nolineartexture
#numproc visualization npv5 -scalarfunction=l -subdivision=2 -nolineartexture

