/* NumProc to calculate primal error problem
 * Extension module for netgen/ngsolve.
 * Copyright (C) 2015-2016 Navid Rahimi <RahimiA@cardiff.ac.uk>,
 *                         Cardiff University
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <solve.hpp>

using namespace ngsolve;

class NumProcPrimalError : public NumProc
{
protected:
  // grid function provides the solution vector of the featured PDE
  GridFunction * gfphi;

  // grid function provides the solutoin vector of the defeatured PDE
  GridFunction * gfphibar;

  // Difference of featured and defeatured
  GridFunction * diff;

  // bilinear form for gfphi
  BilinearForm * bfa;

  //bilinear form for gfphibar
  BilinearForm * bfabar;

  // Feature domain
  int domain_feature;

  // Verbose output
  bool verbose;

public:

  NumProcPrimalError (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the featured and defeatured gridfunction as
    // -featured=u -defeatured=phi -feature=N [-verbose]

    gfphi = pde.GetGridFunction (flags.GetStringFlag ("featured", "phi"));
    gfphibar = pde.GetGridFunction (flags.GetStringFlag ("defeatured", "phibar"));
    diff = pde.GetGridFunction (flags.GetStringFlag ("error", "diff"));
    bfa = pde.GetBilinearForm(flags.GetStringFlag ("bff", "a"));
    bfabar = pde.GetBilinearForm(flags.GetStringFlag ("bfdf", "abar"));
    domain_feature = flags.GetNumFlag ("feature", 1) - 1; // -1 as outside domain is 0, so domain 1 is really domain 0, etc.
    verbose = flags.GetDefineFlag("verbose");
  }

  virtual ~NumProcPrimalError()
  { ; }

  // creates a solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcPrimalError (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {

    if (verbose) {
      cout << "Featured GridFunction info:\n";
      gfphi->PrintReport (cout);
      cout << "Defeatured GridFunction info:\n";
      gfphibar->PrintReport (cout);
    }

    FESpace *fes = const_cast<FESpace*> (&gfphi->GetFESpace ()); //gfphibar should have same space
    if (verbose) {
      cout << "FESpace info:\n";
      fes->PrintReport (cout);
    }


    // Setup phi, phibar and err vector

    //FlatVector<double> phi_vec = gfphi->GetVector ().FV<double> ();
    //BaseVector *phi_base = bfa->CreateVector ();
      //for (int l = 0; l < gfphi->GetVector () . FV<double> () . Size (); l++) {
        //phi_base->FV<double>  () [l] = gfphi->GetVector (). FV<double> () [l];
      //}
    //FlatVector<double> phi_vec = phi_base->FV<double> ();
    FlatVector<double> phi_vec = gfphi->GetVector () . FV<double> ();

    //FlatVector<double> phibar_vec = gfphibar->GetVector ().FV<double> ();
    //BaseVector *phibar_base = bfabar->CreateVector ();
      //for (int l = 0; l < gfphibar->GetVector () . FV<double> () . Size (); l++) {
        //phibar_base->FV<double>  () [l] = gfphibar->GetVector (). FV<double> () [l];
      //}
    //FlatVector<double> phibar_vec = phibar_base->FV<double> ();
    FlatVector<double> phibar_vec = gfphibar->GetVector () . FV<double> ();

    //FlatVector<double> err_vec = gfphi->GetVector ().FV<double> ();
    //BaseVector *error_base = bfa->CreateVector ();
    FlatVector<double> error_vec = diff->GetVector () . FV<double> ();

    BaseVector *Kerror_base = bfa->CreateVector ();
    FlatVector<double> Kerror_vec = Kerror_base->FV<double> ();

/*
      // Find domain in the mesh
      // Get MeshAccess object
      MeshAccess & ma = pde.GetMeshAccess ();

      // DOFs
      int ndof = 0;
      int feature_ndof = 0;
      int n_feature = 0;
     // Number of Elements in mesh
      int ne = ma.GetNE ();
      // Degree of freedoms array for DoFs of individual mesh element.
      Array<int> dnums;
      //Final NU result
      double Res = 0.0;
      // Check domains of all elements
      for (int l = 0; l < ne; ++l) {
        lh.CleanUp ();
        const FiniteElement & fel = fes->GetFE (l, lh);
        ndof += fel.GetNDof ();
        fes->GetDofNrs (l, dnums);

        // Get Element domain:
        int eldom = ma.GetElIndex (l); // Domain of element l
        if (verbose) {
          std::cout << l << ": " << eldom;
        }
        if (eldom == domain_feature) {
          feature_ndof += fel.GetNDof ();
          if (verbose) {
            std::cout << " [Feature domain]";
          }
          n_feature++;
        } else {
          for (int l = 0; l < dnums.Size ();++l) {
            if (dnums[l] != -1) {
              phibar_vec[dnums[l]] = 0;
            }
          }
        }
        if (verbose) {
          std::cout << "\n";
        }
      }
*/
    //BaseMatrix *bf_A=  &(bfa->GetMatrix ());
    //std::printf("featured stiffness matrix: %p\n", bf_A);
    //BaseVector *bv_A = &(bf_A->AsVector ());
    //std::printf("vectorised featured stiffness: %p\n", bv_A);
    //FlatVector<double> fv_A = bv_A->FV<double> ();

    //BaseMatrix *bf_Abar = &(bfabar->GetMatrix ());
    //std::printf("defeatured stiffness matrix: %p\n", bf_Abar);
    //BaseVector *bv_Abar = &(bf_Abar->AsVector ());
    //std::printf("vectorized defeatured stiffness: %p\n", bv_Abar);
    //FlatVector<double> fv_Abar = bv_Abar->FV<double> ();

    int n = bfabar->GetMatrix () . Width ();
    //std::printf("Size of bilinear form matrix: %d\n", n);
    // Get diagonal of A and Abar matrix by multiplying matrix with standard basis vectors
    //BaseVector *helper = bfa->CreateVector ();
    //BaseVector *out = bfa->CreateVector ();
    //BaseVector *diaga = bfa->CreateVector ();
    //BaseVector *diagabar = bfabar->CreateVector ();
    //for (int l = 0; l < n; ++l) {
      //for (int k = 0; k < helper->Size (); k++) {
        //helper->FVDouble () [k] = (k == l) ? 1 : 0;
      //}
      //bfa->ApplyMatrix (*helper, *out);
      //diaga->FVDouble () [l] = out->FVDouble () [l];
      //bfabar->ApplyMatrix (*helper, *out);
      //diagabar->FVDouble () [l] = out->FVDouble () [l];
    //}

    // DEBUG CODE
    //std::printf("A matrix:\n");
    //bf_A->Print (std::cout);
    //std::printf("A diagonal:\n");
    //diaga->Print (std::cout);
    //std::printf("Abar matrix:\n");
    //bf_Abar->Print (std::cout);
    //std::printf("Abar diagonal:\n");
    //diagabar->Print (std::cout);
    //

    for (int l = 0; l < error_vec.Size (); ++l) {
      // std::printf("error vector[%d]: phi:%g - (%g / %g) * phibar:%g\n", l, phi_vec[l], (diagabar->FVDouble ())[l], (diaga->FVDouble ())[l], phibar_vec[l] );
      // err_vec[l] = u_vec[l] - phi_vec[l];
      //std::printf ("l: %d\n", l);
      error_vec[l] = phi_vec[l] - /*(diagabar->FVDouble ())[l] / (diaga->FVDouble ())[l] * */  phibar_vec[l];
      //error_vec[l] = phi_vec[l] - (diagabar->FVDouble ())[l] / (diaga->FVDouble ())[l] *  phibar_vec[l];
      //std::printf(" = %g\n", l, error_vec[l]);
    }

    /*for (int z = 0; z < error_vec.Size (); ++z){
        std::printf("The Factor[%d]: %g / %g = %g\n", z, (diagabar->FVDouble ())[z], (diaga->FVDouble ())[z], (diagabar->FVDouble ())[z] / (diaga->FVDouble ())[z] );
      }
      std::printf("\n"); */



    //for( int l=0; l<error_vec.Size (); ++l) {
      //std::printf("error vector %g\n", error_vec[l]);
    //}
    //std::printf("\n");

    //double Res = 0.0;
    //BaseVector *bphi_vec = bfa->CreateVector ();
    //bphi_vec->FVDouble () = phi_vec;
    //BaseVector *bphibar_vec = bfb->CreateVector ();
    //bphibar_vec->FVDouble () = phibar_vec;
    //BaseVector *berr_vec = bfb->CreateVector ();
    //berr_vec ->FVDouble () = err_vec;
    //BaseVector *bKerr_vec = bfa->CreateVector ();
    //bKerr_vec->FVDouble () = u_vec; // Initialisatoin might not be needed?
    bfa->ApplyMatrix (diff->GetVector (), *Kerror_base);
    //FlatVector<double> Kerr_vec = bKerr_vec->FVDouble ();
    double res = 0.0;
    for (int l = 0; l < error_vec.Size (); ++l) {
       res += error_vec[l] * Kerror_vec[l];
     }
     //Res = sqrt(res);

    cout << "a(e_phi , e_phi): " << res << "\n";
    cout << "End Compute Primal Error" << endl;
  }


  virtual string GetClassName () const
  {
    return "Error of Primal Problem";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl
	<< "Domain    = " << domain_feature << endl
	<< "featured    = " << gfphi->GetName() << endl
  << "defeatured  = " << gfphibar->GetName() << endl;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << "\n\nQOI:\n"
        << "------------------------\n"
        << "Calculate error of bilinear form for primal problem\n"
        << "Required flags:\n"
        << "-featured=<gfname>\n"
        << " featured grid-function\n"
        << "-defeatured=<gfname>\n"
        << " defeatured grid-function\n"
        << "-feature=<domain-number>\n"
        << " number of domain from geometry file of feature\n"
        << "-bf=<bilinear-form>\n"
        << " bilinear form for featured grid function\n"
        << "-verbose\n"
        << "    verbose output\n"
        << endl;
  }
};


namespace primal_error_cpp
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("primal_error", NumProcPrimalError::Create, NumProcPrimalError::PrintDoc);
  }

  Init init;
}

