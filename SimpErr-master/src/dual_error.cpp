/* NumProc to calculate dual error problem
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

class NumProcDualError : public NumProc
{
protected:
  // grid function provides the solution vector of the featured PDE
  GridFunction * gfpsi;

  // grid function provides the solutoin vector of the defeatured PDE
  GridFunction * gfpsibar;

  // Solution for linear form
  //GridFunction * gfL;

  // Difference of featured and defeatured
  GridFunction * diff;

  // bilinear form for featured
  BilinearForm * bfA;

  //bilinear form for defeatured
  BilinearForm * bfAbar;

  // Feature domain
  int domain_feature;

  // Verbose output
  bool verbose;

public:

  NumProcDualError (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the featured and defeatured gridfunction as
    // -featured=u -defeatured=phi -feature=N [-verbose]

    gfpsi = pde.GetGridFunction (flags.GetStringFlag ("featured", "psi"));
    gfpsibar = pde.GetGridFunction (flags.GetStringFlag ("defeatured", "psibar"));
    diff = pde.GetGridFunction (flags.GetStringFlag ("error", "diff"));
    bfA = pde.GetBilinearForm(flags.GetStringFlag ("bf", "A"));
    bfAbar = pde.GetBilinearForm(flags.GetStringFlag ("bfdf", "Abar"));
    //gfL = pde.GetGridFunction (flags.GetStringFlag ("lin", "L"));
    domain_feature = flags.GetNumFlag ("feature", 1) - 1; // -1 as outside domain is 0, so domain 1 is really domain 0, etc.
    verbose = flags.GetDefineFlag("verbose");
  }

  virtual ~NumProcDualError()
  { ; }

  // creates a solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcDualError (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {

    if (verbose) {
      cout << "Featured GridFunction info:\n";
      gfpsi->PrintReport (cout);
      cout << "Defeatured GridFunction info:\n";
      gfpsibar->PrintReport (cout);
    }

    FESpace *fes = const_cast<FESpace*> (&gfpsi->GetFESpace ()); //gfpsibar should have same space
    if (verbose) {
      cout << "FESpace info:\n";
      fes->PrintReport (cout);
    }


    // Setup psi, psibar and diff vector

    //FlatVector<double> psi_vec = gfpsi->GetVector ().FV<double> ();
    //BaseVector *psi_base = bfA->CreateVector ();
      //for (int l = 0; l < gfpsi->GetVector () . FV<double> () . Size (); l++) {
        //psi_base->FV<double>  () [l] = gfpsi->GetVector (). FV<double> () [l];
      //}
    //FlatVector<double> psi_vec = psi_base->FV<double> ();
    FlatVector<double> psi_vec = gfpsi->GetVector () . FV<double> ();

    //FlatVector<double> psibar_vec = gfpsibar->GetVector ().FV<double> ();
    //BaseVector *psibar_base = bfAbar->CreateVector ();
      //for (int l = 0; l < gfpsibar->GetVector () . FV<double> () . Size (); l++) {
        //psibar_base->FV<double>  () [l] = gfpsibar->GetVector (). FV<double> () [l];
      //}
    //FlatVector<double> psibar_vec = psibar_base->FV<double> ();
    FlatVector<double> psibar_vec = gfpsibar->GetVector () . FV<double> ();

    //FlatVector<double> err_vec = gfpsi->GetVector ().FV<double> ();
    FlatVector<double> error_vec = diff->GetVector () . FV<double> ();

    BaseVector *Kerror_base = bfA->CreateVector ();
    FlatVector<double> Kerror_vec = Kerror_base->FV<double> ();

    //BaseMatrix *bf_K=  &(bfA->GetMatrix ());
    //std::printf("featured stiffness matrix: %p\n", bf_K);
    //BaseVector *bv_K = &(bf_K->AsVector ());
    //std::printf("vectorised featured stiffness: %p\n", bv_K);
    //FlatVector<double> fv_K = bv_K->FV<double> ();

    //BaseMatrix *bf_Kbar = &(bfAbar->GetMatrix ());
    //std::printf("defeatured stiffness matrix: %p\n", bf_Kbar);
    //BaseVector *bv_Kbar = &(bf_Kbar->AsVector ());
    //std::printf("vectorized defeatured stiffness: %p\n", bv_Kbar);
    //FlatVector<double> fv_Kbar = bv_Kbar->FV<double> ();

    int n = bfAbar->GetMatrix () . Width ();
    std::printf("Size of bilinear form matrix: %d\n", n);

    // Get diagonal of K and Kbar matrix by multiplying matrix with standard basis vectors
    //BaseVector *helper = bfA->CreateVector ();
    //BaseVector *out = bfA->CreateVector ();
    //BaseVector *diagA = bfA->CreateVector ();
    //BaseVector *diagAbar = bfAbar->CreateVector ();
    //for (int l = 0; l < n; ++l) {
      //for (int k = 0; k < helper->Size (); k++) {
        //helper->FVDouble () [k] = (k == l) ? 1 : 0;
      //}
      //bfA->ApplyMatrix (*helper, *out);
      //diagA->FVDouble () [l] = out->FVDouble () [l];
      //bfAbar->ApplyMatrix (*helper, *out);
      //diagAbar->FVDouble () [l] = out->FVDouble () [l];
    //}

    // DEBUG CODE
    //std::printf("K matrix:\n");
    //bf_K->Print (std::cout);
    //std::printf("A diagonal:\n");
    //diagA->Print (std::cout);
    //std::printf("Kbar matrix:\n");
    //bf_Kbar->Print (std::cout);
    //std::printf("Abar diagonal:\n");
    //diagAbar->Print (std::cout);
    //

    for (int l = 0; l < error_vec.Size (); ++l) {
      error_vec[l] = psi_vec[l] - /* (diagAbar->FVDouble ())[l] / (diagA->FVDouble ())[l] * */ psibar_vec[l];
    }
    //for( int l=0; l<err_vec.Size (); ++l) {
      //std::printf("error vector: %f\n", err_vec[l]);
    //}
    //std::printf("\n");

    //double Res = 0.0;
    //BaseVector *bpsi_vec = bfA->CreateVector ();
    //bpsi_vec->FVDouble () = psi_vec;
    //BaseVector *bpsibar_vec = bfAbar->CreateVector ();
    //bpsibar_vec->FVDouble () = psibar_vec;
    //BaseVector *berr_vec = bfAbar->CreateVector ();
    //berr_vec ->FVDouble () = err_vec;
    //BaseVector *bKerr_vec = bfA->CreateVector ();
    bfA->ApplyMatrix (diff->GetVector (), *Kerror_base);
    //FlatVector<double> Kerr_vec = bKerr_vec->FVDouble ();
    double res = 0.0;
    for (int l = 0; l < error_vec.Size (); ++l) {
       res += error_vec[l] * Kerror_vec[l];
     }
     //Res = sqrt(res);

    cout << "a(e_psi, e_psi) " << res << "\n";
    cout << "End Compute Dual Error" << endl;
  }


  virtual string GetClassName () const
  {
    return "Error of Dual Problem";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl
	<< "Domain    = " << domain_feature << endl
	<< "featured    = " << gfpsi->GetName() << endl
	<< "defeatured  = " << gfpsibar->GetName() << endl;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << "\n\nQOI:\n"
        << "------------------------\n"
        << "Calculate error of bilinear form for dual problem\n"
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


namespace dual_error_cpp
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("dual_error", NumProcDualError::Create, NumProcDualError::PrintDoc);
  }

  Init init;
}

