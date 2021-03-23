/* NumProc to calculate the Energy Form of the bounded QUantity of INterest
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

class NumProcQOIErrorBounded : public NumProc
{
protected:
  // grid function provides the solution vector of the featured primal PDE
  GridFunction * gfphi;

  // grid function provides the solutoin vector of the defeatured primal PDE
  GridFunction * gfphibar;

  // bilinear form for featured primal
  BilinearForm * bfa;

  // bilinear form for defeatured primal
  BilinearForm * bfabar;


  // Feature domain
  int domain_feature;

  // Verbose output
  bool verbose;

public:

  NumProcQOIErrorBounded (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the primal and dual gridfunction as
    // -primal=u -dual=phi -feature=N [-verbose]

    gfphi = pde.GetGridFunction (flags.GetStringFlag ("featured", "phi"));
    gfphibar = pde.GetGridFunction (flags.GetStringFlag ("defeatured", "phibar"));
    bfa = pde.GetBilinearForm(flags.GetStringFlag ("bf", "a"));
    bfabar = pde.GetBilinearForm(flags.GetStringFlag ("bfdf", "abar"));
    domain_feature = flags.GetNumFlag ("feature", 1) - 1; // -1 as outside domain is 0, so domain 1 is really domain 0, etc.
    verbose = flags.GetDefineFlag("verbose");
  }

  virtual ~NumProcQOIErrorBounded()
  { ; }

  // creates a solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcQOIErrorBounded (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "Compute Energy in Quantity of Interest for feature domain on the Bounded Term " << domain_feature << endl;

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


    if (fes->DefinedOn (domain_feature)) {
      if (verbose) {
        cout << "Is defined on Domain " << domain_feature << "\n";
      }

      // Setup phi and phibar vector
      FlatVector<double> phi_vec = gfphi->GetVector ().FV<double> ();
      FlatVector<double> phibar_vec = gfphibar->GetVector ().FV<double> ();

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
              phi_vec[dnums[l]] = 0;
              phibar_vec[dnums[l]] = 0;
            }
          }
        }
        if (verbose) {
          std::cout << "\n";
        }
      }


      BaseVector *bv_phi = bfa->CreateVector ();
      bv_phi->FVDouble () = phi_vec;
      BaseVector *bv_phibar = bfabar->CreateVector ();
      bv_phibar->FVDouble () = phibar_vec;
      //BaseVector *bv_Kphibar = bfabar->CreateVector ();
      BaseVector *bv_Kphi = bfa->CreateVector ();
      //bfabar->ApplyMatrix (*bv_phibar, *bv_Kphibar);
      bfa->ApplyMatrix (*bv_phibar, *bv_Kphi);
      //FlatVector<double> Kphibar_vec = bv_Kphibar->FVDouble ();
      FlatVector<double> Kphi_vec = bv_Kphi->FVDouble ();

      for (int z = 0; z < Kphi_vec.Size (); ++z){
        std::printf("The Multiplication Vector: %f\n", Kphi_vec[z] );
      }
      std::printf("\n");

      //int n = bfabar->GetMatrix () . Width ();
      int n = bfa->GetMatrix () . Width ();
      BaseVector *helper = bfa->CreateVector ();
      BaseVector *out = bfa->CreateVector ();
      BaseVector *diaga = bfa->CreateVector ();
      BaseVector *diagabar = bfabar->CreateVector ();
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < helper->Size (); j++) {
          helper->FVDouble () [j] = (j == i) ? 1 : 0;
        }
        bfa->ApplyMatrix (*helper, *out);
        diaga->FVDouble () [i] = out->FVDouble () [i];
        bfabar->ApplyMatrix (*helper, *out);
        diagabar->FVDouble () [i] = out->FVDouble () [i];
      }

      double res = 0.0;
      for (int l = 0; l < phi_vec.Size (); ++l) {
        //res += (diaga->FVDouble ())[l] / (diagabar->FVDouble ())[l] * phi_vec[l] * Kphibar_vec[l];
        //res += (diaga->FVDouble ())[l] / (diagabar->FVDouble ())[l] * phi_vec[l] * Kphi_vec[l];
        //res += phi_vec[l] * Kphibar_vec[l];
        res += phi_vec[l] * Kphi_vec[l];
      }
      cout << "DOFs: " << feature_ndof << " DOFs of " << ndof << " in " << n_feature << " features of domain\n";
      cout << "Result: " << res << "\n";
    } else {
      cout << "ERROR: grid function not defined on Domain " << domain_feature << "\n";
    }

    cout << "End Compute Bounded Energy in Quantity of Interest" << endl;
  }


  virtual string GetClassName () const
  {
    return "Bounded Energy in Quantity Of Interest";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl
	<< "Domain    = " << domain_feature << endl
	<< "Featured    = " << gfphi->GetName() << endl
	<< "Defeatured      = " << gfphibar->GetName() << endl	;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << "\n\nQOI:\n"
        << "------------------------\n"
        << "Calculate bounded quantity of interest for given feature domain\n"
        << "Required flags:\n"
        << "-featured=<gfname>\n"
        << "    featured grid-function\n"
        << "-defeatured=<gfname>\n"
        << "    defeatured grid-function\n"
        << "-feature=<domain-number>\n"
        << "    number of domain from geometry file of feature\n"
        << "-bfa=<bilinear-form>\n"
        << "    bilinear form for featured grid function\n"
	<< "-bfabar=<bilinear-form>\n"
        << "    bilinear form for defeatured grid function\n"
        << "-verbose\n"
        << "    verbose output\n"
        << endl;
  }
};


namespace qoi_error_bounded_cpp
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("qoi_error_bounded", NumProcQOIErrorBounded::Create, NumProcQOIErrorBounded::PrintDoc);
  }

  Init init;
}

