/* NumProc to calculate Quantity of Interest in Interested Domain
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

class NumProcQOIError : public NumProc
{
protected:
  // grid function provides the solution vector of the primal PDE
  GridFunction * gfphibar;

  // bilinear form for defeatured primal
  BilinearForm * bfabar;

  // Feature domain
  int domain_feature;

  // Verbose output
  bool verbose;

public:

  NumProcQOIError (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the primal and dual gridfunction as
    // -primal=u -dual=phi -feature=N [-verbose]

    gfphibar = pde.GetGridFunction (flags.GetStringFlag ("defeatured", "phibar"));
    bfabar = pde.GetBilinearForm(flags.GetStringFlag ("bfdf", "abar"));
    domain_feature = flags.GetNumFlag ("feature", 1) - 1; // -1 as outside domain is 0, so domain 1 is really domain 0, etc.
    verbose = flags.GetDefineFlag("verbose");
  }

  virtual ~NumProcQOIError()
  { ; }

  // creates a solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcQOIError (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "Compute Eenergy in the Interested Domain " << domain_feature << endl;

    if (verbose) {
      cout << "Primal GridFunction info:\n";
      gfphibar->PrintReport (cout);
    }

    FESpace *fes = const_cast<FESpace*> (&gfphibar->GetFESpace ());
    if (verbose) {
      cout << "FESpace info:\n";
      fes->PrintReport (cout);
    }


    if (fes->DefinedOn (domain_feature)) {
      if (verbose) {
        cout << "Is defined on Domain " << domain_feature << "\n";
      }

      // Setup phibar vector
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
              phibar_vec[dnums[l]] = 0;
            }
          }
        }
        if (verbose) {
          std::cout << "\n";
        }
      }

      BaseVector *bphibar_vec = bfabar->CreateVector ();
      bphibar_vec->FVDouble () = phibar_vec;
      BaseVector *bKphibar_vec = bfabar->CreateVector ();
      bfabar->ApplyMatrix (*bphibar_vec, *bKphibar_vec);
      FlatVector<double> Kphibar_vec = bKphibar_vec->FVDouble ();
      double res = 0.0;
      for (int l = 0; l < phibar_vec.Size (); ++l) {
        res += phibar_vec[l] * Kphibar_vec[l];
      }
      cout << "DOFs: " << feature_ndof << " DOFs of " << ndof << " in " << n_feature << " features of domain\n";
      cout << "Result: " << res << "\n";
    } else {
      cout << "ERROR: grid function not defined on Domain " << domain_feature << "\n";
    }

    cout << "End Compute Energy in Quantity of Interest" << endl;
  }


  virtual string GetClassName () const
  {
    return "Energy in Quantity Of Interest";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl
	<< "Domain    = " << domain_feature << endl
	<< "Primal    = " << gfphibar->GetName() << endl ;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << "\n\nQOI:\n"
        << "------------------------\n"
        << "Calculate error in quantity of interest for given feature domain\n"
        << "Required flags:\n"
        << "-Primal=<gfname>\n"
        << "    Primal grid-function\n"
        << "-feature=<domain-number>\n"
        << "    number of domain from geometry file of feature\n"
        << "-bf=<bilinear-form>\n"
        << "    bilinear form for primal grid function\n"
        << "-verbose\n"
        << "    verbose output\n"
        << endl;
  }
};


namespace qoi_error_cpp
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("qoi_error", NumProcQOIError::Create, NumProcQOIError::PrintDoc);
  }

  Init init;
}

