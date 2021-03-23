/* NumProc to calculate the Energy As Quantity Of INterest for Defeatured Model
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

class NumProcQOIDefeaturedModel : public NumProc
{
protected:
  // grid function provides the solutoin vector of the defeatured primal PDE
  GridFunction * gfphibar;

  // solution vector for the restricted domain in the feature
  GridFunction * gfphibar_S;

  //solution vector for right hand side
  GridFunction * gfKbar;

  // bilinear form for defeatured primal
  BilinearForm * bfabar;

  // QOI domain
  int domain_interest;

  // Verbose output
  bool verbose;

public:

  NumProcQOIDefeaturedModel (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the primal and dual gridfunction as
    // -primal=u -dual=phi -feature=N [-verbose]

    gfphibar = pde.GetGridFunction (flags.GetStringFlag ("defeatured", "phibar"));
    gfphibar_S = pde.GetGridFunction (flags.GetStringFlag ("Defeatured.rest", "phibar_S"));
    gfKbar = pde.GetGridFunction (flags.GetStringFlag ("A_phibar", "Kbar"));
    bfabar = pde.GetBilinearForm(flags.GetStringFlag ("pbfdf", "abar"));
    domain_interest = flags.GetNumFlag ("domain", 1) - 1; // -1 as outside domain is 0, so domain 1 is really domain 0, etc.
    verbose = flags.GetDefineFlag("verbose");
  }

  virtual ~NumProcQOIDefeaturedModel()
  { ; }

  // creates a solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcQOIDefeaturedModel (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "Compute the Energy in Quantity of Interest Domian for Defeatured Model " << domain_interest << endl;

    if (verbose) {
      cout << "Defeatured GridFunction info:\n";
      gfphibar->PrintReport (cout);
    }

    FESpace *fes = const_cast<FESpace*> (&gfphibar->GetFESpace ());
    if (verbose) {
      cout << "FESpace info:\n";
      fes->PrintReport (cout);
    }


    if (fes->DefinedOn (domain_interest)) {
      if (verbose) {
        cout << "Is defined on Domain " << domain_interest << "\n";
      }


      // Setup phibar vector
      FlatVector<double> phibar_S_vec = gfphibar_S->GetVector () . FV<double> ();

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
        if (eldom == domain_interest) {
          feature_ndof += fel.GetNDof ();
          n_feature++;
        }
         for (int l = 0; l < dnums.Size ();++l) {
           if (dnums[l] != -1) {
             phibar_S_vec[dnums[l]] = (eldom == domain_interest) ? gfphibar->GetVector () . FV<double> () [dnums[l]] : 0;
           }
        }
      }

      BaseVector *Kbar_base = &gfKbar->GetVector ();
      FlatVector<double> Kbar_vec = gfKbar->GetVector () . FV<double> ();
      bfabar->ApplyMatrix (gfphibar_S->GetVector (), *Kbar_base);

      double res = 0.0;
      for (int l = 0; l < phibar_S_vec.Size (); ++l) {
        res +=  phibar_S_vec[l] * Kbar_vec[l];
        //std::printf("stiffness * solution vector: %g\n", Kbar_vec[l]);
      }
      cout << "DOFs: " << feature_ndof << " DOFs of " << ndof << " in " << n_feature << " features of domain\n";
      std::printf ("Q(phibar): %g\n", res);
    } else {
      cout << "ERROR: grid function not defined on Domain " << domain_interest << "\n";
    }

    cout << "End of Compute Defeatured Energy in Quantity of Interest" << endl;
  }


  virtual string GetClassName () const
  {
    return "Defeatured Energy in Quantity Of Interest";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl
  << "Domain    = " << domain_interest << endl
  << "Defeatured = " << gfphibar->GetName() << endl ;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << "\n\nQOI:\n"
        << "------------------------\n"
        << "Calculate subtraction of quantity of interest for given feature domain\n"
        << "Required flags:\n"
        << "-defeatured=<gfname>\n"
        << "    defeatured grid-function\n"
        << "-feature=<domain-number>\n"
        << "    number of domain from geometry file of feature\n"
        << "-bfabar=<bilinear-form>\n"
        << "    bilinear form for defeatured grid function\n"
        << "-verbose\n"
        << "    verbose output\n"
        << endl;
  }
};


namespace qoi_defeatured_model_cpp
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("qoi_defeatured_model", NumProcQOIDefeaturedModel::Create, NumProcQOIDefeaturedModel::PrintDoc);
  }

  Init init;
}

