/* NumProc to calculate bilinear error for primal models
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

class NumProcBilinearErrorPrimal : public NumProc
{
protected:
  // grid function provides the solution vector of the primal PDE
  GridFunction * gfu;
  // grid function provides the solutoin vector of the dual PDE
  // GridFunction * gfphi;
  // bilinear form for gfu
  BilinearForm * bfa;

  // Feature domain
  int domain_feature;

  // Verbose output
  bool verbose;

public:

  NumProcNU (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the primal and dual gridfunction as
    // -primal=u -dual=phi -feature=N [-verbose]

    gfu = pde.GetGridFunction (flags.GetStringFlag ("primal", "u"));
    // gfphi = pde.GetGridFunction (flags.GetStringFlag ("dual", "phi"));
    bfa = pde.GetBilinearForm(flags.GetStringFlag ("bf", "a"));
    domain_feature = flags.GetNumFlag ("feature", 1) - 1; // -1 as outside domain is 0, so domain 1 is really domain 0, etc.
    verbose = flags.GetDefineFlag("verbose");
  }

  virtual ~NumProcNU()
  { ; }

  // creates a solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcNU (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "Compute Upper Bound Error Element For Primal Problem in Quantity of Interest for feature domain " << domain_feature << endl;

    if (verbose) {
      cout << "Primal GridFunction info:\n";
      gfu->PrintReport (cout);
      // cout << "Dual GridFunction info:\n";
      // gfphi->PrintReport (cout);
    }

    FESpace *fes = const_cast<FESpace*> (&gfu->GetFESpace ()); //gfphi should have same space
    if (verbose) {
      cout << "FESpace info:\n";
      fes->PrintReport (cout);
    }

    const BilinearFormIntegrator *bfi = fes->GetIntegrator ();

    if (fes->DefinedOn (domain_feature)) {
      if (verbose) {
        cout << "Is defined on Domain " << domain_feature << "\n";
      }

      // Setup u and phi vector
      FlatVector<double> u_vec = gfu->GetVector ().FV<double> ();
      // FlatVector<double> phi_vec = gfphi->GetVector ().FV<double> ();

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
        //FlatVector<double> elvecx (dnums.Size() * fes->GetDimension(), lh);
        //gfu->GetVector().GetIndirect (dnums, elvecx);

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
              u_vec[dnums[l]] = 0;
              // phi_vec[dnums[l]] = 0;
            }
          }
        }
        if (verbose) {
          std::cout << "\n";
        }
      }
        // Get Finite Element
/*
      /* FIXME - needed?
      // Boundary energy:
      if (verbose) {
        std::cout <<"Boundary energy:\n";
      }
      // Number of Boundary Elements in mesh
      int nse = ma.GetNSE ();
      // Check domains of all elements
      for (int l = 0; l < nse; ++l) {
        // Get Element domain:
        int eldom = ma.GetSElIndex (l); // Domain of element l
        if (verbose) {
          std::cout << l << ": " << eldom;
          if (eldom == domain) {
            std::cout << " [QOI domain]";
          }
        }
        // Get Finite Element
        lh.CleanUp ();
        const FiniteElement & fel = fes->GetSFE (l, lh);
        ElementTransformation & eltrans = ma.GetTrafo (l, true, lh);
        fes->GetSDofNrs (l, dnums);
        FlatVector<double> elvecx (dnums.Size() * fes->GetDimension(), lh);
        gfu->GetVector().GetIndirect (dnums, elvecx);
        fes->TransformVec (l, true, elvecx, TRANSFORM_SOL);
        // ...and calculate Energy
        if (bfi->BoundaryForm ()) {
          double el_energy = bfi->Energy (fel, eltrans, elvecx, lh);
          // Or call CalcFlux here? (TODO?)
          if (verbose) {
            cout << " Energy: " << el_energy << "\n";
          }
          if (eldom == domain) {
            sdomain_energy += el_energy;
          }
          stotal_energy += el_energy;
        }
      }

      // Special energy:
      if (verbose) {
        std::cout << "Special energy:\n";
      }
      lh.CleanUp();
      for (int l = 0; l < fes->specialelements.Size(); l++) {
        const SpecialElement *el = fes->specialelements[l];
        el->GetDofNrs (dnums);
        FlatVector<double> elvecx (dnums.Size() * fes->GetDimension(), lh);
        gfu->GetVector().GetIndirect (dnums, elvecx);
        double el_energy = el->Energy (elvecx, lh);
        if (verbose) {
          cout << " Energy: " << el_energy << "\n";
        }
        special_energy += el_energy;
      }
      lh.CleanUp();
*/
      BaseVector *bu_vec = bfa->CreateVector ();
      bu_vec->FVDouble () = u_vec;
      BaseVector *bKu_vec = bfa->CreateVector ();
      bfa->ApplyMatrix (*bu_vec, *bKu_vec);
      FlatVector<double> Ku_vec = bKu_vec->FVDouble ();
      double res = 0.0;
      for (int l = 0; l < u_vec.Size (); ++l) {
        res += u_vec[l] * Ku_vec[l];
      }
      Res = 4 * res;
      cout << "DOFs: " << feature_ndof << " DOFs of " << ndof << " in " << n_feature << " features of domain\n";
      cout << "Result: " << Res << "\n";
    } else {
      cout << "ERROR: grid function not defined on Domain " << domain_feature << "\n";
    }

    cout << "End Compute Primal Upper Bound Error in Quantity of Interest" << endl;
  }


  virtual string GetClassName () const
  {
    return "Upper Bound Error of Primal Problem in Quantity Of Interest";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl
	<< "Domain    = " << domain_feature << endl
	<< "Primal    = " << gfu->GetName() << endl ;
  // << "Dual      = " << gfphi->GetName() << endl	;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << "\n\nQOI:\n"
        << "------------------------\n"
        << "Calculate error in quantity of interest for given feature domain\n"
        << "Required flags:\n"
        << "-primal=<gfname>\n"
        << "    primal grid-function\n"
        // << "-dual=<gfname>\n"
        // << "    dual grid-function\n"
        << "-feature=<domain-number>\n"
        << "    number of domain from geometry file of feature\n"
        << "-bf=<bilinear-form>\n"
        << "    bilinear form for primal grid function\n"
        << "-verbose\n"
        << "    verbose output\n"
        << endl;
  }
};


namespace nu_cpp
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("nu", NumProcNU::Create, NumProcNU::PrintDoc);
  }

  Init init;
}