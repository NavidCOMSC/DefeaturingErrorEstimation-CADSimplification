/* NumProc to calculate Residual Error for Positive Feature
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

class NumProcResidualPosError : public NumProc
{
protected:
  // grid function provides the solution vector of the primal defeatured PDE
  GridFunction * gfphibar;

  // grid function provides the solution vector of the primal defeatured PDE
  GridFunction * gfpsibar;

  // solution vector for the restricted domain in the feature
  GridFunction * gfphibar_f;

  // solution vector for the restricted domain in the feature
  GridFunction * gfpsibar_f;

  // solution vector for the restricted domain in the interface between feature and capacitor
  GridFunction * gfphibar_I;

  // solution vector for the restricted domain in the interface between feature and capacitor
  GridFunction * gfpsibar_I;

  // solution vector for the restricted domain in the Capacitor
  GridFunction * gfphibar_C;

  // solution vector for the restricted domain in the Capacitor
  GridFunction * gfpsibar_C;

  // bilinear form for primal defeatured
  BilinearForm * bf;

  // Feature domain
  int domain_interest;

  //Capacitor domain
  int domain_capacitor;

  // Factor (times eps0)
  double factor;

  // Verbose output
  bool verbose;

public:

  NumProcResidualPosError (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the primal and dual gridfunction as
    // -primal=u -dual=phi -feature=N [-verbose]

    gfphibar = pde.GetGridFunction (flags.GetStringFlag ("primal", "phibar"));
    gfpsibar = pde.GetGridFunction (flags.GetStringFlag ("dual", "psibar"));
    gfphibar_f = pde.GetGridFunction (flags.GetStringFlag ("primal.rest", "phibar_f"));
    gfpsibar_f = pde.GetGridFunction (flags.GetStringFlag ("dual.rest", "psibar_f"));
    gfphibar_I = pde.GetGridFunction (flags.GetStringFlag ("interface", "phibar_I"));
    gfpsibar_I = pde.GetGridFunction (flags.GetStringFlag ("interface", "psibar_I"));
    gfphibar_C = pde.GetGridFunction (flags.GetStringFlag ("primal.Cap", "phibar_C"));
    gfpsibar_C = pde.GetGridFunction (flags.GetStringFlag ("dual.Cap", "psibar_C"));
    bf = pde.GetBilinearForm(flags.GetStringFlag ("bf", "a_inetgral"));
    factor = flags.GetNumFlag("factor", 1);
    domain_interest = flags.GetNumFlag ("domain", 1) - 1; // -1 as outside domain is 0, so domain 1 is really domain 0, etc.
    domain_capacitor = flags.GetNumFlag ("capacitor", 1) - 1;
    verbose = flags.GetDefineFlag("verbose");
  }

  virtual ~NumProcResidualPosError()
  { ; }

  // creates a solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcResidualPosError (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "Compute Residual Error Problem Bounded to Positive Feature Domain " << domain_interest << endl;

    if (verbose) {
      cout << "Dual GridFunction info:\n";
      gfpsibar->PrintReport (cout);
    }

    FESpace *fes = const_cast<FESpace*> (&gfphibar->GetFESpace ());
    if (verbose) {
      cout << "FESpace info:\n";
      fes->PrintReport (cout);
    }

    if (fes->DefinedOn (domain_interest) /*&& fes->DefinedOn(domain_capacitor)*/) {
      if (verbose) {
        cout << "Is defined on Domain " << domain_interest << "\n";
      }

      // Setup phibar and psibar vector bound to feature domain
      FlatVector<double> phibar_I_vec = gfphibar_I->GetVector () . FV<double> ();
      FlatVector<double> phibar_f_vec = gfphibar_f->GetVector () . FV<double> ();
      FlatVector<double> phibar_C_vec = gfphibar_C->GetVector () . FV<double> ();
      FlatVector<double> psibar_I_vec = gfpsibar_I->GetVector () . FV<double> ();
      FlatVector<double> psibar_f_vec = gfpsibar_f->GetVector () . FV<double> ();
      FlatVector<double> psibar_C_vec = gfpsibar_C->GetVector () . FV<double> ();

      // Find domain in the mesh
      // Get MeshAccess object
      MeshAccess & ma = pde.GetMeshAccess ();

      // DOFs
      int ndof = 0;
      int feature_ndof = 0;
      int n_feature = 0;
     // Number of Elements in mesh
      int ne = ma.GetNE ();
      //int bcnr = ma.GetBCNumBCName ();
      // Degree of freedoms array for DoFs of individual mesh element.
      Array<int> dnums;
      //Final NU result
      double Res = 0.0;
      // Check domains of all elements
      for (int i = 0; i < ne; ++i) {
        lh.CleanUp ();
        const FiniteElement & fel = fes->GetFE (i, lh);
        ndof += fel.GetNDof ();
        fes->GetDofNrs (i, dnums);

        // Get Element domain:
        int eldom = ma.GetElIndex (i); // Domain of element l
        if (eldom == domain_interest ) {
          //if (eldom == domain_capacitor) {
            feature_ndof += fel.GetNDof ();
            n_feature++;
          //}
        }
        for (int l = 0; l < dnums.Size ();++l) {
          if (dnums[l] != -1) {
            phibar_I_vec[dnums[l]] = 0.0;
            psibar_I_vec[dnums[l]] = 0.0;
            phibar_f_vec[dnums[l]] = (eldom == domain_interest ) ? gfphibar->GetVector () . FV<double> () [dnums[l]] : 0;
            psibar_f_vec[dnums[l]] = (eldom == domain_interest ) ? gfpsibar->GetVector () . FV<double> () [dnums[l]] : 0;
            phibar_C_vec[dnums[l]] = (eldom == domain_capacitor ) ? gfphibar->GetVector () . FV<double> () [dnums[l]] : 0;
            psibar_C_vec[dnums[l]] = (eldom == domain_capacitor ) ? gfpsibar->GetVector () . FV<double> () [dnums[l]] : 0;
          }
        }
        if (eldom == domain_interest) {
          Array<int> ednums;
          ma.GetElEdges (i, ednums);
          Array<int> neighbours;
          ma.GetEdgeElements (i, neighbours);
          for (int ll = 0; ll < neighbours.Size (); ++ll) {
            if (neighbours[ll] != i && ma.GetElIndex (neighbours[ll]) == domain_capacitor) {
              for (int l = 0; l <dnums.Size(); ++l) {
                if (dnums[l] != -1) {
                  phibar_I_vec[dnums[l]] = gfphibar->GetVector () . FV<double> () [dnums[l]];
                  psibar_I_vec[dnums[l]] = gfpsibar->GetVector () . FV<double> () [dnums[l]];
                }
              }
            }
          }
        }
      }

      BaseVector *Kphibar_base = bf->CreateVector ();
      FlatVector<double> Kphibar_vec = Kphibar_base->FV<double> ();
      //bfabar->ApplyMatrix (gfphibar_f->GetVector (), *Kphibar_base);
      bf->ApplyMatrix (gfphibar_I->GetVector (), *Kphibar_base);

      double res = 0.0;
      for (int l = 0; l < phibar_f_vec.Size (); ++l) {
        //res += gfphibar_f->GetVector () . FV<double> ()[l] * Kphibar_vec[l]; // / e0
        res += gfpsibar_I->GetVector () . FV<double> ()[l]  * Kphibar_vec[l];
      }
      Res = factor * res; // factor should not include e0
      cout << "DOFs: " << feature_ndof << " DOFs of " << ndof << " in " << n_feature << " features of domain\n";
      cout << "Residual(psi): " << Res << "\n";
    } else {
      cout << "ERROR: grid function not defined on Domain " << domain_interest << "\n";
    }

    cout << "End Compute Residual Error in Quantity of Interest for Positive Feature" << endl;
  }


  virtual string GetClassName () const
  {
    return "Residual Error in Quantity Of Interest for Positve Feature";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl
  << "Domain    = " << domain_interest << endl
  << "Capacitor    = " << domain_capacitor << endl
  << "Primal    = " << gfphibar->GetName() << endl
  << "Dual    = " << gfpsibar->GetName() << endl;
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
        << "-feature=<domain-number>\n"
        << "    number of domain from geometry file of feature\n"
        << "-bf=<bilinear-form>\n"
        << "    bilinear form for defeatured primal grid function\n"
        << "-verbose\n"
        << "    verbose output\n"
        << endl;
  }
};


namespace residual_pos_error_cpp
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("residual_pos_error", NumProcResidualPosError::Create, NumProcResidualPosError::PrintDoc);
  }

  Init init;
}

