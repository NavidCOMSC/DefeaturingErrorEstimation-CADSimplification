/* NumProc to calculate upper bound error of primal problem
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

class NumProcNU : public NumProc
{
protected:
  // grid function provides the solution vector of the primal defeatured PDE
  GridFunction * gfphibar;

  // solution vector for the restricted domain in the feature
  GridFunction * gfphibar_f;

  // bilinear form (for primal defeatured)
  BilinearForm * bf;

  // Feature domain
  int domain_interest;

  // Factor (times eps0)
  double factor;

  // Verbose output
  bool verbose;

public:

  NumProcNU (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the primal and dual gridfunction as
    // -primal=u -dual=phi -feature=N [-verbose]

    gfphibar = pde.GetGridFunction (flags.GetStringFlag ("primal", "phibar"));
    gfphibar_f = pde.GetGridFunction (flags.GetStringFlag ("primal.rest", "phibar_f"));
    bf = pde.GetBilinearForm(flags.GetStringFlag ("bfdfp", "a_integral"));
    factor = flags.GetNumFlag("factor", 1);
    domain_interest = flags.GetNumFlag ("domain", 1) - 1; // -1 as outside domain is 0, so domain 1 is really domain 0, etc.
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
    cout << "Compute Upper Bound Error Element For Primal Problem Bounded to Feature Domain " << domain_interest << endl;

    if (verbose) {
      cout << "Primal GridFunction info:\n";
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

      // Setup phibar and psibar vector bound to feature domain
      FlatVector<double> phibar_f_vec = gfphibar_f->GetVector () . FV<double> ();

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
        if (eldom == domain_interest) {
          feature_ndof += fel.GetNDof ();
          n_feature++;
        }
        for (int l = 0; l < dnums.Size ();++l) {
          if (dnums[l] != -1) {
            phibar_f_vec[dnums[l]] = (eldom == domain_interest) ? gfphibar->GetVector () . FV<double> () [dnums[l]] : 0;
          }
        }
      }

      BaseVector *Kphibar_base = bf->CreateVector ();
      FlatVector<double> Kphibar_vec = Kphibar_base->FV<double> ();

      bf->ApplyMatrix (gfphibar_f->GetVector (), *Kphibar_base);

      double res = 0.0;
      for (int l = 0; l < phibar_f_vec.Size (); ++l) {
        res += phibar_f_vec[l] * Kphibar_vec[l]; // / e0
      }
      Res = factor * res; // factor should not include e0
      cout << "DOFs: " << feature_ndof << " DOFs of " << ndof << " in " << n_feature << " features of domain\n";
      cout << "NuPhi^2: " << Res << "\n";
    } else {
      cout << "ERROR: grid function not defined on Domain " << domain_interest << "\n";
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
	<< "Domain    = " << domain_interest << endl
	<< "Primal    = " << gfphibar->GetName() << endl ;
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


namespace nu_cpp
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("nu_phi", NumProcNU::Create, NumProcNU::PrintDoc);
  }

  Init init;
}

