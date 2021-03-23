/* NumProc to calculate the right hand side for dual problem
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

class NumProcRHS : public NumProc
{
protected:
  // grid function provides the solution vector of the primal PDE
  GridFunction * gfphibar;

  // solution vector for the restricted domain in the feature
  GridFunction * gfphibar_S;

  //solution vector for right hand side
  GridFunction * gfKbar;

  // bilinear form for primal defeatured PDE
  BilinearForm * bfabar;

  // Output for linear form
  //GridFunction * gfl;

  // QOI domain
  int domain_interest;

  // Epsilon
  //double epsilon;

  // Verbose output
  bool verbose;

public:

  NumProcRHS (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the primal and dual gridfunction as
    // -primal=phi -dual=psi -feature=N [-verbose]

    gfphibar = pde.GetGridFunction (flags.GetStringFlag ("defeatured", "phibar"));
    gfphibar_S = pde.GetGridFunction (flags.GetStringFlag ("Defeatured.rest", "phibar_S"));
    gfKbar = pde.GetGridFunction (flags.GetStringFlag ("A_phibar", "Kbar"));
    //gfl = pde.GetGridFunction (flags.GetStringFlag ("lin", "lf"));
    bfabar = pde.GetBilinearForm(flags.GetStringFlag ("pbfdf", "abar"));
    domain_interest = flags.GetNumFlag ("domain", 1) - 1; // -1 as outside domain is 0, so domain 1 is really domain 0, etc.
    //epsilon = flags.GetNumFlag ("epsilon", 1.0);
    verbose = flags.GetDefineFlag("verbose");
  }

  virtual ~NumProcRHS()
  { ; }

  // creates a solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcRHS (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "Compute Right Hand Side for Defeatured Dual " << domain_interest << endl;

    FESpace *fes = const_cast<FESpace*> (&gfphibar->GetFESpace ());
    if (verbose) {
      cout << "FESpace info:\n";
      fes->PrintReport (cout);
    }

    //if (fes->DefinedOn (domain_interest)) {
      //if (verbose) {
        //cout << "Is defined on Domain " << domain_interest << "\n";
      //}

    //FlatVector<double> l_vec = gfl->GetVector () . FV<double> ();
    FlatVector<double> phibar_S_vec = gfphibar_S->GetVector () . FV<double> ();

    //Find the Domain for QOI
    //Requirement for MeshAccess Object
    MeshAccess & ma = pde.GetMeshAccess ();

    //Degree of Freedoms
    int ndof = 0;
    int feature_ndof = 0;
    int n_feature = 0;
    //Object for number of elements in the mesh
    int ne = ma.GetNE ();
    //Degree of freedoms array for DoFs of individual mesh element
    Array<int> dnums;
    // Check each domain tagged for the each elmenet
    for (int i = 0; i < ne; ++i) {
        lh.CleanUp ();
        const FiniteElement & fel = fes->GetFE (i, lh);
        ndof += fel.GetNDof ();
        fes->GetDofNrs (i, dnums);

        //Related domain for each element l
        int eldom = ma.GetElIndex (i);
        if (eldom == domain_interest) {
          feature_ndof += fel.GetNDof ();
          n_feature++;
        }

        for (int i = 0; i < dnums.Size ();++i) {
          if (dnums[i] != -1) {
            phibar_S_vec[dnums[i]] = (eldom == domain_interest) ?  gfphibar->GetVector () . FV <double> () [dnums[i]] : 0.0;
          }
        }
    }

    BaseVector *Kbar_base = &gfKbar->GetVector ();
    FlatVector<double> Kbar_vec = gfKbar->GetVector () . FV<double> ();
    //bf->ApplyMatrix (gfphibar_S->GetVector (), *Kbar_base);
    bfabar->ApplyMatrix (gfphibar_S->GetVector (), *Kbar_base);

    //lin = 0;

    /*
    for (int i = 0; i < phibar_S_vec.Size (); ++i) {
      lin[i] = Kbar_vec[i]
    }
    */

    cout << "End Compute Right Hand Side for Defeatured Dual Problem" << endl;
  }


  virtual string GetClassName () const
  {
    return "Right Hand Side Defeatured Dual";
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
        << "Calculate upper bound error in quantity of interest for given feature domain in dual problem\n"
        << "Required flags:\n"
        << "-dual=<gfname>\n"
        << "dual grid-function\n"
        << "-feature=<domain-number>\n"
        << "number of domain from geometry file of feature\n"
        << "-bf=<bilinear-form>\n"
        << "bilinear form for primal grid function\n"
        << "-verbose\n"
        << "verbose output\n"
        << endl;
  }
};


namespace rhs_cpp
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("rhs", NumProcRHS::Create, NumProcRHS::PrintDoc);
  }

  Init init;
}

