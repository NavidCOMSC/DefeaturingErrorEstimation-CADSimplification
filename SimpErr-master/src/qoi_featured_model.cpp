/* NumProc to calculate the Energy Form in Quantity Of INterest for Featured Model
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

class NumProcQOIFeaturedModel : public NumProc
{
protected:
  // grid function provides the solution vector of the featured primal PDE
  GridFunction * gfphi;

  // grid function provides the solutoin vector of the defeatured primal PDE
  //GridFunction * gfphibar;

  // solution vector for the restricted domain in the feature
  //GridFunction * gfphibar_S;

  // solution vector for the restricted domain in the feature
  GridFunction * gfphi_S;

  // Difference of featured and defeatured
  //GridFunction * diff;

  // bilinear form for featured primal
  BilinearForm * bfa;

  // bilinear form for defeatured primal
  //BilinearForm * bfabar;

  // QOI domain
  int domain_interest;

  // Verbose output
  bool verbose;

public:

  NumProcQOIFeaturedModel (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the primal and dual gridfunction as
    // -primal=u -dual=phi -feature=N [-verbose]

    gfphi = pde.GetGridFunction (flags.GetStringFlag ("featured", "phi"));
    //gfphibar = pde.GetGridFunction (flags.GetStringFlag ("defeatured", "phibar"));
    //gfphibar_S = pde.GetGridFunction (flags.GetStringFlag ("Defeatured.rest", "phibar_S"));
    gfphi_S = pde.GetGridFunction (flags.GetStringFlag ("Featured.rest", "phi_S"));
    //diff = pde.GetGridFunction (flags.GetStringFlag ("error", "diff"));
    bfa = pde.GetBilinearForm(flags.GetStringFlag ("pbf", "a"));
    //bfabar = pde.GetBilinearForm(flags.GetStringFlag ("pbfdf", "abar"));
    domain_interest = flags.GetNumFlag ("domain", 1) - 1; // -1 as outside domain is 0, so domain 1 is really domain 0, etc.
    verbose = flags.GetDefineFlag("verbose");
  }

  virtual ~NumProcQOIFeaturedModel()
  { ; }

  // creates a solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcQOIFeaturedModel (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "Compute the Energy in Quantity of Interest Domian for Featured Model " << domain_interest << endl;

    if (verbose) {
      cout << "Featured GridFunction info:\n";
      gfphi->PrintReport (cout);
      //cout << "Defeatured GridFunction info:\n";
      //gfphibar->PrintReport (cout);
    }

    FESpace *fes = const_cast<FESpace*> (&gfphi->GetFESpace ()); //gfphibar should have the same space
    if (verbose) {
      cout << "FESpace info:\n";
      fes->PrintReport (cout);
    }


    if (fes->DefinedOn (domain_interest)) {
      if (verbose) {
        cout << "Is defined on Domain " << domain_interest << "\n";
      }


      // Setup phi vector
      //FlatVector<double> phibar_S_vec = gfphibar_S->GetVector () . FV<double> ();

      FlatVector<double> phi_S_vec = gfphi_S->GetVector () . FV<double> ();

      //FlatVector<double> error_vec = diff->GetVector () . FV<double> ();

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
             phi_S_vec[dnums[l]] = (eldom == domain_interest) ? gfphi->GetVector () . FV<double> () [dnums[l]] : 0;
             //phibar_S_vec[dnums[l]] = (eldom == domain_interest) ? gfphibar->GetVector () . FV<double> () [dnums[l]] : 0;
           }
        }
      }

      //for (int l = 0; l < error_vec.Size (); ++l){
        //error_vec[l] = phi_S_vec[l] - phibar_S_vec[l];
      //}

      BaseVector *Kphi_S_base = bfa->CreateVector ();
      FlatVector<double> Kphi_S_vec = Kphi_S_base->FV<double> ();
      bfa->ApplyMatrix (gfphi_S->GetVector (), *Kphi_S_base);

      double res = 0.0;
      for (int l = 0; l < phi_S_vec.Size (); ++l) {
        res +=  phi_S_vec[l] * Kphi_S_vec[l];
      }
      cout << "DOFs: " << feature_ndof << " DOFs of " << ndof << " in " << n_feature << " features of domain\n";
      std::printf ("Q(phi) : %g\n", res);
    } else {
      cout << "ERROR: grid function not defined on Domain " << domain_interest << "\n";
    }

    cout << "End of Computation for Featured Energy in Quantity of Interest" << endl;
  }


  virtual string GetClassName () const
  {
    return "Featured Energy in Quantity Of Interest";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl
  << "Domain    = " << domain_interest << endl
  << "Featured  = " << gfphi->GetName() << endl;
  //<< "Defeatured = " << gfphibar->GetName() << endl ;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << "\n\nQOI:\n"
        << "------------------------\n"
        << "Calculate subtraction of quantity of interest for given feature domain\n"
        << "Required flags:\n"
        << "-featured=<gfname>\n"
        << "    featured grid-function\n"
        << "-feature=<domain-number>\n"
        << "    number of domain from geometry file of feature\n"
        << "-bfa=<bilinear-form>\n"
        << "    bilinear form for featured grid function\n"
        << "-verbose\n"
        << "    verbose output\n"
        << endl;
  }
};


namespace qoi_featured_model_cpp
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("qoi_featured_model", NumProcQOIFeaturedModel::Create, NumProcQOIFeaturedModel::PrintDoc);
  }

  Init init;
}

