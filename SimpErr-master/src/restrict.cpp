/* NumProc to calculate the Energy Form of the subtraction Quantity Of INterest
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

class NumProcRestrict : public NumProc
{
protected:
  GridFunction * in;
  GridFunction * out;
  int domain;

public:

  NumProcRestrict (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    in = pde.GetGridFunction (flags.GetStringFlag ("in", "phi"));
    out = pde.GetGridFunction (flags.GetStringFlag ("out", "phiout"));
    domain = flags.GetNumFlag ("domain", 1) - 1; // -1 as outside domain is 0, so domain 1 is really domain 0, etc.
  }

  virtual ~NumProcRestrict ()
  { ; }

  // creates a solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcRestrict (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "Restrict Grid Function to domain " << domain << endl;

    FESpace *fes = const_cast<FESpace*> (&in->GetFESpace ()); //gfphibar should have the same space
    if (fes->DefinedOn (domain)) {
      // Setup error, phi and phibar vector
      FlatVector<double> in_vec = in->GetVector () . FV<double> ();
      FlatVector<double> out_vec = out->GetVector () . FV<double> ();

      // Find domain in the mesh
      // Get MeshAccess object
      MeshAccess & ma = pde.GetMeshAccess ();

     // Number of Elements in mesh
      int ne = ma.GetNE ();
      // Degree of freedoms array for DoFs of individual mesh element.
      Array<int> dnums;
      // Check domains of all elements
      for (int l = 0; l < ne; ++l) {
        lh.CleanUp ();
        const FiniteElement & fel = fes->GetFE (l, lh);
        fes->GetDofNrs (l, dnums);
        // Get Element domain:
        int eldom = ma.GetElIndex (l); // Domain of element l
        for (int l = 0; l < dnums.Size ();++l) {
          if (dnums[l] != -1) {
            out_vec[dnums[l]] = (eldom == domain) ? in_vec[dnums[l]] : 0;
          }
        }
      }
    } else {
      cout << "ERROR: grid function not defined on Domain " << domain << "\n";
    }

    cout << "End of Restrict" << endl;
  }


  virtual string GetClassName () const
  {
    return "Bounded Energy in Quantity Of Interest";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl;
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


namespace qoi_restrict
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("restrict", NumProcRestrict::Create, NumProcRestrict::PrintDoc);
  }

  Init init;
}

