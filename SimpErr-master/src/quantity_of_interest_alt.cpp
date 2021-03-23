/* NumProc to calculate QUantity of INterest Alt
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

class NumProcQuantityOfInterestAlt : public NumProc
{
protected:
  // grid function provides the solution vector of the first PDE
  GridFunction * gfu;

  // Area
  int domain;

  // Verbose output
  bool verbose;

public:

  NumProcQuantityOfInterestAlt (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the linear-form and the gridfunction as
    // -linearform=f -gridfunction=u

    gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));
    domain = flags.GetNumFlag ("domain", 1) - 1; // -1 as outside domain is 0, so domain 1 is really domain 0, etc.
    verbose = flags.GetDefineFlag("verbose");
  }

  virtual ~NumProcQuantityOfInterestAlt()
  { ; }


  // creates a solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcQuantityOfInterestAlt (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "Compute General Quantity of Interest in domain " << domain << endl;

    if (verbose) {
      cout << "GridFunction info:\n";
      gfu->PrintReport (cout);
    }

    FESpace *fes = const_cast<FESpace*> (&gfu->GetFESpace ());
    if (verbose) {
      cout << "FESPace info:\n";
      fes->PrintReport (cout);
    }

    const BilinearFormIntegrator *bfi = fes->GetIntegrator ();
    if (verbose) {
      cout << "GetEvaluator (BilinearFormIntegrator):\n";
      //bfi->PrintReport (cout);
    }

    if (fes->DefinedOn (domain)) {
      if (verbose) {
        cout << "Is defined on Domain " << domain << "\n";
      }
      FiniteElement *e = const_cast<FiniteElement*> (&fes->GetFE (0, lh));
    } else {
      cout << "ERROR: grid function not defined on Domain " << domain << "\n";
    }

    // QOI Energy
    double qoi_energy = 0.0;
    double total_energy = 0.0;

    // Find domain in the mesh
    // Get MeshAccess object
    MeshAccess & ma = pde.GetMeshAccess ();
    // Number of Elements in mesh
    int ne = ma.GetNE ();
    // Degree of freedoms array for DoFs of individual mesh element.
    Array<int> dnums;

    bool hasbound = false;
    bool hasinner = false;

    for (int j = 0; j < NumIntegrators(); j++)
      {
	const BilinearFormIntegrator & bfi = *GetIntegrator(j);
	if (bfi.BoundaryForm())
	  hasbound = true;
	else
	  hasinner = true;
      }

    if (hasinner)
    // Check domains of all elements
    for (int l = 0; l < ne; ++l) {
      // Get Element domain:
      int eldom = ma.GetElIndex (l); // Domain of element l
      if (verbose) {
        std::cout << l << ": " << eldom;
        if (eldom == domain) {
          std::cout << " [QOI domain]";
        }
      }
      // Get Finite Element
      lh.CleanUp ();
      const FiniteElement & fel = fes->GetFE (l, lh);
      // ma.GetElementTransformation (l, eltrans, lh);
      ElementTransformation & eltrans = ma.GetTrafo (l, false, lh);
      fes->GetDofNrs (l, dnums);
      FlatVector<double> elvecx (dnums.Size() * fes->GetDimension(), lh);
      gfu->GetVector().GetIndirect (dnums, elvecx);
      fes->TransformVec (l, false, elvecx, TRANSFORM_SOL);

      for (int j = 0; j < NumIntegrators(); j++)
      {
	const BilinearFormIntegrator & bfi = *parts[j];

	if (bfi.BoundaryForm()) continue;
	// ...and calculate Energy
	double el_energy = bfi->Energy (fel, eltrans, elvecx, lh);

      }
      // Or call CalcFlux here? (TODO?)
      if (verbose) {
        cout << " Energy: " << el_energy << "\n";
      }
      if (eldom == domain) {
        qoi_energy += el_energy;
      }
      total_energy += el_energy;
    }

    int nse = ma.GetNSE();
    if (hasbound)
      // Check domains of all elements
    for (int l = 0; l < ne; ++l) {
      // Get Element domain:
      int eldom = ma.GetElIndex (l); // Domain of element l
      if (verbose) {
        std::cout << l << ": " << eldom;
        if (eldom == domain) {
          std::cout << " [QOI domain]";
        }
      }
      // Get Finite Element
      lh.CleanUp ();
      const FiniteElement & fel = fes->GetFE (l, lh);
      // ma.GetElementTransformation (l, eltrans, lh);
      ElementTransformation & eltrans = ma.GetTrafo (l, true, lh);
      fes->GetDofNrs (l, dnums);
      FlatVector<double> elvecx (dnums.Size() * fes->GetDimension(), lh);
      gfu->GetVector().GetIndirect (dnums, elvecx);
      fes->TransformVec (l, true, elvecx, TRANSFORM_SOL);

      for (int j = 0; j < NumIntegrators(); j++)
      {
	const BilinearFormIntegrator & bfi = *parts[j];

	if (bfi.BoundaryForm()) continue;
	// ...and calculate Energy
	double el_energy = bfi->Energy (fel, eltrans, elvecx, lh);

      }
      // Or call CalcFlux here? (TODO?)
      if (verbose) {
        cout << " Energy: " << el_energy << "\n";
      }
      if (eldom == domain) {
        qoi_energy += el_energy;
      }
      total_energy += el_energy;
    }

    cout << "Energy Total: " << total_energy << "\n";
    cout << "QOI Energy Total: " << qoi_energy << "\n";

    cout << "End Compute Quantity of Interest" << endl;
  }


  virtual string GetClassName () const
  {
    return "Quantity Of Interest";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl
	<< "Domain         = " << domain << endl
	<< "Gridfunction    = " << gfu->GetName() << endl;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << "\n\nQOI:\n"
        << "------------------------\n"
        << "Calculate quantity of interest as integral of gridfunction over given domain\n"
        << "Required flags:\n"
        << "-gridfunction=<gfname>\n"
        << "    grid-function to integrate over\n"
        << "-domain=<domain-number>\n"
        << "    number of domain from geometry file for integration\n"
        << "-verbose\n"
        << "    verbose output\n"
        << endl;
  }
};


namespace qoi_cpp
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("quantity_of_interest_alt", NumProcQuantityOfInterestAlt::Create, NumProcQuantityOfInterestAlt::PrintDoc);
  }

  Init init;
}

