/* Linear form
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

#ifndef FILE_MYLF_HPP
#define FILE_MYLF_HPP

namespace ngfem
{

  // integrator for \int f v dx
  class MyLFIntegrator : public LinearFormIntegrator
  {
  protected:

    //GridFunction * gfKbar;

    CoefficientFunction * coef_f;
  public:
    //gfKbar = pde.GetGridFunction (flags.GetStringFlag ("A_phibar", "Kbar"));

    MyLFIntegrator (CoefficientFunction * acoef)
      : coef_f(acoef)
    { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new MyLFIntegrator (coeffs[0]);
    }

    virtual ~MyLFIntegrator ()  { ; }

    virtual int DimElement () const { return 2; }
    virtual int DimSpace () const { return 2; }

    virtual string Name () const { return "MyLF"; }

    virtual bool BoundaryForm () const { return false; }

    // Calculates the right hand side element vector
    virtual void
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & eltrans,
		       FlatVector<double> & elvec,
		       LocalHeap & lh) const;

    /*
    AssembleElementVector (const FiniteElement & fel,
			   const ElementTransformation & eltrans,
			   FlatVector<double> & elvec,
			   LocalHeap & lh) const;
*/
  };


}

#endif