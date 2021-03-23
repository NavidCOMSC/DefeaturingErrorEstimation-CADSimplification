/* MyLF
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

#include <fem.hpp>
#include <comp.hpp>
#include "myLF.hpp"
//#include <cassert>

using namespace ngcomp;

namespace ngfem
{

  /*
    Calculates the element vector.

    Input is:
    the finite element: fel
    the geometry of the element: eltrans

    Output is:
    the element vector: elvec

    Efficient memorymanagement is provided my locheap
  */

  /*
  void MyLFIntegrator ::
  AssembleElementVector (const FiniteElement & base_fel,
                         //const ElementTransformation & eltrans,
                         FlatVector<double> & elvec,
                         LocalHeap & lh) const
   */

  void MyLFIntegrator ::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans,
		     FlatVector<double> & elvec,
		     LocalHeap & lh) const
   {
    //GridFunctionCoefficientFunction *gfcf = dynamic_cast<GridFunctionCoefficientFunction*> (coef_f);
    //if (!gfcf) {
    //  cerr << "myLF needs a GrifFunctionCoefficientFunction\n";
    //  exit (1);
    //}
    //GridFunction const * gf = &gfcf->GetGridFunction ();

    //double epsilon = 8.854e-12; // Must be epsilon in quantity of interest domain (in defeatured primal model)
    //double epsilon = 1.00; // Must be epsilon in quantity of interest domain (in defeatured primal model)

    //elvec = 0;

    //solution vector for right hand side
    GridFunctionCoefficientFunction *gfcf = dynamic_cast<GridFunctionCoefficientFunction*> (coef_f);
    if (!gfcf) {
      cerr << "myLF needs a GrifFunctionCoefficientFunction\n";
      exit (1);
    }
    GridFunction const * gfKbar = &gfcf->GetGridFunction ();
    //FlatVector<double> Kbar_vec = gfKbar->GetVector () . FV<double> ();

    //const ScalarFiniteElement<2> & fel =
    //  dynamic_cast<const ScalarFiniteElement<2> &> (base_fel);

    //int elnr = eltrans.GetElementIndex ();
    // number of element basis functions:
    //int ndof = fel.GetNDof();


    //int ne = ma.GetNE ();
    //int ndof = 0;
    //for (int i=0; i < elnr; ++i) {
     //int ndof = fel.GetNDof ();

    /*
    Matrix<> dshape_ref(ndof, 2); // gradient on reference element
    Matrix<> dshape(ndof, 2);     // gradient on mapped element

    // Get vector elnr to elnr+ndof-1 from gf and multiply with matrix?
    FlatVector<double> gf_el (ndof, lh);
    Array<int> dnums;
    //fes.GetNodeDofNrs (nt, i,  dnums);
    gf->GetFESpace().GetDofNrs (elnr, dnums);
    gf->GetElementVector (dnums, gf_el);

    for (int l = 0; l < dnums.Size ();++l) {
             elvec[dnums[l]] = (sizeof(Kbar_vec) == sizeof(dnums)) ? gfKbar->GetVector () . FV<double> () [dnums[l]] : 0;
        }

    //for (int l = 0; l < ndof; ++l) {
    //  gf_el[l] = gf->GetVector().FV<double> () [elnr+l];
   // }

    const IntegrationRule & ir =
      SelectIntegrationRule (fel.ElementType(), 2*fel.Order());

    //Related domain for each element l
    for (int i = 0; i < ndof; ++i) {
      for (int j = 0 ; j < ir.GetNIP(); j++){
        // calculate Jacobi matrix in the integration point
        //SpecificIntegrationPoint<2,2> sip(ir[j], eltrans); //, lh);

        MappedIntegrationPoint<2,2> mip(ir[i], eltrans);
         //
        //gradient on reference element
        //the i-th row of the matrix is the grad of the i-th shape function
        //
        fel.CalcDShape (ir[j], dshape_ref);
        // transform it for the mapped element
        //dshape = dshape_ref * sip.GetJacobianInverse();

        dshape = dshape_ref * mip.GetJacobianInverse();

        // integration weight and Jacobi determinant
        //double fac = sip.IP().Weight() * fabs (sip.GetJacobiDet());

        double fac = mip.IP().Weight() * mip.GetMeasure();
        elvec += (fac * epsilon) * dshape * Trans(dshape) * gf_el ;
      }
    }
    */

    Array<int> dnums;
    int elnr = eltrans.GetElementNr ();
    gfKbar->GetFESpace().GetDofNrs (elnr, dnums);
    for (int l = 0; l < dnums.Size (); ++l) {
      elvec[l] = gfKbar->GetVector () . FV<double> () [dnums[l]];
    }

    //for (int l = 0; l < ndof; ++l) {
      //cout << "El " << elnr << ", dof: " << l << ", index: " << dnums[l] << ", Value: " << elvec[l] << "\n";
    //}

    return;
  }


  namespace init_myLF
  {
    class Init
    {
    public:
      Init ();
    };

    Init::Init()
    {
      GetIntegrators().AddLFIntegrator ("myLF", 2, 1,
                                        MyLFIntegrator::Create);
    }

    Init init;
  }
}

