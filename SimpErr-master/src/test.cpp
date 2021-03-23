/* NumProc to calculate primal error problem
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

class NumProcTest : public NumProc
{
protected:
  // grid function
  GridFunction * gf1;
  GridFunction * gf2;
  GridFunction * gf3;

  // bilinear form
  BilinearForm * bf;

  // Verbose output
  bool verbose;

public:

  NumProcTest (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    gf1 = pde.GetGridFunction (flags.GetStringFlag ("gf1", "phi"));
    gf2 = pde.GetGridFunction (flags.GetStringFlag ("gf2", "phibar"));
    gf3 = pde.GetGridFunction (flags.GetStringFlag ("gf3", "diff"));
    bf = pde.GetBilinearForm(flags.GetStringFlag ("bf", "a"));
    verbose = flags.GetDefineFlag("verbose");
  }

  virtual ~NumProcTest()
  { ; }

  // creates a solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcTest (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "Start Test\n";
    if (verbose) {
      cout << "GridFunction info:\n";
      gf1->PrintReport (cout);
      cout << "BilinearForm info:\n";
      bf->PrintReport (cout);
    }

    FESpace *fes = const_cast<FESpace*> (&gf1->GetFESpace ()); //gfphibar should have same space
    if (verbose) {
      cout << "FESpace info:\n";
      fes->PrintReport (cout);
    }

    FlatVector<double> vec_1 = gf1->GetVector().FV<double> ();
    FlatVector<double> vec_2 = gf2->GetVector().FV<double> ();
    FlatVector<double> vec_3 = gf3->GetVector().FV<double> ();
    for (int l = 0; l < vec_1.Size (); ++l) {
      vec_3 [l] = fabs(vec_2[l] - vec_1[l]);
    }

    cout << "End Test" << endl;
  }


  virtual string GetClassName () const
  {
    return "Test Numproc";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << "Test" << endl;
  }
};


namespace test_cpp
{
  class Init
  {
  public:
    Init ();
  };

  Init::Init()
  {
    GetNumProcs().AddNumProc ("test", NumProcTest::Create, NumProcTest::PrintDoc);
  }

  Init init;
}

