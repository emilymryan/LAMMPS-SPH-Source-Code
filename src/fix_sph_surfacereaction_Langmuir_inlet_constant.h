/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(sph/surfacereaction/Langmuir/inletconstant,FixSPHSurfaceReactionLangmuirInletConstant)

#else

#ifndef LMP_FIX_SPH_SURFACEREACTION_LANGMUIR_INLET_CONSTANT_H
#define LMP_FIX_SPH_SURFACEREACTION_LANGMUIR_INLET_CONSTANT_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixSPHSurfaceReactionLangmuirInletConstant : public Fix {
  public:
    FixSPHSurfaceReactionLangmuirInletConstant(class LAMMPS *, int, char **);
    int setmask();
    virtual void init();
    virtual void initial_integrate(int);

  private:
    class NeighList *list;
    // region to keep the inlet concentration
    int iregion, maxattempt, scaleflag;
    char *idregion;
    double xlo, xhi, ylo, yhi, zlo, zhi;

    void options(int, char **);

  protected:
    double *xA;
    // concentration to keep the inlet at
    double xAin;

    class Pair *pair;
  };

}

#endif
#endif
