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

FixStyle(sph/constant/cP,FixSPHConstantcP)

#else

#ifndef LMP_FIX_SPH_CONSTANT_cP_H
#define LMP_FIX_SPH_CONSTANT_cP_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixSPHConstantcP : public Fix {
  public:
    FixSPHConstantcP(class LAMMPS *, int, char **);
    int setmask();
    virtual void initial_integrate(int);

  private:
    class NeighList *list;

  protected:
    double *cP, *dcP;
    // Concentration to keep at constant
    double constcP;

    class Pair *pair;
  };

}

#endif
#endif
