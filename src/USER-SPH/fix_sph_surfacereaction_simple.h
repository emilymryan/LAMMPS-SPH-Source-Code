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

FixStyle(sph/surfacereaction/simple,FixSPHSurfaceReactionSimple)

#else

#ifndef LMP_FIX_SPH_SURFACEREACTION_SIMPLE_H
#define LMP_FIX_SPH_SURFACEREACTION_SIMPLE_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixSPHSurfaceReactionSimple : public Fix {
  public:
    FixSPHSurfaceReactionSimple(class LAMMPS *, int, char **);
    virtual ~FixSPHSurfaceReactionSimple();
    int setmask();
    virtual void init();
    virtual void final_integrate();
    void reset_dt();

  private:
    class NeighList *list;
  protected:
    double dtcA;
    double *step_respa;
    double *cA, *dcA, *mA, *dmA;
    int mass_require;

    class Pair *pair;
  };

}

#endif
#endif
