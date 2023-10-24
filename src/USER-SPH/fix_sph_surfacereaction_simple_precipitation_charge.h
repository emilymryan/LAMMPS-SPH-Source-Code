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

FixStyle(sph/surfacereaction/simple/precipitation/charge,FixSPHSurfaceReactionSimplePrecipitationCharge)

#else

#ifndef LMP_FIX_SPH_SURFACEREACTION_SIMPLE_PRECIPITATION_CHARGE_H
#define LMP_FIX_SPH_SURFACEREACTION_SIMPLE_PRECIPITATION_CHARGE_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixSPHSurfaceReactionSimplePrecipitationCharge : public Fix {
  public:
    FixSPHSurfaceReactionSimplePrecipitationCharge(class LAMMPS *, int, char **);
    virtual ~FixSPHSurfaceReactionSimplePrecipitationCharge();
    int setmask();
    virtual void init();
    virtual void init_list(int, class NeighList *);
    virtual void end_of_step();

    void unpack_reverse_comm(int, int *, double *);
    int pack_reverse_comm(int, int, double *);

  private:
    class NeighList *list;

  protected:
    int is_periodic;
    double mAthres, cAeq;
    double *step_respa;
    double *cA, *mA;
    double *ischangecA;
    int mass_require;

    class Pair *pair;
  };

}

#endif
#endif
