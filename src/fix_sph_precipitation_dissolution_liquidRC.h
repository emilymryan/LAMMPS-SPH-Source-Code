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

FixStyle(sph/precipitation/dissolution/liquidRC,FixSPHPrecipitationDissolutionLiquidRC)

#else

#ifndef LMP_FIX_SPH_PRECIPITATION_DISSOLUTION_LIQUIDRC_H
#define LMP_FIX_SPH_PRECIPITATION_DISSOLUTION_LIQUIDRC_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixSPHPrecipitationDissolutionLiquidRC : public Fix {
  public:
    FixSPHPrecipitationDissolutionLiquidRC(class LAMMPS *, int, char **);
    virtual ~FixSPHPrecipitationDissolutionLiquidRC();
    int setmask();
    virtual void init();
    virtual void init_list(int, class NeighList *);
    virtual void end_of_step();

    int pack_forward_comm(int, int *, double *, int, int *);
    void unpack_forward_comm(int, int, double *);

    int pack_reverse_comm(int, int, double *);
    void unpack_reverse_comm(int, int *, double *);

  private:
    class NeighList *list;

  protected:
    int is_periodic;
    double mMthres, cCeq, cAeq, dtcA;
    double *step_respa;
    double *cC, *cA, *dcA, *mM, *local_pot;
    double *ischangecC, *distdcA;
    int  nearLiquidcount;

    class Pair *pair;
  };

}

#endif
#endif
