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

FixStyle(sph/electropotential,FixSPHElectroPotential)

#else

#ifndef LMP_FIX_SPH_ELECTROPOTENTIAL_H
#define LMP_FIX_SPH_ELECRTOPOTENTIAL_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixSPHElectroPotential : public Fix {
  public:
    FixSPHElectroPotential(class LAMMPS *, int, char **);
    int setmask();
    virtual void init();
    virtual void init_list(int, class NeighList *);
    virtual void initial_integrate(int);

  private:
    class NeighList *list;

  protected:
    // Kernal size
    double h;
    double conc_to_charge, applied_pot, electrolyte_start, electrolyte_end;
    double *cC, *cA, *local_pot;
    class Pair *pair;
  };

}

#endif
#endif
