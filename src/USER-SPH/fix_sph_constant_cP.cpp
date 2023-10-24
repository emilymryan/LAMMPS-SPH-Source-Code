/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include "fix_sph_constant_cP.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "pair.h"
#include <cstdio>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSPHConstantcP::FixSPHConstantcP(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->esph_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
  	       "fix sph/constant/cP command requires atom_style with both energy and density, e.g. meso");

  if (narg != 4)
    error->all(FLERR,"Illegal number of arguments for fix sph/constant/cP command");

  time_integrate = 0;

  // Required args
  int m = 3;
  constcP = atof(arg[m++]);

  // find the concentration property
  int fcP;
  int icP = atom->find_custom("cP", fcP);
  if (icP < 0)
    error->all(FLERR,
	       "Can't find property cP for fix sph/constant/cP");
  cP = atom->dvector[icP];

  // find the local concentration property
  int fdcP;
  int idcP = atom->find_custom("dcP", fdcP);
  if (idcP < 0)
    error->all(FLERR,
	       "Can't find property dcP for fix sph/constant/cP");
  dcP = atom->dvector[idcP];
}

/* ---------------------------------------------------------------------- */

int FixSPHConstantcP::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHConstantcP::initial_integrate(int vflag) {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;

  double **x = atom->x;
  
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      //      printf("y: %.6f, cA: %3.3f\n", x[i], constcA);
      cP[i] = constcP; // keep the concentration at the constant level
      dcP[i] = 0.0; // and the change is 0.0
    }
  }
}

/* ---------------------------------------------------------------------- */
