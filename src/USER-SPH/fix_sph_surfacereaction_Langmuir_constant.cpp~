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
#include "fix_sph_surfacereaction_Langmuir_constant.h"
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

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSPHSurfaceReactionLangmuirConstant::FixSPHSurfaceReactionLangmuirConstant(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
	       "fix sph/surfacereaction/Langmuir/constant command requires atom_style with both energy and density, e.g. meso");

  if (narg != 4)
    error->all(FLERR,"Illegal number of arguments for fix sph/surfacereaction/Langmuir/constant command");

  time_integrate = 0;

  // Required args
  int m = 3;
  constxA = atof(arg[m++]);

  // find the concentration property
  int fxA;
  int ixA = atom->find_custom("xA", fxA);
  if (ixA < 0)
    error->all(FLERR,
	       "Can't find property xA for fix sph/surfacereaction/Langmuir/constant");
  xA = atom->dvector[ixA];

  // find the local concentration property
  int fdxA;
  int idxA = atom->find_custom("dxA", fdxA);
  if (idxA < 0)
    error->all(FLERR,
	       "Can't find property dxA for fix sph/surfacereaction/Langmuir/constant");
  dxA = atom->dvector[idxA];
}

/* ---------------------------------------------------------------------- */

int FixSPHSurfaceReactionLangmuirConstant::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirConstant::initial_integrate(int vflag) {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xA[i] = constxA; // keep the concentration at the constant level
      dxA[i] = 0.0; // and the change is 0.0
    }
  }
}

/* ---------------------------------------------------------------------- */
