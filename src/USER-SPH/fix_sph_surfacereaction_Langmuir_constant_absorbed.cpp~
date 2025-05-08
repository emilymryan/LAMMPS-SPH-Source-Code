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
#include "fix_sph_surfacereaction_Langmuir_constant_absorbed.h"
#include "sph_kernel_quintic.h"
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
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSPHSurfaceReactionLangmuirConstantAbsorbed::FixSPHSurfaceReactionLangmuirConstantAbsorbed(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
               "fix sph/surfacereaction/Langmuir/constant/absorbed command requires atom_style with both energy and density, e.g. meso");

  if (narg != 4)
    error->all(FLERR,"Illegal number of arguments for fix sph/surfacereaction/Langmuir/constant/absorbed command");

  // required args
  int m = 3;
  thetaAconst = atof(arg[m++]);

  time_integrate = 0;

  // find the aqueous mass fraction property
  int fxA;
  int ixA = atom->find_custom("xA", fxA);
  if (ixA < 0)
    error->all(FLERR,
               "Can't find property xA for fix sph/surfacereaction/Langmuir/constant/absorbed");
  xA = atom->dvector[ixA];

  // find the change in aqueous mass fraction concentration property
  int fdxA;
  int idxA = atom->find_custom("dxA", fdxA);
  if (idxA < 0)
    error->all(FLERR,
               "Can't find property dxA for fix sph/surfacereaction/Langmuir/constant/absorbed");
  dxA = atom->dvector[idxA];

  // find the absorbed mass fraction property
  int fyA;
  int iyA = atom->find_custom("yA", fyA);
  if (iyA < 0)
    error->all(FLERR,
               "Can't find property yA for fix sph/surfacereaction/Langmuir/constant/absorbed");
  yA = atom->dvector[iyA];

  // find the change in the absorbed mass fraction
  int fdyA;
  int idyA = atom->find_custom("dyA", fdyA);
  if (idyA < 0)
    error->all(FLERR,
               "Can't find property dyA for fix sph/surfacereaction/Langmuir/constant/absorbed");
  dyA = atom->dvector[idyA];

  // find the normalised concentration of adsorbed species
  int fthetaA;
  int ithetaA = atom->find_custom("thetaA", fthetaA);
  if (ithetaA < 0)
    error->all(FLERR,
               "Can't find property thetaA for fix sph/surfacereaction/Langmuir/constant/absorbed");
  thetaA = atom->dvector[ithetaA];
}

/* ---------------------------------------------------------------------- */

int FixSPHSurfaceReactionLangmuirConstantAbsorbed::setmask() {
  int mask = 0;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirConstantAbsorbed::init() {
  dtxA = update->dt;

  // need a full neighbor list, built whenever re-neighboring occurs
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirConstantAbsorbed::init_list(int, NeighList *ptr) {
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirConstantAbsorbed::initial_integrate(int /*vflag*/)
{
  int i, ii, inum, itype;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int *type = atom->type;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    if (itype == 2)
      thetaA[i] = thetaAconst;
  }
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirConstantAbsorbed::final_integrate() {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xA[i] += dtxA*dxA[i];
      if (type[i] == 2) // Only update mass fraction for solid particels
        yA[i] += dtxA*dyA[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirConstantAbsorbed::reset_dt() {
  dtxA = update->dt;
}

/* ---------------------------------------------------------------------- */
