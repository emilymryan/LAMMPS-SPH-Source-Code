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
#include "fix_sph_surfacereaction_Langmuir_porous_simple.h"
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
#include <float.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSPHSurfaceReactionLangmuirPorousSimple::FixSPHSurfaceReactionLangmuirPorousSimple(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
               "fix sph/surfacereaction/Langmuir/porous/simple command requires atom_style with both energy and density, e.g. meso");

  if (narg != 4)
    error->all(FLERR,"Illegal number of arguments for fix sph/surfacereaction/Langmuir/porous/simple command");

  // required args
  int m = 3;
  yAmax = atof(arg[m++]);

  time_integrate = 0;

  // find the aqueous mass fraction property
  int fxA;
  int ixA = atom->find_custom("xA", fxA);
  if (ixA < 0)
    error->all(FLERR,
               "Can't find property xA for fix sph/surfacereaction/Langmuir/porous/simple");
  xA = atom->dvector[ixA];

  // find the change in aqueous mass fraction concentration property
  int fdxA;
  int idxA = atom->find_custom("dxA", fdxA);
  if (idxA < 0)
    error->all(FLERR,
               "Can't find property dxA for fix sph/surfacereaction/Langmuir/porous/simple");
  dxA = atom->dvector[idxA];

  // find the absorbed mass fraction property
  int fyA;
  int iyA = atom->find_custom("yA", fyA);
  if (iyA < 0)
    error->all(FLERR,
               "Can't find property yA for fix sph/surfacereaction/Langmuir/porous/simple");
  yA = atom->dvector[iyA];

  // find the change in the absorbed mass fraction
  int fdyA;
  int idyA = atom->find_custom("dyA", fdyA);
  if (idyA < 0)
    error->all(FLERR,
               "Can't find property dyA for fix sph/surfacereaction/Langmuir/porous/simple");
  dyA = atom->dvector[idyA];

  // find the normalised concentration of adsorbed species
  int fthetaA;
  int ithetaA = atom->find_custom("thetaA", fthetaA);
  if (ithetaA < 0)
    error->all(FLERR,
               "Can't find property thetaA for fix sph/surfacereaction/Langmuir/porous/simple");
  thetaA = atom->dvector[ithetaA];
}

/* ---------------------------------------------------------------------- */

int FixSPHSurfaceReactionLangmuirPorousSimple::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirPorousSimple::init() {
  dtxA = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirPorousSimple::initial_integrate(int /*vflag*/)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (type[i] == 2) // Only update mass fraction for solid particels
        thetaA[i] = yA[i]/yAmax;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirPorousSimple::final_integrate() {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (type[i] <= 2) {
        double ddxA = (isnan(dxA[i]) || (std::abs(dxA[i]) < DBL_EPSILON)) ? 0.0 : dxA[i];
        xA[i] = std::max(0.0, xA[i] + dtxA*ddxA);
        if (type[i] == 2) { // Only update mass fraction for solid particels
          double ddyA = (isnan(dyA[i]) || (std::abs(dyA[i]) < DBL_EPSILON)) ? 0.0 : dyA[i];
          yA[i] = std::max(0.0, yA[i] + dtxA*ddyA);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirPorousSimple::reset_dt() {
  dtxA = update->dt;
}

/* ---------------------------------------------------------------------- */
