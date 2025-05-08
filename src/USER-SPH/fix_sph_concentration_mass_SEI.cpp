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
#include "fix_sph_concentration_mass_SEI.h"
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

FixSPHConcentrationMassSEI::FixSPHConcentrationMassSEI(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {
  if ((atom->esph_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
	       "fix sph/concentration/mass command requires atom_style with both energy and density, e.g. meso");

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix sph/concentration/mass/SEI command");

  time_integrate = 0;

  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property cA for fix sph/concentration/mass/SEI");
  cA = atom->dvector[icA];

  // find the concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property dcA for fix sph/concentration/mass/SEI");
  dcA = atom->dvector[idcA];

  // find the concentration property
  int fcC;
  int icC = atom->find_custom("cC", fcC);
  if (icC < 0)
    error->all(FLERR,
	       "Can't find property cC for fix sph/concentration/mass/SEI");
  cC = atom->dvector[icC];

  // find the concentration property
  int fdcC;
  int idcC = atom->find_custom("dcC", fdcC);
  if (icC < 0)
    error->all(FLERR,
	       "Can't find property dcC for fix sph/concentration/mass/SEI");
  dcC = atom->dvector[idcC];

  // find the mass of A property
  int fmM;
  int imM = atom->find_custom("mM", fmM);
  if (imM < 0)
    error->all(FLERR,
	       "Can't find property mM for fix sph/concentration/mass/SEI");
  mM = atom->dvector[imM];

  // find the solid-liquid interaction
  int fdmM;
  int idmM = atom->find_custom("dmM", fdmM);
  if (idmM < 0)
    error->all(FLERR,
	       "Can't find property dmM for fix sph/concentration/mass/SEI");
  dmM = atom->dvector[idmM];

  // find the SEI thickness
  int fhsei;
  int ihsei = atom->find_custom("hsei", fhsei);
  if (ihsei < 0)
    error->all(FLERR,
	       "Can't find property hsei for fix sph/concentration/mass/SEI");
  hsei = atom->dvector[ihsei];

  // find the SEI change in thickness
  int fdhsei;
  int idhsei = atom->find_custom("dhsei", fdhsei);
  if (idhsei < 0)
    error->all(FLERR,
	       "Can't find property dhsei for fix sph/concentration/mass/SEI");
  dhsei = atom->dvector[idhsei];
  
}

/* ---------------------------------------------------------------------- */

FixSPHConcentrationMassSEI::~FixSPHConcentrationMassSEI() {
}

/* ---------------------------------------------------------------------- */

int FixSPHConcentrationMassSEI::setmask() {
  int mask = 0;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHConcentrationMassSEI::init() {
  dtcA = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixSPHConcentrationMassSEI::final_integrate() {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  double **x = atom->x;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;
 
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      cA[i] += dtcA*dcA[i];
      cC[i] += dtcA*dcC[i];
      hsei[i] += dtcA*dhsei[i];
      // Only update mass for solid particles
      if (type[i] == 2){
	mM[i] += dtcA*dmM[i];
	hsei[i] += dtcA*dhsei[i];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPHConcentrationMassSEI::reset_dt() {
  dtcA = update->dt;
}
