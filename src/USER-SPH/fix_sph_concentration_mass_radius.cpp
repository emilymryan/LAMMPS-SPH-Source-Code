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
#include "fix_sph_concentration_mass_radius.h"
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

FixSPHConcentrationMassRadius::FixSPHConcentrationMassRadius(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {
  if ((atom->esph_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
	       "fix sph/concentration/mass/radius command requires atom_style with both energy and density, e.g. meso");

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix sph/concentration/mass/radius command");

  time_integrate = 0;

  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property cA for fix sph/concentration/mass/radius");
  cA = atom->dvector[icA];

  // find the concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property dcA for fix sph/concentration/mass/radius");
  dcA = atom->dvector[idcA];

  // find the concentration property
  int fcC;
  int icC = atom->find_custom("cC", fcC);
  if (icC < 0)
    error->all(FLERR,
	       "Can't find property cC for fix sph/concentration/mass/radius");
  cC = atom->dvector[icC];

  // find the concentration property
  int fdcC;
  int idcC = atom->find_custom("dcC", fdcC);
  if (icC < 0)
    error->all(FLERR,
	       "Can't find property dcC for fix sph/concentration/mass/radius");
  dcC = atom->dvector[idcC];

  // find the radius  property                                                                                                                                                  
  int frR;
  int irR = atom->find_custom("rR", frR);
  if (irR < 0)
    error->all(FLERR,
               "Can't find property rR for fix sph/concentration/mass/radius");
  rR = atom->dvector[irR];

  // find the solid-liquid interaction                                                                                                                                          
  int fdrR;
  int idrR = atom->find_custom("drR", fdrR);
  if (idrR < 0)
    error->all(FLERR,
               "Can't find property drR for fix sph/concentration/mass/radius");
  drR = atom->dvector[idrR];

  //find the radius  property                                                                                                         
  int fnN;
  int inN = atom->find_custom("nN", fnN);
  if (inN < 0)
    error->all(FLERR,
           "Can't find property nN for fix sph/concentration/mass/radius");
  nN = atom->dvector[inN];
  
  // find the solid-liquid interaction
  int fdnN;
  int idnN = atom->find_custom("dnN", fdnN);
  if (idnN < 0)
    error->all(FLERR,
           "Can't find property dnN for fix sph/concentration/mass/radius");
  dnN = atom->dvector[idnN];

  // find the mass of A property
  int fmM;
  int imM = atom->find_custom("mM", fmM);
  if (imM < 0)
    error->all(FLERR,
	       "Can't find property mM for fix sph/concentration/mass/radius");
  mM = atom->dvector[imM];

  // find the solid-liquid interaction
  int fdmM;
  int idmM = atom->find_custom("dmM", fdmM);
  if (idmM < 0)
    error->all(FLERR,
	       "Can't find property dmM for fix sph/concentration/mass/radius");
  dmM = atom->dvector[idmM];

}

/* ---------------------------------------------------------------------- */

FixSPHConcentrationMassRadius::~FixSPHConcentrationMassRadius() {
}

/* ---------------------------------------------------------------------- */

int FixSPHConcentrationMassRadius::setmask() {
  int mask = 0;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHConcentrationMassRadius::init() {
  dtcA = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixSPHConcentrationMassRadius::final_integrate() {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  double **x = atom->x;
  double xtmp, ytmp;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      cA[i] += dtcA*dcA[i];
      cC[i] += dtcA*dcC[i];
      nN[i] += dtcA*dnN[i];
      //printf("cC: %.10f \n",cC);
      // Only update mass, radius, and number of nuclei for solid particles
      if (type[i] == 2){
	mM[i] += dtcA*dmM[i];
	rR[i] += dtcA*drR[i];
	//nN[i] += dtcA*dnN[i];
	//printf("nN: %.10f \n",nN);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPHConcentrationMassRadius::reset_dt() {
  dtcA = update->dt;
}
