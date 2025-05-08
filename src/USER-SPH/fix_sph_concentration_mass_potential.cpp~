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
#include "fix_sph_migration_surfacereaction_simple.h"
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

FixSPHMigrationSurfaceReactionSimple::FixSPHMigrationSurfaceReactionSimple(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {
  if ((atom->esph_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
	       "fix sph/migration/surfacereaction/simple command requires atom_style with both energy and density, e.g. meso");

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix sph/migration/surfacereaction/simple command");

  time_integrate = 0;

  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property cA for fix sph/migration/surfacereaction/simple");
  cA = atom->dvector[icA];

  // find the concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property dcA for fix sph/migration/surfacereaction/simple");
  dcA = atom->dvector[idcA];

  // find the concentration property
  int fcC;
  int icC = atom->find_custom("cC", fcC);
  if (icC < 0)
    error->all(FLERR,
	       "Can't find property cC for fix sph/migration/surfacereaction/simple");
  cC = atom->dvector[icC];

  // find the concentration property
  int fdcC;
  int idcC = atom->find_custom("dcC", fdcC);
  if (icC < 0)
    error->all(FLERR,
	       "Can't find property dcC for fix sph/migration/surfacereaction/simple");
  dcC = atom->dvector[idcC];

  // find the mass of A property
  int fmM;
  int imM = atom->find_custom("mM", fmM);
  if (imM < 0)
    error->all(FLERR,
	       "Can't find property mM for fix sph/migration/surfacereaction/simple");
  mM = atom->dvector[imM];

  // find the solid-liquid interaction
  int fdmM;
  int idmM = atom->find_custom("dmM", fdmM);
  if (idmM < 0)
    error->all(FLERR,
	       "Can't find property dmM for fix sph/migration/surfacereaction/simple");
  dmM = atom->dvector[idmM];

  // find the local potential
  int flocal_pot;
  int ilocal_pot = atom->find_custom("local_pot", flocal_pot);
  if (ilocal_pot < 0)
    error->all(FLERR,
               "Can't find property local_pot for pair_style sph/surfacereaction/simple/charge");
  local_pot = atom->dvector[ilocal_pot];

  // find the local next potential
  int fnext_local_pot;
  int inext_local_pot = atom->find_custom("next_local_pot", fnext_local_pot);
  if (inext_local_pot < 0)
    error->all(FLERR,
               "Can't find property next_local_pot for pair_style sph/surfacereaction/simple/charge");
  next_local_pot = atom->dvector[inext_local_pot];
   // find the x-component normal
  int fnx;
  int inx = atom->find_custom("nx", fnx);
  if (inx < 0)
    error->all(FLERR,
	       "Can't find property nx for pair_style sph/migration/surfacereaction/simple/charge");
  nx = atom->dvector[inx];
  
  // find the next x-component normal
  int fdnx;
  int idnx = atom->find_custom("dnx", fdnx);
  if (idnx < 0)
    error->all(FLERR,
	       "Can't find property dnx for pair_style sph/migration/surfacereaction/simple/charge");
  dnx = atom->dvector[idnx];
  // find the y-component normal
  int fny;
  int iny = atom->find_custom("ny", fny);
  if (iny < 0)
    error->all(FLERR,
	       "Can't find property ny for pair_style sph/migration/surfacereaction/simple/charge");
  ny = atom->dvector[iny];
  
  // find the next y-component normal
  int fdny;
  int idny = atom->find_custom("dny", fdny);
  if (idny < 0)
    error->all(FLERR,
	       "Can't find property dny for pair_style sph/migration/surfacereaction/simple/charge");
  dny = atom->dvector[idny];
  
  // find the z-component normal
  int fnz;
  int inz = atom->find_custom("nz", fnz);
  if (inz < 0)
    error->all(FLERR,
	       "Can't find property nz for pair_style sph/migration/surfacereaction/simple/charge");
  nz = atom->dvector[inz];
  
  // find the next z-component normal
  int fdnz;
  int idnz = atom->find_custom("dnz", fdnz);
  if (idnz < 0)
    error->all(FLERR,
	       "Can't find property dnz for pair_style sph/migration/surfacereaction/simple/charge");
  dnz = atom->dvector[idnz];

}

/* ---------------------------------------------------------------------- */

FixSPHMigrationSurfaceReactionSimple::~FixSPHMigrationSurfaceReactionSimple() {
}

/* ---------------------------------------------------------------------- */

int FixSPHMigrationSurfaceReactionSimple::setmask() {
  int mask = 0;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHMigrationSurfaceReactionSimple::init() {
  dtcA = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixSPHMigrationSurfaceReactionSimple::final_integrate() {
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
      local_pot[i] = next_local_pot[i];
      //      dcA[i] = 0;
      //      dcC[i] = 0;
      nx[i] = dnx[i];
      ny[i] = dny[i];
      nz[i] = dnz[i];
      // dnx[i] = 0.0;
      // dny[i] = 0.0;
      // dnz[i] = 0.0;
      // Only update mass for solid particles and set potential to 0
      if (type[i] == 2){
	mM[i] += dtcA*dmM[i];
	//	dmM[i] = 0;
	local_pot[i] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPHMigrationSurfaceReactionSimple::reset_dt() {
  dtcA = update->dt;
}
