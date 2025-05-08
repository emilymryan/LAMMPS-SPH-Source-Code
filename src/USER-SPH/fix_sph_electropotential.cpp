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
#include "fix_sph_electropotential.h"
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

FixSPHElectroPotential::FixSPHElectroPotential(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->esph_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
               "fix sph/electropotential command requires atom_style with both energy and density, e.g. meso");

  if (narg != 8)
    error->all(FLERR,"Illegal number of arguments for fix sph/electropotential command");

  // required args
  int m = 3;

  // Get the kernal length
  h = atof(arg[m++]);
  // Get the conversion from concentration to local charge
  conc_to_charge = atof(arg[m++]);
  // Get the applied potential
  applied_pot= atof(arg[m++]);
  // Get the position of the anode end
  electrolyte_start  = atof(arg[m++]);
  // Get the position of the boundary condition
  electrolyte_end  = atof(arg[m++]);
  
  time_integrate = 0;

    // find the anion concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property cA for fix sph/electropotential");
  cA = atom->dvector[icA];

    // find the cation concentration property
  int fcC;
  int icC = atom->find_custom("cC", fcC);
  if (icC < 0)
    error->all(FLERR,
	       "Can't find property cC for fix sph/electropotential");
  cC = atom->dvector[icC];

    // find the local potential
  int flocal_pot;
  int ilocal_pot = atom->find_custom("local_pot", flocal_pot);
  if (ilocal_pot < 0)
    error->all(FLERR,
	       "Can't find property local_pot for fix sph/electropotential");
  local_pot = atom->dvector[ilocal_pot];
}

/* ---------------------------------------------------------------------- */

int FixSPHElectroPotential::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHElectroPotential::init() {

  // need a full neighbor list, built whenever re-neighboring occurs
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void FixSPHElectroPotential::init_list(int, NeighList *ptr) {
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixSPHElectroPotential::initial_integrate(int /*vflag*/)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;

  double delx, dely, delz;
  
  int *ilist, *jlist, *numneigh, **firstneigh;
  double ih;
  double r;

  double chargedensity;
  double pi;
  double cC_value, cA_value, x_value, y_value;
  double **x = atom->x;

  double *rmass = atom->rmass;
  double *rho = atom->rho;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  int *mask = atom->mask;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  pi = 3.1415926535;
  
  for (ii = 0; ii < inum; ii++)
    {
      local_pot[ii] = 0;
    }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    if (itype==1 || itype==3){
      // Calculate charge density at i
      chargedensity = conc_to_charge * (cC[i] - cA[i]);
      local_pot[i] =  chargedensity + applied_pot * ( x[i][1] - electrolyte_start)/(electrolyte_end-electrolyte_start);

      // Check neighbouring atoms
      jnum = numneigh[i];
      jlist = firstneigh[i];
      
      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];
	j &= NEIGHMASK;
	jtype = type[j];

	// Calculate the distance between the particles
	delx = x[i][0] - x[j][0];
	dely = x[i][1] - x[j][1];
	delz = x[i][2] - x[j][2];
	r = sqrt(delx * delx + dely * dely + delz * delz);
	
	// Check if j is within the support kernel and a liquid
	if (r < h && jtype==1) {
	  // Calculating the charge density and local potential for the  time step
	  chargedensity = conc_to_charge * (cC[j] - cA[j]);
	 // Net charge density
	  local_pot[i] -= chargedensity / (4 * pi) / r;
	  
	} // loop inside support kernel
      }// for loop jj
    } else if (itype==2) { // Potential is 0 in anode
      local_pot[i] = 0;
    }
  } // loop through i
}

/* ---------------------------------------------------------------------- */
