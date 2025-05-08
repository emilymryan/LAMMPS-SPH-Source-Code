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
#include "fix_sph_precipitation_dissolution_liquidRC.h"
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

FixSPHPrecipitationDissolutionLiquidRC::FixSPHPrecipitationDissolutionLiquidRC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->esph_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
  	       "fix sph/precipitation/dissolution/liquidRC command requires atom_style with both energy and density, e.g. meso");

  if (narg != 6)
    error->all(FLERR,"Illegal number of arguments for fix sph/precipitation/dissolution/liquidRC command");

  time_integrate = 0;

  // Obtain the neighbor cutoff distance
  int m = 3;
  mMthres = atof(arg[m++]);
  cCeq = atof(arg[m++]);
  cAeq = atof(arg[m++]);
  if (mMthres <= 0)
    error->all(FLERR,"Illegal value for mass threshold");
  if (cCeq <= 0)
    error->all(FLERR,"Illegal value for equilibrium cation concentration");
  if (cAeq <= 0)
    error->all(FLERR,"Illegal value for equilibrium anion concentration");

  // find the concentration property
  int fcC;
  int icC = atom->find_custom("cC", fcC);
  if (icC < 0)
    error->all(FLERR,
	       "Can't find property cC for fix sph/precipitation/dissolution/liquidRC");
  cC = atom->dvector[icC];

  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property cA for fix sph/precipitation/dissolution/liquidRC");
  cA = atom->dvector[icA];

    // find the concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (idcA < 0)
    error->all(FLERR,
	       "Can't find property dcA for fix sph/precipitation/dissolution/liquidRC");
  dcA = atom->dvector[idcA];

  // find the mass of A property
  int fmM;
  int imM = atom->find_custom("mM", fmM);
  if (imM < 0)
    error->all(FLERR,
	       "Can't find property mM for fix sph/precipitation/dissolution/liquidRC");
  mM = atom->dvector[imM];

  // find the potential property
  int flocal_pot;
  int ilocal_pot = atom->find_custom("local_pot", flocal_pot);
  if (ilocal_pot < 0)
    error->all(FLERR,
	       "Can't find property local_pot for fix sph/precipitation/dissolution/liquidRC");
  local_pot = atom->dvector[ilocal_pot];

  // find the potential property                                                                                                                                 
  // int fisAnode;
  // int iisAnode = atom->find_custom("isAnode", fisAnode);
  //if (ilocal_pot < 0)
  //error->all(FLERR,
  //           "Can't find property isAnode for fix sph/precipitation/dissolution/liquidRC");
  //isAnode = atom->dvector[iisAnode];
  
  // set comm size by this fix
  comm_forward = 2;
  comm_reverse = 2;
}

/* ---------------------------------------------------------------------- */

FixSPHPrecipitationDissolutionLiquidRC::~FixSPHPrecipitationDissolutionLiquidRC() {
}

/* ---------------------------------------------------------------------- */

int FixSPHPrecipitationDissolutionLiquidRC::setmask() {
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHPrecipitationDissolutionLiquidRC::init() {
  // need a full neighbor list, built whenever re-neighboring occurs
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  dtcA = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixSPHPrecipitationDissolutionLiquidRC::init_list(int, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixSPHPrecipitationDissolutionLiquidRC::end_of_step()
{
  int i, j, ii, jj, itype, jtype, jshortest, inum, jnum;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double delx, dely, delz;
  double shortest, rsq, r, cAdist;

  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  bool foundShortest;

  // Create memory for ischangecC
  memory->create(ischangecC, nall, "fix:ischangecC");
  // Create memory for distributed cA from precipitation
  memory->create(distdcA, nall, "fix:distdcA");

  // Initilise all phase change variables
  for (i=0; i < nall; i++) {
    ischangecC[i] = 0.0;
    distdcA[i] = 0.0;
  }
  
  for (ii=0; ii<inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    // Only deal with local solid particles
    // and check for periodicity
    if ((i<nlocal) && (itype == 2)) {
      if (mM[i] > mMthres) { // precipitation
	// Check neighbouring atoms
	jnum = numneigh[i];
        jlist = firstneigh[i];

        // Then need to find the closest fluid particles
        shortest = 100000.0;
        foundShortest = false;
	nearLiquidcount = 0;
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          jtype = type[j];
	  delx = x[i][0] - x[j][0];
	  dely = x[i][1] - x[j][1];
	  delz = x[i][2] - x[j][2];
	  r = sqrt(delx*delx + dely*dely + delz*delz);

	  if (jtype==1) {nearLiquidcount++;}
	  // find and store the shortest distance liquid particle that won't be precipitated yet
	  if ((r<shortest) && (jtype==1) && ischangecC[j]==0.0) {
	    // check that we are not triggering precipitation in ghost zone
	    // Also, to avoid one precipitation triggering many
	    foundShortest = true;
	    shortest = r;
	    //printf("shortest distance: %d\n", shortest);
	    jshortest = j;
	  } // find and store the shortest distance
	} // for loop to find closest fluid
        // if there is a closest liquid particle
	if ((foundShortest)) { // In electrolyte
	  // Trigger phase change for i (solid) and j (liquid) 
	  ischangecC[i] += 1.0;
	  ischangecC[jshortest] += 1.0;
	  cAdist = cA[jshortest]/(nearLiquidcount-1); // Calculate the cA that will be distributed to near by liquid particles
	  // Loop thru the near by liquid particles to distribute cA of the new solid
	  for (jj = 0; jj < jnum; jj++) { 
	    j = jlist[jj];
	    j &= NEIGHMASK;
	    jtype = type[j];
	    if (jtype==1 and j!=jshortest) {
	      cA[j]+= cAdist;
	      distdcA[j] += cAdist;
	    }
	  }
	}
      }
      else if ((mM[i] <= -mMthres)) { // Dissolution in electrolyte
	ischangecC[i] += 1.0;
      }
    }
  } // end of all particle loop
  
  // Send the phase change status
  comm->reverse_comm_fix(this);

  // Then loop through the local atoms to see if there is a trigger
  for (i = 0; i < nlocal; i++) {
    if (ischangecC[i] > 0.01) {
      // Phase change happening
      if (type[i] == 1) {
        // Changing liquid into solid
        mM[i] = 0.0;
        type[i] = 2; // convert the liquid to the solid
        cC[i] = 0.0; // concentration is 0.0
	cA[i] = 0.0; // concentration is 0.0
	local_pot[i] = 0.0;
        v[i][0] = 0.0; // set velocity to 0.0
        v[i][1] = 0.0;
        v[i][2] = 0.0;
	//isAnode[i] = 0.0;
      }
      else if (type[i] == 2) { // In electrolyte
        // changing solid into liquid
        if (mM[i] <= 0.0) { //0.0) {
          mM[i] = 0.0; // mass becomes 0.0 for dissolution
          type[i] = 1; // convert solid to liquid
          cC[i] = cCeq; // concentration set to equilibrium
	  cA[i] = cAeq; // concentration set to equilibrium
	  local_pot[i] = 0.0;
	  //isAnode[i] = 0.0;
        }
	else if (mM[i] > mMthres) { // In electrolyte
          // Only causing the particle to reduce the mass because of precipitation
          mM[i] = mM[i] - mMthres;
        }
      }
    }
  }

    // Send the new cC and cA to the ghost atoms
   comm->forward_comm_fix(this);
  // Destroy the memory created for ischangecC
  memory->destroy(ischangecC);
  memory->destroy(distdcA);
}
/* ---------------------------------------------------------------------- */

int FixSPHPrecipitationDissolutionLiquidRC::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = cC[j];
    buf[m++] = cA[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixSPHPrecipitationDissolutionLiquidRC::unpack_forward_comm(int n, int first, double *buf) 
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    cC[i] += buf[m++];
    cA[i] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int FixSPHPrecipitationDissolutionLiquidRC::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  int *type = atom->type;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = ischangecC[i];
    buf[m++] = distdcA[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixSPHPrecipitationDissolutionLiquidRC::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  int *type = atom->type;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    ischangecC[j] += buf[m++];
    cA[j] += buf[m++];
  }
}
