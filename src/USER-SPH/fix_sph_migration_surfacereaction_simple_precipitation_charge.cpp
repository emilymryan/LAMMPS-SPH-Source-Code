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
#include "fix_sph_migration_surfacereaction_simple_precipitation_charge.h"
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

FixSPHMigrationSurfaceReactionSimplePrecipitationCharge::FixSPHMigrationSurfaceReactionSimplePrecipitationCharge(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->esph_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
  	       "fix sph/migration/surfacereaction/simple/precipitation/charge command requires atom_style with both energy and density, e.g. meso");

  if (narg != 7)
    error->all(FLERR,"Illegal number of arguments for fix sph/migration/surfacereaction/simple/precipitation/charge command");

  time_integrate = 0;

  // Obtain the neighbor cutoff distance
  int m = 3;
  mMthres = atof(arg[m++]);
  cCeq = atof(arg[m++]);
  cAeq = atof(arg[m++]);
  is_periodic = atoi(arg[m++]);
  if (mMthres <= 0)
    error->all(FLERR,"Illegal value for mass threshold");
  if (cCeq <= 0)
    error->all(FLERR,"Illegal value for equilibrium cation concentration");
  if (cAeq <= 0)
    error->all(FLERR,"Illegal value for equilibrium anion concentration");
  if ((is_periodic!=0) && (is_periodic!=1))
    error->all(FLERR,"Illegal value for setting periodicity of reaction");

  // find the concentration property
  int fcC;
  int icC = atom->find_custom("cC", fcC);
  if (icC < 0)
    error->all(FLERR,
	       "Can't find property cC for fix sph/migration/surfacereaction/simple/precipitation/charge");
  cC = atom->dvector[icC];

  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property cA for fix sph/migration/surfacereaction/simple/precipitation/charge");
  cA = atom->dvector[icA];

  // find the mass of A property
  int fmM;
  int imM = atom->find_custom("mM", fmM);
  if (imM < 0)
    error->all(FLERR,
	       "Can't find property mM for fix sph/migration/surfacereaction/simple/precipitation/charge");
  mM = atom->dvector[imM];

  // find the potential property
  int flocal_pot;
  int ilocal_pot = atom->find_custom("local_pot", flocal_pot);
  if (ilocal_pot < 0)
    error->all(FLERR,
	       "Can't find property local_pot for fix sph/migration/surfacereaction/simple/precipitation/charge");
  local_pot = atom->dvector[ilocal_pot];

  // find the reaction rate coefficient for the cations                                     
  int fRC;
  int iRC = atom->find_custom("RC", fRC);
  if (iRC < 0)
    error->all(FLERR,
               "Can't find property RC for pair_style sph/migration/surfacereaction/simple/\
   charge");
  RC = atom->dvector[iRC];

  // set comm size by this fix
  comm_forward = 3;
  comm_reverse = 2;
}

/* ---------------------------------------------------------------------- */

FixSPHMigrationSurfaceReactionSimplePrecipitationCharge::~FixSPHMigrationSurfaceReactionSimplePrecipitationCharge() {
}

/* ---------------------------------------------------------------------- */

int FixSPHMigrationSurfaceReactionSimplePrecipitationCharge::setmask() {
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHMigrationSurfaceReactionSimplePrecipitationCharge::init() {
  // need a full neighbor list, built whenever re-neighboring occurs
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void FixSPHMigrationSurfaceReactionSimplePrecipitationCharge::init_list(int, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixSPHMigrationSurfaceReactionSimplePrecipitationCharge::end_of_step()
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
  // Create memory for newRC
  memory->create(newRC, nall, "fix:newRC");

  // Initilise all phase change variables
  for (i=0; i < nall; i++) {
    ischangecC[i] = 0.0;
    newRC[i] = 0.0;
  }

  for (ii=0; ii<inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    // Only deal with local solid particles
    // and check for periodicity
    if (((i<nlocal) && (itype == 2)) &&
        ((is_periodic==1) || (not (x[i][0] < domain->boxlo[0] || x[i][0] > domain->boxhi[0] ||
                                   x[i][1] < domain->boxlo[1] || x[i][1] > domain->boxhi[1] ||
                                   x[i][2] < domain->boxlo[2] || x[i][2] > domain->boxhi[2]))) ) {
      if (mM[i] > mMthres) { // precipitation
        // Check neighbouring atoms
        jnum = numneigh[i];
        jlist = firstneigh[i];

        // Then need to find the closest fluid particles
        shortest = 100000.0;
        foundShortest = false;

        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          jtype = type[j];

          // Check for periodicity
          if ((is_periodic==1) || (not (x[j][0] < domain->boxlo[0] || x[j][0] > domain->boxhi[0] ||
                                        x[j][1] < domain->boxlo[1] || x[j][1] > domain->boxhi[1] ||
                                        x[j][2] < domain->boxlo[2] || x[j][2] > domain->boxhi[2]))) {
            delx = x[i][0] - x[j][0];
            dely = x[i][1] - x[j][1];
            delz = x[i][2] - x[j][2];
            r = sqrt(delx*delx + dely*dely + delz*delz);

            // find and store the shortest distance liquid particle that won't be precipitated yet
            if ((r<shortest) && (jtype==1) && ischangecC[j]==0.0) {
              // check that we are not triggering precipitation in ghost zone
              // Also, to avoid one precipitation triggering many
              foundShortest = true;
              shortest = r;
              jshortest = j;
            } // find and store the shortest distance
          } // check for periodicity
        } // for loop to find closest fluid
        // if there is a closest liquid particle
        if ((foundShortest)) { // In electrolyte
          // Trigger phase change for i (solid) and j (liquid) and store RC for new solid
          ischangecC[i] += 1.0;
          ischangecC[jshortest] += 1.0;
	  newRC[jshortest] = RC[i];
	  cAdist = cA[jshortest]/(jnum - 1); // switch to count of local liquid particles, not all jnum
	  for (jj = 0; jj < jnum; jj++) { 

	    j = jlist[jj];
	    j &= NEIGHMASK;
	    jtype = type[j];
	    
	    // Check for periodicity
	    if ((is_periodic==1) || (not (x[j][0] < domain->boxlo[0] || x[j][0] > domain->boxhi[0] ||
                                        x[j][1] < domain->boxlo[1] || x[j][1] > domain->boxhi[1] ||
                                        x[j][2] < domain->boxlo[2] || x[j][2] > domain->boxhi[2]))) {
	      if (jtype==1) {cA[j] += cAdist;}
	    }
	  }
	  
        }
      }
      else if ((mM[i] <= -mMthres)) { // Dissolution in electrolyte
        ischangecC[i] += 1.0;
      }
    }
  } // end of all particle loop
  
  // Send the phase change status and parent RC in reverse first
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
	RC[i] = newRC[i]; // The new solid takes the reaction rate coefficient from the closest solid
	local_pot[i] = 0.0;
        v[i][0] = 0.0; // set velocity to 0.0
        v[i][1] = 0.0;
        v[i][2] = 0.0;
      }
      else if (type[i] == 2) { // In electrolyte
        // changing solid into liquid
        if (mM[i] < -mMthres) {
          mM[i] = 0.0; // mass becomes 0.0 for dissolution
          type[i] = 1; // convert solid to liquid
          cC[i] = cCeq; // concentration set to equilibrium
	  cA[i] = cAeq; // concentration set to equilibrium
	  RC[i] = 0; // Dissolved solids (now liquid) have 0 RC
	  local_pot[i] = 0.0;
        }
        else if (mM[i] > mMthres) { // In electrolyte
          // Only causing the particle to reduce the mass because of precipitation
          mM[i] = mM[i] - mMthres;
        }
      }
    }
  }

    // Send the new RC to the ghost atoms
  comm->forward_comm_fix(this);
  // Destroy the memory created for ischangecC and newRC
  memory->destroy(ischangecC);
  memory->destroy(newRC);
}
/* ---------------------------------------------------------------------- */

int FixSPHMigrationSurfaceReactionSimplePrecipitationCharge::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = RC[j];
    buf[m++] = cC[j];
    buf[m++] = cA[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixSPHMigrationSurfaceReactionSimplePrecipitationCharge::unpack_forward_comm(int n, int first, double *buf) 
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    RC[i] = buf[m++];
    cC[i] += buf[m++];
    cA[i] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int FixSPHMigrationSurfaceReactionSimplePrecipitationCharge::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  int *type = atom->type;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = ischangecC[i];
    if(newRC[i] != 0.0){
      buf[m++] = newRC[i];
    }
    else {
      buf[m++] += newRC[i];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixSPHMigrationSurfaceReactionSimplePrecipitationCharge::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  int *type = atom->type;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    ischangecC[j] += buf[m++];
    if(buf[m] != 0.0){
      newRC[j] = buf[m++];
    }
    else {
      newRC[j] += buf[m++];
    }
  }
}
