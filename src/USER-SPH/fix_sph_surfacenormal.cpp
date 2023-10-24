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
#include "fix_sph_surfacenormal.h"
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

FixSPHSurfaceNormal::FixSPHSurfaceNormal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->esph_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
               "fix sph/surfacenormal command requires atom_style with both energy and density, e.g. meso");

  if (narg != 4)
    error->all(FLERR,"Illegal number of arguments for fix sph/surfacenormal command");

  // required args
  int m = 3;
  d = atof(arg[m++]);

  time_integrate = 0;

  // find the normal vector
  int fnx;
  int inx = atom->find_custom("nx", fnx);
  if (inx < 0)
    error->all(FLERR,
	       "Can't find property nx for fix sph/surfacenormal");
  nx = atom->dvector[inx];

  int fny;
  int iny = atom->find_custom("ny", fny);
  if (iny < 0)
    error->all(FLERR,
	       "Can't find property ny for fix sph/surfacenormal");
  ny = atom->dvector[iny];

  int fnz;
  int inz = atom->find_custom("nz", fnz);
  if (inz < 0)
    error->all(FLERR,
	       "Can't find property nz for fix sph/surfacenormal");
  nz = atom->dvector[inz];
}

/* ---------------------------------------------------------------------- */

int FixSPHSurfaceNormal::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceNormal::init() {

  // need a full neighbor list, built whenever re-neighboring occurs
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceNormal::init_list(int, NeighList *ptr) {
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceNormal::initial_integrate(int /*vflag*/)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;

  double delx, dely, delz, rdelx, rdely, rdelz;
  double nmag;
  
  int *ilist, *jlist, *numneigh, **firstneigh;
  double ih;
  double r, wf, wfd;
  double wfdx, nwfdx, wfdy, nwfdy, wfdz, nwfdz;
  
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

  for (ii = 0; ii < inum; ii++)
    {
      nx[ii] = 0;
      ny[ii] = 0;
      nz[ii] = 0;
    }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

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
      rdelx = fabs(delx);
      rdely = fabs(dely);
      rdelz = fabs(delz);

      r = sqrt(delx * delx + dely * dely + delz * delz);
	
      // Check if j is within the support kernel
      if (r < d) {
	ih = 1.0/d;

      	// kernel function
	if (domain->dimension == 3) {
	  // printf("3D sims not vetted");
	  wfd = sph_dw_quintic3d(r*ih);
	  wfd = wfd*ih*ih*ih*ih;
	  wfdx = sph_dw_quintic3d(rdelx*ih);
	  wfdx = wfdx*ih*ih*ih*ih;
	  nwfdx=wfdx;
	  if(delx<0){nwfdx=wfdx*-1;}
	  wfdy = sph_dw_quintic3d(rdely*ih);
	  wfdy = wfdy*ih*ih*ih*ih;
	  nwfdy = wfdy;
	  if(dely<0){nwfdy=wfdy*-1;}
	  wfdz = sph_dw_quintic3d(rdelz*ih);
	  wfdz = wfdz*ih*ih*ih*ih;
	  nwfdz = wfdz;
	  if(delz<0){nwfdz=wfdz*-1;}
	  wf = sph_kernel_quintic3d(r*ih)*ih*ih*ih;
	} else {
	  wfd = sph_dw_quintic2d(r*ih);
	  wfd = wfd*ih*ih*ih;
	  wfdx = sph_dw_quintic2d(rdelx*ih);
	  wfdx = wfdx*ih*ih*ih;
	  nwfdx = wfdx;
	  if(delx<0){nwfdx=wfdx*-1;}
	  wfdy = sph_dw_quintic2d(rdely*ih);
	  wfdy = wfdy*ih*ih*ih;
	  nwfdy = wfdy;
	  if(dely<0){nwfdy=wfdy*-1;}
	  wf = sph_kernel_quintic2d(r*ih)*ih*ih;
	}
	if ((itype==2) && (jtype==1)||(itype==1) && (jtype==2)) { // solid-fluid interaction only
	  nx[i] += nwfdx;
	  ny[i] += nwfdy;
	} 
      } // loop inside support kernel
    } // for loop jj
    nmag = sqrt(nx[i]*nx[i] + ny[i]*ny[i]);// + nz[i]*nz[i]);
    if(nmag==0) {
      nx[i] = 0;
      ny[i] = 0;
    }
    else {
      nx[i] = nx[i]/nmag;
      ny[i] = ny[i]/nmag;
    }
  } // loop through i
}

/* ---------------------------------------------------------------------- */
