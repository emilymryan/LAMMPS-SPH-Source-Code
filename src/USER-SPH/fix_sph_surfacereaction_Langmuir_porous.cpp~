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
#include "fix_sph_surfacereaction_Langmuir_porous.h"
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

FixSPHSurfaceReactionLangmuirPorous::FixSPHSurfaceReactionLangmuirPorous(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
               "fix sph/surfacereaction/Langmuir/porous command requires atom_style with both energy and density, e.g. meso");

  if (narg != 5)
    error->all(FLERR,"Illegal number of arguments for fix sph/surfacereaction/Langmuir/porous command");

  // required args
  int m = 3;
  h = atof(arg[m++]);
  sAmax = atof(arg[m++]);

  time_integrate = 0;

  // find the normal vector
  int fnx;
  int inx = atom->find_custom("nx", fnx);
  if (inx < 0)
    error->all(FLERR,
	       "Can't find property nx for fix sph/surfacereaction/Langmuir/porous");
  nx = atom->dvector[inx];

  int fny;
  int iny = atom->find_custom("ny", fny);
  if (iny < 0)
    error->all(FLERR,
	       "Can't find property ny for fix sph/surfacereaction/Langmuir/porous");
  ny = atom->dvector[iny];

  int fnz;
  int inz = atom->find_custom("nz", fnz);
  if (inz < 0)
    error->all(FLERR,
	       "Can't find property nz for fix sph/surfacereaction/Langmuir/porous");
  nz = atom->dvector[inz];

  // find the aqueous mass fraction property
  int fxA;
  int ixA = atom->find_custom("xA", fxA);
  if (ixA < 0)
    error->all(FLERR,
               "Can't find property xA for fix sph/surfacereaction/Langmuir/porous");
  xA = atom->dvector[ixA];

  // find the change in aqueous mass fraction concentration property
  int fdxA;
  int idxA = atom->find_custom("dxA", fdxA);
  if (idxA < 0)
    error->all(FLERR,
               "Can't find property dxA for fix sph/surfacereaction/Langmuir/porous");
  dxA = atom->dvector[idxA];

  // find the absorbed mass fraction property
  int fyA;
  int iyA = atom->find_custom("yA", fyA);
  if (iyA < 0)
    error->all(FLERR,
               "Can't find property yA for fix sph/surfacereaction/Langmuir/porous");
  yA = atom->dvector[iyA];

  // find the change in the absorbed mass fraction
  int fdyA;
  int idyA = atom->find_custom("dyA", fdyA);
  if (idyA < 0)
    error->all(FLERR,
               "Can't find property dyA for fix sph/surfacereaction/Langmuir/porous");
  dyA = atom->dvector[idyA];

  // find the normalised concentration of adsorbed species
  int fthetaA;
  int ithetaA = atom->find_custom("thetaA", fthetaA);
  if (ithetaA < 0)
    error->all(FLERR,
               "Can't find property thetaA for fix sph/surfacereaction/Langmuir/porous");
  thetaA = atom->dvector[ithetaA];
}

/* ---------------------------------------------------------------------- */

int FixSPHSurfaceReactionLangmuirPorous::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirPorous::init() {
  dtxA = update->dt;

  // need a full neighbor list, built whenever re-neighboring occurs
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirPorous::init_list(int, NeighList *ptr) {
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirPorous::initial_integrate(int /*vflag*/)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double di, dj;

  double delx, dely, delz;
  double xtmp, ytmp, ztmp;
  double xNij, yNij, zNij, Nij;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, ih, ihsq;
  double rsq, wf, wfd;
  double yAmax;

  double **x = atom->x;
  double **v = atom->v;

  double *rmass = atom->rmass;
  double *rho = atom->rho;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  int *mask = atom->mask;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    // Only update maximum mass fraction for solid particles
    if (itype == 2) {
      imass = rmass[i];
      // Keep the position of the solid particle
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      // Check neighbouring atoms
      jnum = numneigh[i];
      jlist = firstneigh[i];

      // Init value of yAmax to 0
      yAmax = 0.0;

      // Then need to find the closest fluid particles
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        jtype = type[j];

        // Calculate the distance between the particles
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;

        // Check if j is within the support kernel
        if (rsq < h) {
          ih = 1.0/h;

          // kernel function
          if (domain->dimension == 3) {
            wfd = sph_dw_quintic3d(sqrt(rsq)*ih);
            wfd = wfd*ih*ih*ih*ih;
            wf = sph_kernel_quintic3d(sqrt(rsq)*ih)*ih*ih*ih;
          } else {
            wfd = sph_dw_quintic2d(sqrt(rsq)*ih);
            wfd = wfd*ih*ih*ih;
            wf = sph_kernel_quintic2d(sqrt(rsq)*ih)*ih*ih;
          }

          // Perform interaction calculation
          if (jtype == 1) { // only for fluid particles
            jmass = rmass[j];
            // Calculate the normal vector
            xNij = nx[i] + ny[j];
            yNij = ny[i] + ny[j];
            zNij = nz[i] + nz[j];
            // Dot product with position vector
            double nij_rij = xNij*delx + yNij*dely + zNij*delz;
            // Calculate the exchange in concentration
            di = rho[i] / imass;
            dj = rho[j] / jmass;
            yAmax = yAmax + (itype - jtype)*(nij_rij*abs(wfd)*sAmax)/(di*dj*imass);
          } // jtype fluid
        } // loop inside support kernel
      } // for loop jj
      // Calculate the absorbed concentration
      if ((nx[i]==0.0) && (ny[i]==0.0) && (nz[i]==0.0))
        thetaA[i] = yA[i]/(sAmax/rho[i]);
      else
        thetaA[i] = (yAmax == 0.0) ? 0.0 : (yA[i]/yAmax);
    } // if i type is solid
  } // loop through i
}

/* ---------------------------------------------------------------------- */

void FixSPHSurfaceReactionLangmuirPorous::final_integrate() {
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

void FixSPHSurfaceReactionLangmuirPorous::reset_dt() {
  dtxA = update->dt;
}

/* ---------------------------------------------------------------------- */
