#include <iostream>
#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "pair_sph_surfacereaction_normal.h"
#include "sph_kernel_quintic.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHSurfaceReactionNormal::PairSPHSurfaceReactionNormal(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // find the normal vector
  int fnx;
  int inx = atom->find_custom("nx", fnx);
  if (inx < 0)
    error->all(FLERR,
	       "Can't find property nx for pair_style sph/surfacereaction/normal");
  nx = atom->dvector[inx];

  int fny;
  int iny = atom->find_custom("ny", fny);
  if (iny < 0)
    error->all(FLERR,
	       "Can't find property ny for pair_style sph/surfacereaction/normal");
  ny = atom->dvector[iny];

  int fnz;
  int inz = atom->find_custom("nz", fnz);
  if (inz < 0)
    error->all(FLERR,
	       "Can't find property nz for pair_style sph/surfacereaction/normal");
  nz = atom->dvector[inz];

  // set comm size needed by this Pair
  comm_forward = 3;
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHSurfaceReactionNormal::~PairSPHSurfaceReactionNormal() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
   ------------------------------------------------------------------------- */

void PairSPHSurfaceReactionNormal::init_style() {
  // need a full neighbor list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHSurfaceReactionNormal::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype;
  double delx, dely, delz;
  double r, h, ih;
  int *jlist;
  double di, dj;
  double wfd;

  const int ndim = domain->dimension;
  double eij[ndim];
  // neighbor list variables
  int inum, *ilist, *numneigh, **firstneigh;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double *rho = atom->rho;
  int *type = atom->type;
  double *rmass = atom->rmass;

  // check consistency of pair coefficients
  if (first) {
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = 1; i <= atom->ntypes; i++) {
        if (cutsq[i][j] > 0.0) {
          if (!setflag[i][i] || !setflag[j][j]) {
            if (comm->me == 0) {
              printf(
                     "SPH particle types %d and %d interact, but not all of their single particle properties are set.\n",
                     i, j);
            }
          }
        }
      }
    }
    first = 0;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // recompute density
  // we use a full neighborlist here

  // initialize color gradient with zeros
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    nx[i] = 0.0;
    ny[i] = 0.0;
    nz[i] = 0.0;
  } // ii loop

  // add density at each atom via kernel function overlap
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    di = rho[i]/rmass[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      // Calculate the vector rij
      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      r = sqrt(delx*delx + dely*dely + delz*delz);

      if (r < cut[itype][jtype]) {
        h = cut[itype][jtype];
        ih = 1.0/h;

        // Calculate kernel function value
        if (domain->dimension == 3) {
          // Quintic spline
          wfd = sph_dw_quintic3d(r*ih);
          wfd = wfd*ih*ih*ih*ih;
        } else {
          wfd = sph_dw_quintic2d(r*ih);
          wfd = wfd*ih*ih*ih;
        }
        dj = rho[j]/rmass[j];
        double dphi = dj*(itype - jtype)*wfd;
        nx[i] += dphi*delx;
        ny[i] += dphi*dely;
        nz[i] += dphi*delz;
      } // rsq < cutsq
    } // jj loop
    // Return the normalised value
    double nmag = sqrt(nx[i]*nx[i] + ny[i]*ny[i] + nz[i]*nz[i]);
    nx[i] /= nmag;
    ny[i] /= nmag;
    nz[i] /= nmag;
  } // ii loop

  // communicate surface normal
  comm->forward_comm_pair(this);
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairSPHSurfaceReactionNormal::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairSPHSurfaceReactionNormal::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
               "Illegal number of setting arguments for pair_style sph/surfacereaction/normal");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairSPHSurfaceReactionNormal::coeff(int narg, char **arg) {
  if (narg != 3)
    error->all(FLERR,"Incorrect number of args for sph/surfacereaction/normal coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);

  double cut_one = force->numeric(FLERR, arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

double PairSPHSurfaceReactionNormal::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/surfacereaction/normal coeffs are not set");
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHSurfaceReactionNormal::single(int i, int j, int itype, int jtype, double rsq,
                                            double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHSurfaceReactionNormal::pack_comm(int n, int *list, double *buf, int pbc_flag,
                                            int *pbc) {
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = nx[j];
    buf[m++] = ny[j];
    buf[m++] = nz[j];
  }
  return 3;
}

/* ---------------------------------------------------------------------- */

void PairSPHSurfaceReactionNormal::unpack_comm(int n, int first, double *buf) {
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    nx[i] = buf[m++];
    ny[i] = buf[m++];
    nz[i] = buf[m++];
  }
}
