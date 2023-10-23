#include <math.h>
#include <stdlib.h>
#include "pair_sph_surfacereaction_simple_discharge.h"
#include "sph_kernel_quintic.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "domain.h"
#include <iostream>
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHSurfaceReactionSimpleDischarge::PairSPHSurfaceReactionSimpleDischarge(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property cA for pair_style sph/surfacereaction/simple/discharge");
  cA = atom->dvector[icA];

  // find the local concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (idcA < 0)
    error->all(FLERR,
	       "Can't find property dcA for pair_style sph/surfacereaction/simple/discharge");
  dcA = atom->dvector[idcA];

  // find the local diffusivity constant
  int fDAx;
  int iDAx = atom->find_custom("DAx", fDAx);
  if (iDAx < 0)
    error->all(FLERR,
	       "Can't find property DAx for pair_style sph/surfacereaction/simple/discharge");
  DAx = atom->dvector[iDAx];

  // find the local diffusivity constant y
  int fDAy;
  int iDAy = atom->find_custom("DAy", fDAy);
  if (iDAy < 0)
    error->all(FLERR,
	       "Can't find property DAy for pair_style sph/surfacereaction/simple/discharge");
  DAy = atom->dvector[iDAy];

  // find the change in mass of A property
  int fdmA;
  int idmA = atom->find_custom("dmA", fdmA);
  if (idmA < 0)
    error->all(FLERR,
	       "Can't find property dmA for pair_style sph/surfacereaction/simple/discharge");
  dmA = atom->dvector[idmA];

  // find the mass of A property
  int fmA;
  int imA = atom->find_custom("mA", fmA);
  if (imA < 0)
    error->all(FLERR,
	       "Can't find property mA for pair_style sph/surfacereaction/simple/discharge");
  mA = atom->dvector[imA];

  // set comm size needed by this pair
  comm_forward = 3;
  comm_reverse = 2;
}

/* ---------------------------------------------------------------------- */

PairSPHSurfaceReactionSimpleDischarge::~PairSPHSurfaceReactionSimpleDischarge() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(phase_support);
  }
}

/* ---------------------------------------------------------------------- */
void PairSPHSurfaceReactionSimpleDischarge::init_style()
{
  // need a full neighbour list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHSurfaceReactionSimpleDischarge::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double delx, dely, delz;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, h, ih, ihsq;
  double r, wf, wfd, D, K, deltacA, deltacAixx, deltacAjxx, deltacAiyy, deltacAjyy, deltacAy, deltacAx;

  double ni, nj;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double *rmass = atom->rmass;
  double *rho = atom->rho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Communicate the local cA to the ghost atoms
  comm->forward_comm_pair(this);

  // loop over neighbors of my atoms and do heat diffusion
  for (ii = 0; ii < inum; ii++)
    {
      dmA[ii] = 0.0;
      dcA[ii] = 0.0;
    }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    // check that we are only doing local and ghost atoms only
    itype = type[i];
    // check if the i particles is within the domain
    if ((is_periodic==1) || (not (x[i][0] < domain->boxlo[0] || x[i][0] > domain->boxhi[0] ||
                                  x[i][1] < domain->boxlo[1] || x[i][1] > domain->boxhi[1] ||
                                  x[i][2] < domain->boxlo[2] || x[i][2] > domain->boxhi[2]))) {
      jlist = firstneigh[i];
      jnum = numneigh[i];

      imass = rmass[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        // check that we are only doing local and ghost atoms only
        jtype = type[j];
        jmass = rmass[j];

        // check if the j particles is within the domain
        // also check if the periodicity of reaction is required
        if ((is_periodic==1) || (not (x[j][0] < domain->boxlo[0] || x[j][0] > domain->boxhi[0] ||
                                      x[j][1] < domain->boxlo[1] || x[j][1] > domain->boxhi[1] ||
                                      x[j][2] < domain->boxlo[2] || x[j][2] > domain->boxhi[2]))) {
          delx = x[i][0] - x[j][0];
          dely = x[i][1] - x[j][1];
          delz = x[i][2] - x[j][2];
          r = sqrt(delx * delx + dely * dely + delz * delz);

          if (r < cut[itype][jtype]) {
            h = cut[itype][jtype];
            ih = 1.0/h;

            // kernel function
            if (domain->dimension == 3) {
              wfd = sph_dw_quintic3d(r*ih);
              wfd = wfd*ih*ih*ih*ih;
              wf = sph_kernel_quintic3d(r*ih)*ih*ih*ih;
            } else {
              wfd = sph_dw_quintic2d(r*ih);
              wfd = wfd*ih*ih*ih;
              wf = sph_kernel_quintic2d(r*ih)*ih*ih;
            }

            if ((itype==1) && (jtype==1)) { // fluid-fluid interaction in electrolyte
              // Calculating the particle exchange
              // Reference: Tartakovsky(2007) - Simulations of reactive transport
              // and precipitation with sph
              // The constants give better results...
              ni = rho[i] / imass;
              nj = rho[j] / jmass;

	      deltacAixx = (DAx[i]*ni*imass) * (delx * delx) / (r*r) * 0.5;
              deltacAjxx = (DAx[j]*nj*jmass) * (delx * delx) / (r*r) * 0.5;
              deltacAx = deltacAixx + deltacAjxx;

              deltacAiyy = (DAy[i]*ni*imass) * (dely * dely) / (r*r) * 0.5;
              deltacAjyy = (DAy[j]*nj*jmass) * (dely * dely) / (r*r) * 0.5;
              deltacAy = deltacAiyy + deltacAjyy;

              deltacA = (1.0/(imass*r))*(deltacAx + deltacAy)/(nj*nj) * (cA[i] - cA[j])*wfd;

              // deltacA = (1.0/(imass*r))*
              //   ((DA[i]*ni*imass + DA[j]*nj*jmass)/(ni*nj))*(cA[i] - cA[j])*wfd;
              // dcA[i] = dcA[i] + deltacA;
            } // fluid-fluid interaction
	    else if ((itype==1) && (jtype==2)) { // fluid-solid interaction
	      if (r <= phase_support[itype][jtype]) {
		deltacA = 1.0*RA*(cA[i] - cAeq);     //
		dcA[i] = dcA[i] + deltacA;
	      }
	    } // fluid-solid interaction
	    else if ((itype==2) && (jtype==1)) { // solid-fluid interaction
	      if (r <= phase_support[itype][jtype]) {
		dmA[i] = dmA[i] - (jmass)*RA*(cA[j] - cAeq);//
	      }
	    } // solid-fluid interaction
	  } // check if j particle is inside kernel
        } // check if j particle is inside domain
      } // jj loop
    } // check i atom is inside domain
  } // ii loop

  // Communicate the ghost dcA and dmA to the locally owned atoms
  comm->reverse_comm_pair(this);
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairSPHSurfaceReactionSimpleDischarge::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(phase_support, n + 1, n + 1, "pair:phase_support");
}

/* e----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairSPHSurfaceReactionSimpleDischarge::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
	       "Illegal number of setting arguments for pair_style sph/concprecipitation");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairSPHSurfaceReactionSimpleDischarge::coeff(int narg, char **arg) {
  if (narg != 7)
    error->all(FLERR,"Incorrect number of args for pair_style sph/surfacereaction/simple/discharge coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR,arg[0],1, atom->ntypes, ilo, ihi,error);
  utils::bounds(FLERR,arg[1],1, atom->ntypes, jlo, jhi,error);

  // Get kernel size
  double kernel_one = utils::numeric(FLERR,arg[2],false,lmp);
  double phase_one = utils::numeric(FLERR,arg[3],false,lmp);

  // Get the reaction parameters
  cAeq = utils::numeric(FLERR,arg[4],false,lmp);
  RA = utils::numeric(FLERR,arg[5],false,lmp);

  // Check if periodicity for transport is allowed
  is_periodic = utils::numeric(FLERR,arg[6],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = kernel_one;
      phase_support[i][j] = phase_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients sph/surfacereaction/simple/discharge");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

double PairSPHSurfaceReactionSimpleDischarge::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/surfacereaction/simple/discharge coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  phase_support[j][i] = phase_support[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHSurfaceReactionSimpleDischarge::single(int i, int j, int itype, int jtype,
					    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHSurfaceReactionSimpleDischarge::pack_forward_comm(int n, int *list, double *buf,
						    int pbc_flag, int *pbc)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = cA[j];
    buf[m++] = mA[j];
    buf[m++] = atom->type[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHSurfaceReactionSimpleDischarge::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    cA[i] = buf[m++];
    mA[i] = buf[m++];
    atom->type[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int PairSPHSurfaceReactionSimpleDischarge::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = dcA[i];
    buf[m++] = dmA[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHSurfaceReactionSimpleDischarge::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    dcA[j] += buf[m++];
    dmA[j] += buf[m++];
  }
}
