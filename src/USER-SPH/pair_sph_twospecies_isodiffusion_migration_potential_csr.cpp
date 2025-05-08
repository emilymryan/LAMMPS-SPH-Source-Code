#include <math.h>
#include <stdlib.h>
#include "pair_sph_twospecies_isodiffusion_migration_potential_csr.h"
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
#include <cstdio>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::PairSPHTwospeciesIsodiffusionMigrationPotentialCSR(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // find the anion concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property cA for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  cA = atom->dvector[icA];

  // find the local anion concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (idcA < 0)
    error->all(FLERR,
	       "Can't find property dcA for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  dcA = atom->dvector[idcA];

  // find the cation concentration property
  int fcC;
  int icC = atom->find_custom("cC", fcC);
  if (icC < 0)
    error->all(FLERR,
	       "Can't find property cC for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  cC = atom->dvector[icC];

  // find the local cation concentration property
  int fdcC;
  int idcC = atom->find_custom("dcC", fdcC);
  if (idcC < 0)
    error->all(FLERR,
	       "Can't find property dcC for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  dcC = atom->dvector[idcC];

  // find the anion mobility property
  int fmuA;
  int imuA = atom->find_custom("muA", fmuA);
  if (imuA < 0)
    error->all(FLERR,
	       "Can't find property muA for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  muA = atom->dvector[imuA];

    // find the cation mobility property
  int fmuC;
  int imuC = atom->find_custom("muC", fmuC);
  if (imuC < 0)
    error->all(FLERR,
	       "Can't find property muC for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  muC = atom->dvector[imuC];

  // find the local anion diffusivity constant
  int fDA;
  int iDA = atom->find_custom("DA", fDA);
  if (iDA < 0)
    error->all(FLERR,
	       "Can't find property DA for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  DA = atom->dvector[iDA];

    // find the local cation diffusivity constant
  int fDC;
  int iDC = atom->find_custom("DC", fDC);
  if (iDC < 0)
    error->all(FLERR,
	       "Can't find property DC for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  DC = atom->dvector[iDC];

  // find the reaction rate coefficient for the cations
  int fRC;
  int iRC = atom->find_custom("RC", fRC);
  if (iRC < 0)
    error->all(FLERR,
	       "Can't find property RC for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  RC = atom->dvector[iRC];

  // find the mass of metal property
  int fmM;
  int imM = atom->find_custom("mM", fmM);
  if (imM < 0)
    error->all(FLERR,
	       "Can't find property mM for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  mM = atom->dvector[imM];

  // find the change in mass of metal property
  int fdmM;
  int idmM = atom->find_custom("dmM", fdmM);
  if (idmM < 0)
    error->all(FLERR,
	       "Can't find property dmM for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  dmM = atom->dvector[idmM];

  // find the local potential
  int flocal_pot;
  int ilocal_pot = atom->find_custom("local_pot", flocal_pot);
  if (ilocal_pot < 0)
    error->all(FLERR,
	       "Can't find property local_pot for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  local_pot = atom->dvector[ilocal_pot];

  // find the next local potential
  int fnext_local_pot;
  int inext_local_pot = atom->find_custom("next_local_pot", fnext_local_pot);
  if (inext_local_pot < 0)
    error->all(FLERR,
	       "Can't find property next_local_pot for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  next_local_pot = atom->dvector[inext_local_pot];

  // find the x-component normal
  int fnx;
  int inx = atom->find_custom("nx", fnx);
  if (inx < 0)
    error->all(FLERR,
	       "Can't find property nx for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  nx = atom->dvector[inx];
  
  // find the next x-component normal
  int fdnx;
  int idnx = atom->find_custom("dnx", fdnx);
  if (idnx < 0)
    error->all(FLERR,
	       "Can't find property dnx for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  dnx = atom->dvector[idnx];
  // find the y-component normal
  int fny;
  int iny = atom->find_custom("ny", fny);
  if (iny < 0)
    error->all(FLERR,
	       "Can't find property ny for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  ny = atom->dvector[iny];
  
  // find the next y-component normal
  int fdny;
  int idny = atom->find_custom("dny", fdny);
  if (idny < 0)
    error->all(FLERR,
	       "Can't find property dny for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  dny = atom->dvector[idny];
  
  // find the z-component normal
  int fnz;
  int inz = atom->find_custom("nz", fnz);
  if (inz < 0)
    error->all(FLERR,
	       "Can't find property nz for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  nz = atom->dvector[inz];
  
  // find the next z-component normal
  int fdnz;
  int idnz = atom->find_custom("dnz", fdnz);
  if (idnz < 0)
    error->all(FLERR,
	       "Can't find property dnz for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
  dnz = atom->dvector[idnz];

  // set comm size needed by this pair
  comm_forward = 9;
  comm_reverse = 6;
}

/* ---------------------------------------------------------------------- */

PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::~PairSPHTwospeciesIsodiffusionMigrationPotentialCSR() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(phase_support);
  }
}

/* ---------------------------------------------------------------------- */
void PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::init_style()
{
  // need a full neighbour list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double delx, dely, delz, rdelx, rdely, rdelz;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, h, ih, ihsq;
  double r, wf, wfd, wfdx, wfdy, wfdz, nwfdx, nwfdy, nwfdz, D, K;

  double deltaDcA, deltaDcC;
  double deltaMcA, deltaMcAi, deltaMcAj;
  double deltaMcC, deltaMcCi, deltaMcCj;
  double deltadmM;
  double chargedensityi, chargedensityj;
  double  pi;
  double ni, nj;
  double dn_mag;

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

  pi = 3.1415927;
  dn_mag = 0.0;
  
  // Communicate the local cA to the ghost atoms
  comm->forward_comm_pair(this);

  for (ii = 0; ii < inum; ii++)
    {
      dmM[ii] = 0.0;
      dcA[ii] = 0.0;
      dcC[ii] = 0.0;
      next_local_pot[ii] = 0.0;
      dnx[ii] = 0.0;
      dny[ii] = 0.0;
      dnz[ii] = 0.0;
    }
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    // check that we are only doing local and ghost atoms only
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    imass = rmass[i];
    // Calculate charge density &&at i
    chargedensityi = conc_to_charge * (cC[i] - cA[i]);
    next_local_pot[i] =  chargedensityi + applied_pot * ( x[i][1] - electrolyte_start)/(electrolyte_end-electrolyte_start);
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      // check that we are only doing local and ghost atoms only
      jtype = type[j];
      jmass = rmass[j];

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      rdelx = fabs(delx);
      rdely = fabs(dely);
      rdelz = fabs(delz);

      r = sqrt(delx * delx + dely * dely + delz * delz);

      
      if (r+0.0000000001 < cut[itype][jtype]) {
	h = cut[itype][jtype];
	ih = 1.0/h;
	
	// kernel function
	if (domain->dimension == 3) {
	  printf("3D sims not vetted");
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
	
	if ((itype==1) && (jtype==1)) { // fluid-fluid interaction outside of flux
	  // Calculating the particle exchange
	  // Reference: Tartakovsky(2007) - Simulations of reactive transport
	  // and precipitation with sph
	  // The constants give better results...
	  ni = rho[i] / imass;
	  nj = rho[j] / jmass;
	  
	  // Anion mass transport
	  // Anion diffusion
	  deltaDcA = (1.0/(2.0*imass*r))*(DA[i]*ni*imass + DA[j]*nj*jmass)/(ni*nj) * (cA[i] - cA[j])*wfd;
	  
	  // Migration
	  deltaMcAi = muA[i]*cA[i]*ni*imass;
	  deltaMcAj = muA[i]*cA[j]*nj*jmass; // TODO: muA[j] is always 0 when it should be equal to muA[i]
	  // Total migration
	  deltaMcA = (1.0/(2.0*imass*r))*(deltaMcAi + deltaMcAj)/(ni*nj) * (local_pot[i] - local_pot[j])*wfd;
	  // Total anion transport
	  dcA[i] += deltaDcA - deltaMcA;
	  
	  // Cation mass transport
	  deltaDcC = (1.0/(2.0*imass*r))*(DC[i]*ni*imass + DC[j]*nj*jmass)/(ni*nj) * (cC[i] - cC[j] )*wfd;
	  
	  // Migration
	  deltaMcCi = muC[i]*cC[i]*ni*imass;
	  deltaMcCj = muC[i]*cC[j]*nj*jmass;// TODO: muA[j] is always 0 when it should be equal to muA[i]	      // Total migration
	  deltaMcC = (1.0/(imass*r))*(deltaMcCi + deltaMcCj)/(ni*nj) * (local_pot[i] - local_pot[j])*wfd;  
	  // Total cation transport
	  dcC[i] += deltaDcC + deltaMcC;
	  
	  // Calculating the charge density and local potential for the  time step
	  chargedensityj = conc_to_charge * (cC[j] - cA[j]);
	  // Net charge density for next time step
	  next_local_pot[i] = next_local_pot[i]  - chargedensityj / (4 * pi) / r;
	} // fluid-fluid interaction
	else if ((itype==1) && (jtype==2)) {  // fluid-solid interaction
	  if (r+0.00000000001 <= phase_support[itype][jtype]) {
	    deltaDcC = 1.0*RC[j]*(cC[i] - cCeq);
	    deltaDcC *= fabs((nx[j] + nx[i]) + (ny[j] + ny[i]));
	    dcC[i] -= deltaDcC;
	    // Normals calculation for next step
	    dnx[i] += nwfdx;
	    dny[i] += nwfdy;
	    // j particle (solid) contributes 0v to net charge density
	    next_local_pot[i] = next_local_pot[i];
	  }
	} // fluid-solid interaction   
	else if ((itype==2) && (jtype==1)) { // solid-fluid interaction               
	  if (r+0.00000000001 <= phase_support[itype][jtype]) {
	    deltadmM = (jmass)*RC[i]*(cC[j] - cCeq);
	    deltadmM *= fabs((nx[j] + nx[i]) + (ny[j] + ny[i]));
	    dmM[i] += deltadmM;
	    
	    // Normals calculation for next step
	    dnx[i] -= nwfdx;
	    dny[i] -= nwfdy;
	    // Potential in solid is 0v
	    next_local_pot[i] = 0;
	  }
	} // solid-fluid interaction                                             
      } // check if j particle is inside kernel
    } // jj loop
    dn_mag = sqrt(dnx[i]*dnx[i] + dny[i]*dny[i]);// + dnz[i]*dnz[i]);
    if(dn_mag==0) {
      dnx[i] = 0;
      dny[i] = 0;
    }
    else {
      dnx[i] = dnx[i]/dn_mag;
      dny[i] = dny[i]/dn_mag;
    }
  } // ii loop

  // Communicate the ghost dcA, dcC and dmM to the locally owned atoms
  comm->reverse_comm_pair(this);
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::allocate() {
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

void PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
	       "Illegal number of setting arguments for pair_style sph/twospecies/isodiffusion/migration/potential/csr");
}

/* ----------------------------------------------------------------------
   set coefs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::coeff(int narg, char **arg) {
  if (narg != 9)
    error->all(FLERR,"Incorrect number of args for pair_style sph/twospecies/isodiffusion/migration/potential/csr  coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  // Get kernel size
  double kernel_one = utils::numeric(FLERR,arg[2],false,lmp);
  double phase_one = utils::numeric(FLERR,arg[3],false,lmp);

  // Get the conversion from concentration to local charge
  conc_to_charge = utils::numeric(FLERR,arg[4],false,lmp);

  // Get the applied potential
  applied_pot = utils::numeric(FLERR,arg[5],false,lmp);

  // Get the position of the anode end
  electrolyte_start = utils::numeric(FLERR,arg[6],false,lmp);

  // Get the position of the boundary condition
  electrolyte_end = utils::numeric(FLERR,arg[7],false,lmp);
  
  // Get the cation equilibrium concentration
  cCeq = utils::numeric(FLERR,arg[8],false,lmp);

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
    error->all(FLERR,"Incorrect args for pair coefficients sph/twospecies/isodiffusion/migration/potential/csr");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

double PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/twospecies/isodiffusion/migration/potential/csr coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  phase_support[j][i] = phase_support[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::single(int i, int j, int itype, int jtype,
					    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::pack_forward_comm(int n, int *list, double *buf,
						    int pbc_flag, int *pbc)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = cA[j];
    buf[m++] = cC[j];
    buf[m++] = mM[j];
    buf[m++] = RC[j];
    buf[m++] = local_pot[j];
    buf[m++] = atom->type[j];
    buf[m++] = nx[i];
    buf[m++] = ny[i];
    buf[m++] = nz[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    cA[i] = buf[m++];
    cC[i] = buf[m++];
    mM[i] = buf[m++];
    RC[i] = buf[m++];
    local_pot[i] = buf[m++];
    atom->type[i] = buf[m++];
    nx[i] = buf[m++];
    ny[i] = buf[m++];
    nz[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = dcA[i];
    buf[m++] = dcC[i];
    buf[m++] = dmM[i];
    buf[m++] = dnx[i];
    buf[m++] = dny[i];
    buf[m++] = dnz[i];

  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHTwospeciesIsodiffusionMigrationPotentialCSR::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    dcA[j] += buf[m++];
    dcC[j] += buf[m++];
    dmM[j] += buf[m++];
    dnx[i] += buf[m++];
    dny[i] += buf[m++];
    dnz[i] += buf[m++];
  }
}
