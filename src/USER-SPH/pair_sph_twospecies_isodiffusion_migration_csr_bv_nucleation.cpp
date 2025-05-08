#include <math.h>
#include <stdlib.h>
#include "pair_sph_twospecies_isodiffusion_migration_csr_bv_nucleation.h"
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

PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // find the anion concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property cA for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  cA = atom->dvector[icA];

  // find the local anion concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (idcA < 0)
    error->all(FLERR,
	       "Can't find property dcA for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  dcA = atom->dvector[idcA];

  // find the cation concentration property
  int fcC;
  int icC = atom->find_custom("cC", fcC);
  if (icC < 0)
    error->all(FLERR,
	       "Can't find property cC for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  cC = atom->dvector[icC];

  // find the local cation concentration property
  int fdcC;
  int idcC = atom->find_custom("dcC", fdcC);
  if (idcC < 0)
    error->all(FLERR,
	       "Can't find property dcC for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  dcC = atom->dvector[idcC];

  // find the anion mobility property
  int fmuA;
  int imuA = atom->find_custom("muA", fmuA);
  if (imuA < 0)
    error->all(FLERR,
	       "Can't find property muA for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  muA = atom->dvector[imuA];

    // find the cation mobility property
  int fmuC;
  int imuC = atom->find_custom("muC", fmuC);
  if (imuC < 0)
    error->all(FLERR,
	       "Can't find property muC for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  muC = atom->dvector[imuC];

  // find the local anion diffusivity constant
  int fDA;
  int iDA = atom->find_custom("DA", fDA);
  if (iDA < 0)
    error->all(FLERR,
	       "Can't find property DA for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  DA = atom->dvector[iDA];

    // find the local cation diffusivity constant
  int fDC;
  int iDC = atom->find_custom("DC", fDC);
  if (iDC < 0)
    error->all(FLERR,
	       "Can't find property DC for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  DC = atom->dvector[iDC];
  
  //find the reaction rate coefficient for cations
  int fRC;
  int iRC = atom->find_custom("RC", fRC);
  if (iRC < 0)
  error->all(FLERR,
             "Can't find property RC for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  RC = atom->dvector[iRC];

  //find the surface energy
  int fsurf_energy;
  int isurf_energy = atom->find_custom("surf_energy",fsurf_energy);
  if (isurf_energy < 0)
    error->all(FLERR,
	       "Can't find property surf_energy sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  surf_energy = atom->dvector[isurf_energy];

  //find the curvature radius property
  int frR;
  int irR = atom->find_custom("rR", frR);
  if (irR < 0)
    error->all(FLERR,
               "Can't find property rR for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  rR = atom->dvector[irR];

  //find the change in radius property                                                                                               
  int fdrR;
  int idrR = atom->find_custom("drR", fdrR);
  if (idrR < 0)
    error->all(FLERR,
               "Can't find property drR for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  drR = atom->dvector[idrR];

  //find the number of nuclei present at the current time step
  int fnN;
  int inN = atom->find_custom("nN", fnN);
  if (inN < 0)
  error->all(FLERR,
  	       "Can't find property nN for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  nN = atom->dvector[inN];

  //find the nucleation rate property
  int fdnN;
  int idnN = atom->find_custom("dnN", fdnN);;
  if (idnN < 0)
  error->all(FLERR,
  	       "Can't find property dnN for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  dnN = atom->dvector[idnN];

  // find the mass of metal property
  int fmM;
  int imM = atom->find_custom("mM", fmM);
  if (imM < 0)
    error->all(FLERR,
	       "Can't find property mM for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  mM = atom->dvector[imM];

  // find the change in mass of metal property
  int fdmM;
  int idmM = atom->find_custom("dmM", fdmM);
  if (idmM < 0)
    error->all(FLERR,
	       "Can't find property dmM for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  dmM = atom->dvector[idmM];

  // find the local potential
  int flocal_pot;
  int ilocal_pot = atom->find_custom("local_pot", flocal_pot);
  if (ilocal_pot < 0)
    error->all(FLERR,
	       "Can't find property local_pot for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  local_pot = atom->dvector[ilocal_pot];

  // find the x-component normal
  int fnx;
  int inx = atom->find_custom("nx", fnx);
  if (inx < 0)
    error->all(FLERR,
	       "Can't find property nx for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  nx = atom->dvector[inx];
  
  // find the y-component normal
  int fny;
  int iny = atom->find_custom("ny", fny);
  if (iny < 0)
    error->all(FLERR,
	       "Can't find property ny for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  ny = atom->dvector[iny];
  
  
  // find the z-component normal
  int fnz;
  int inz = atom->find_custom("nz", fnz);
  if (inz < 0)
    error->all(FLERR,
	       "Can't find property nz for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
  nz = atom->dvector[inz];
 
  // set comm size needed by this pair
  comm_forward = 8;
  comm_reverse = 6;
}

/* ---------------------------------------------------------------------- */

PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::~PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(phase_support);
  }
}

/* ---------------------------------------------------------------------- */
void PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::init_style()
{
  // need a full neighbour list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype; //imol, jmol;
  double delx, dely, delz;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, h, ih, ihsq;
  double r, wf, wfd, D, K;

  double deltaDcA, deltaDcAixx, deltaDcAjxx, deltaDcAiyy, deltaDcAjyy, deltaDcAy, deltaDcAx;
  double deltaDcC, deltaDcCixx, deltaDcCjxx, deltaDcCiyy, deltaDcCjyy, deltaDcCy, deltaDcCx;
  double deltaMcA, deltaMcAi, deltaMcAj;
  double deltaMcC, deltaMcCi, deltaMcCj;
  double deltaDrR;
  double deltaDmM;
  double deltaDnN;
  double ni, nj;
  double xtmp, ytmp;
  double anodePot, overPot;
  double pi, faraday, R, T;
  double transCoefC, transCoefM;
  double bulkcA, bulkcC;
  double mVol, mVol2;
  double Kb, e, Na, nu;
 
  pi = 3.141592654;
  faraday = 96485;
  R = 8.3145;
  T = 300;
  transCoefC = 0.5;
  transCoefM = 0.5;
  anodePot = 0.0;
  overPot = 0.0;
  bulkcA = pow(cA_init, transCoefM);
  bulkcC = pow(cC_init, transCoefC);
  mVol = 13.02*pow(10,6); //um3/umol
  mVol2 = 13.02*pow(10,12);
  e = 1.6e-19;
  Na = 6.02e17; //atoms/umol
  Kb = 8.6e-5; //eV/K

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

  for (ii = 0; ii < inum; ii++)
    {
      dnN[ii] = 0.0;
      drR[ii] = 0.0;
      dmM[ii] = 0.0;
      dcA[ii] = 0.0;
      dcC[ii] = 0.0;
    }
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    // check that we are only doing local and ghost atoms only
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = rmass[i];
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      // check that we are only doing local and ghost atoms only
      jtype = type[j];
      jmass = rmass[j];
    
      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      r = sqrt(delx * delx + dely * dely + delz * delz);
      
      if (r< cut[itype][jtype]) {
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
	
	if ((itype==1) && (jtype==1)) { // fluid-fluid interaction outside of flux
	  // Calculating the particle exchange
	  // Reference: Tartakovsky(2007) - Simulations of reactive transport
	  // and precipitation with sph
	  // The constants give better results...
	  ni = rho[i] / imass;
	  nj = rho[j] / jmass;

	  // Anion mass transport
	  // Anion diffusion
	  deltaDcA = (1.0/(imass*r))*(DA[i]*ni*imass + DA[j]*nj*jmass) \
	    /(ni*nj) * (cA[i] - cA[j])*wfd;
	  
	  // Migration
	  deltaMcAi = muA[i]*cA[i]*ni*imass/0.41615; // 
	  deltaMcAj = muA[i]*cA[j]*nj*jmass/0.41615; // /0.041 // TODO: muA[j] is always 0 when it should be equal to muA[i]
	  // Total migration
	  deltaMcA = (1.0/(imass*r))*(deltaMcAi + deltaMcAj)/(ni*nj) \
	    * (local_pot[i] - local_pot[j])*wfd;
          // Total anion transport
	  dcA[i] += deltaDcA - deltaMcA;
	  
	  // Cation mass transport
	  // Cation diffusion
	  deltaDcC = (1.0/(imass*r))*(DC[i]*ni*imass + DC[j]*nj*jmass) \
	    /(ni*nj) * (cC[i] - cC[j] )*wfd;
	  
	  // Migration
	  deltaMcCi = muC[i]*cC[i]*ni*imass/0.41615;
	  deltaMcCj = muC[i]*cC[j]*nj*jmass/0.41615;// TODO: muA[j] is always 0 when it should be equal to muA[i]	      // Total migration
	  deltaMcC = (1.0/(imass*r))*(deltaMcCi + deltaMcCj)/(ni*nj) \
	    * (local_pot[i] - local_pot[j])*wfd;  
	  // Total cation transport
	  dcC[i] += deltaDcC + deltaMcC;

	  }// solid-fluid interaction
	else if ((itype==1) && (jtype==2)) {  // fluid-solid interaction
	  if (r <= phase_support[itype][jtype]) {
	    
	    // Overpotential calculation
	    overPot = applied_pot - local_pot[i]; //- anodePot
	    //printf("overPot: %.10f \n", overPot); 
	    //Nucleation Rate Equations                                                                                                                        
	    Nj = cC[i]*N; 
            rc = surf_energy[i]*mVol/overPot/(faraday*pow(10,-6)); //Critical Radius of Nucleate
	    //i0 = faraday*RCc*bulkcC*bulkcA*pow(10,-6)*pow(10,-9); //Exchange Current Density
	    i0 = faraday*RC[j]*bulkcC*bulkcA*pow(10,-6)*pow(10,-9);
	    Sc = 2*pi*pow(rc,2);
            omega = Sc*i0/e;
            gc = (32*pi/3)*rc*rc*rc*0.5*Na/mVol;
            dGc = (0.26-applied_pot);
            Z = pow(dGc/(3*pi*Kb*T),0.5)/gc;
            deltaDnN = Nj*omega*Z*exp(-dGc/(Kb*T));
            dnN[i] = deltaDnN; 
	    // Butler - Volmer equations for reaction at solid fluid interface
	    deltaDcC = 1.0*RC[j]*nN[j]*bulkcA*bulkcC*(((cC[i]/cC_init)*exp(transCoefC*(faraday*overPot/R/T + surf_energy[j]*2*mVol2/rR[j]/R/T))) - (exp(-1*transCoefM*(faraday*overPot/R/T + surf_energy[j]*2*mVol2/rR[j]/R/T))));
	    deltaDcC *= fabs(nx[i]) + fabs(nx[j]) + fabs(ny[i]) + fabs(ny[j]);
	    //printf("deltaDcC: %.10f \n", deltaDcC);
	    dcC[i] -= deltaDcC*fabs(wfd);
	  }
	} // fluid-solid interaction   
	else if ((itype==2) && (jtype==1)) { // solid-fluid interaction               
	  if (r<= phase_support[itype][jtype]) {
	   
	    // Overpotential calculation
	    overPot = applied_pot - local_pot[j]; // - anodePot
	   
	    // Radius Calculation
	    deltaDrR = (RC[i]*bulkcA*bulkcC*mVol*pow(10,-15))*(((cC[j]/cC_init)*exp(transCoefC*(faraday*overPot/R/T + surf_energy[i]*2*mVol2/rR[i]/R/T))) - (exp(-1*transCoefM*(faraday*overPot/R/T + surf_energy[i]*2*mVol2/rR[i]/R/T))));
	    drR[i] = deltaDrR;

	    //Nucleation Rate Equations
	    Nj = cC[j]*N;
	    rc = surf_energy[i]*mVol/overPot/(faraday*pow(10,-6)); //Critical Radius of Nucleate
	    //i0 = faraday*RCc*bulkcC*bulkcA*pow(10,-6)*pow(10,-9); //Exchange Current Density
	    i0 = faraday*RC[i]*bulkcC*bulkcA*pow(10,-6)*pow(10,-9);
	    Sc = 2*pi*pow(rc,2);
	    omega = Sc*i0/e;
	    gc = (32*pi/3)*rc*rc*rc*0.5*Na/mVol;
	    dGc = (0.26-applied_pot);
	    Z = pow(dGc/(3*pi*Kb*T),0.5)/gc;
	    deltaDnN = Nj*omega*Z*exp(-dGc/(Kb*T));
	    dnN[i] = deltaDnN;

	    // Bulter-Volmer equations for reaction at solid fluid interface
	    deltaDmM = (jmass)*RC[i]*nN[i]*bulkcA*bulkcC*(((cC[j]/cC_init)*exp(transCoefC*(faraday*overPot/R/T + surf_energy[i]*2*mVol2/rR[i]/R/T))) - (exp(-1*transCoefM*(faraday*overPot/R/T + surf_energy[i]*2*mVol2/rR[i]/R/T))));
	    deltaDmM *= fabs(nx[j]) + fabs(nx[i]) + fabs(ny[j]) + fabs(ny[i]);
	    //printf("deltaDmM: %.10f \n", deltaDmM);
	    dmM[i] += deltaDmM*fabs(wfd);

	  }
	} // solid-fluid interaction 
      } // check if j particle is inside kernel
    } // jj loop
  } // ii loop
  // Communicate the ghost dcA, dcC and dmM to the locally owned atoms
  comm->reverse_comm_pair(this);
}


/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::allocate() {
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

void PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
	       "Illegal number of setting arguments for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation");
}

/* ----------------------------------------------------------------------
   set coefs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::coeff(int narg, char **arg) {
  if (narg != 10)
    error->all(FLERR,"Incorrect number of args for pair_style sph/twospecies/isodiffusion/migration/csr/bv/nucleation coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  // Get kernel size
  double kernel_one = utils::numeric(FLERR,arg[2],false,lmp);
  double phase_one = utils::numeric(FLERR,arg[3],false,lmp);

  // Get the applied potential
  applied_pot = utils::numeric(FLERR,arg[4],false,lmp);

  // Get the initial cation and anion concentrations
  cA_init = utils::numeric(FLERR,arg[5],false,lmp);
  cC_init = utils::numeric(FLERR,arg[6],false,lmp);

  // Get the inital mass 
  mass_init = utils::numeric(FLERR,arg[7],false,lmp);
  
  // Get the value of N (total # of atoms/surface area)
  N = utils::numeric(FLERR,arg[8],false,lmp);
  
  // Get the value of the wetting angle
  theta = utils::numeric(FLERR,arg[9],false,lmp);

  // Get the value of RC
  //RCc = utils::numeric(FLERR,arg[10],false,lmp);

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
    error->all(FLERR,"Incorrect args for pair coefficients sph/twospecies/isodiffusion/migration/csr/constrc/bv/nucleation");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

double PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/twospecies/isodiffusion/migration/csr/bv/nucleation coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  phase_support[j][i] = phase_support[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::single(int i, int j, int itype, int jtype,
					    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::pack_forward_comm(int n, int *list, double *buf,
						    int pbc_flag, int *pbc)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = cA[j];
    buf[m++] = cC[j];
    buf[m++] = mM[j];
    buf[m++] = rR[j];
    buf[m++] = nN[j];
    buf[m++] = local_pot[j];
    buf[m++] = atom->type[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    cA[i] = buf[m++];
    cC[i] = buf[m++];
    mM[i] = buf[m++];
    rR[i] = buf[m++];
    nN[i] = buf[m++];
    local_pot[i] = buf[m++];
    atom->type[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = dcA[i];
    buf[m++] = dcC[i];
    buf[m++] = dmM[i];
    buf[m++] = drR[i];
    buf[m++] = dnN[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    dcA[j] += buf[m++];
    dcC[j] += buf[m++];
    dmM[j] += buf[m++];
    drR[j] += buf[m++];
    dnN[j] += buf[m++];
  }
}
