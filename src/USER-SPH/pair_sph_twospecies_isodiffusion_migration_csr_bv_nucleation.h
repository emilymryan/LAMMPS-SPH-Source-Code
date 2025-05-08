#ifdef PAIR_CLASS

PairStyle(sph/twospecies/isodiffusion/migration/csr/bv/nucleation, PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation)

#else

#ifndef LMP_PAIR_SPH_TWOSPECIES_ISODIFFUSION_MIGRATION_CSR_BV_NUCLEATION_H
#define LMP_PAIR_SPH_TWOSPECIES_ISODIFFUSION_MIGRATION_CSR_BV_NUCLEATION_H

#include "pair.h"

namespace LAMMPS_NS {
  class PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation : public Pair {
  public:
    PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation(class LAMMPS *);
    virtual ~PairSPHTwospeciesIsodiffusionMigrationCSRbvnucleation();
    void init_style();
    virtual void compute(int, int);
    void settings(int, char **);
    void coeff(int, char **);
    virtual double init_one(int, int);
    virtual double single(int, int, int, int, double, double, double, double &);

    int pack_forward_comm(int, int *, double *, int, int *);
    void unpack_forward_comm(int, int, double *);

    int pack_reverse_comm(int, int, double *);
    void unpack_reverse_comm(int, int *, double *);

  protected:
    // cCeq is the equilibrium concentration for the cations at the anode
    double applied_pot, cCeq, cA_init, cC_init, mass_init, rc, i0, Sc, omega, gc, dGc, Z, N, theta, deltajJ, Nj;
    // cut represents the size of kernel support
    double **cut;
    // phase_support is the size of the inter-phase interaction
    double **phase_support;
    double *cA, *dcA ,*cC, *dcC,*mM, *dmM, *rR, *drR, *nN, *dnN, *local_pot;
    double *muA, *muC, *DA, *DC, *RC,*surf_energy, *nx, *ny, *nz;
    void allocate();
  };
}

#endif
#endif
