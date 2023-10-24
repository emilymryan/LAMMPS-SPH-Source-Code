#ifdef PAIR_CLASS

PairStyle(sph/twospecies/isodiffusion/migration/csr/CO2, PairSPHTwospeciesIsodiffusionMigrationCSRCO2)

#else

#ifndef LMP_PAIR_SPH_TWOSPECIES_ISODIFFUSION_MIGRATION_CSR_CO2_H
#define LMP_PAIR_SPH_TWOSPECIES_ISODIFFUSION_MIGRATION_CSR_CO2_H

#include "pair.h"

namespace LAMMPS_NS {
  class PairSPHTwospeciesIsodiffusionMigrationCSRCO2 : public Pair {
  public:
    PairSPHTwospeciesIsodiffusionMigrationCSRCO2(class LAMMPS *);
    virtual ~PairSPHTwospeciesIsodiffusionMigrationCSRCO2();
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
    double keq, applied_pot, cC_init;
    // cut represents the size of kernel support
    double **cut;
    // phase_support is the size of the inter-phase interaction
    double **phase_support;
    double *cA, *dcA ,*cC, *dcC, *local_pot;
    double *muA, *muC, *DA, *DC, *nx, *ny, *nz;
    void allocate();
  };
}

#endif
#endif
