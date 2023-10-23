#ifdef PAIR_CLASS

PairStyle(sph/twospecies/isodiffusion/migration/csr/constRC, PairSPHTwospeciesIsodiffusionMigrationCSRConstRC)

#else

#ifndef LMP_PAIR_SPH_TWOSPECIES_ISODIFFUSION_MIGRATION_CSR_CONSTRC_H
#define LMP_PAIR_SPH_TWOSPECIES_ISODIFFUSION_MIGRATION_CSR_CONSTRC_H

#include "pair.h"

namespace LAMMPS_NS {
  class PairSPHTwospeciesIsodiffusionMigrationCSRConstRC : public Pair {
  public:
    PairSPHTwospeciesIsodiffusionMigrationCSRConstRC(class LAMMPS *);
    virtual ~PairSPHTwospeciesIsodiffusionMigrationCSRConstRC();
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
    double conc_to_charge, applied_pot, cCeq, electrolyte_start, electrolyte_end, RC;
    // cut represents the size of kernel support
    double **cut;
    // phase_support is the size of the inter-phase interaction
    double **phase_support;
    double *cA, *dcA ,*cC, *dcC,*mM, *dmM, *local_pot;
    double *muA, *muC, *DA, *DC, *nx, *ny, *nz;
    void allocate();
  };
}

#endif
#endif
