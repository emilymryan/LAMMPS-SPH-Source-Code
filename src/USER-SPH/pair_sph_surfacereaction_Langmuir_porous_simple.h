#ifdef PAIR_CLASS

PairStyle(sph/surfacereaction/Langmuir/porous/simple, PairSPHSurfaceReactionLangmuirPorousSimple)

#else

#ifndef LMP_PAIR_SPH_SURFACEREACTION_LANGMUIR_POROUS_SIMPLE_H
#define LMP_PAIR_SPH_SURFACEREACTION_LANGMUIR_POROUS_SIMPLE_H

#include "pair.h"

namespace LAMMPS_NS {
  class PairSPHSurfaceReactionLangmuirPorousSimple : public Pair {
  public:
    PairSPHSurfaceReactionLangmuirPorousSimple(class LAMMPS *);
    virtual ~PairSPHSurfaceReactionLangmuirPorousSimple();
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
    double **cut;
    // Mass fraction for aqueous species and absrobed species
    double *xA, *dxA, *yA, *dyA;
    // Normalised mass fraction
    double *thetaA;
    // Binary diffusion coefficient
    double *DA;
    // Adsorption rate coefficient
    double kaA;
    // Micro pore volume
    double micro_pore;
    // Number of adsorption sites needed
    int lambda;
    // phase_support is the size of the inter-phase interaction
    double **phase_support;
    // Check if periodic property is on
    int is_periodic;
    void allocate();
  };
}

#endif
#endif
