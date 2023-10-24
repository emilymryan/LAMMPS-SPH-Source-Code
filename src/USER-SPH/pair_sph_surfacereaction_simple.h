#ifdef PAIR_CLASS

PairStyle(sph/surfacereaction/simple, PairSPHSurfaceReactionSimple)

#else

#ifndef LMP_PAIR_SPH_SURFACEREACTION_SIMPLE_H
#define LMP_PAIR_SPH_SURFACEREACTION_SIMPLE_H

#include "pair.h"

namespace LAMMPS_NS {
  class PairSPHSurfaceReactionSimple : public Pair {
  public:
    PairSPHSurfaceReactionSimple(class LAMMPS *);
    virtual ~PairSPHSurfaceReactionSimple();
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
    // Boolean value to check the periodicity
    int is_periodic;
    // RA is the reaction rate
    double RA;
    // cAeq is the equilibrium concentration
    double cAeq;
    // cut represents the size of kernel support
    double **cut;
    // phase_support is the size of the inter-phase interaction
    double **phase_support;
    double *cA, *dcA, *DA, *mA, *dmA;
    void allocate();
  };
}

#endif
#endif
