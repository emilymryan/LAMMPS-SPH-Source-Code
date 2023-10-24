/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_sph_esph_atom.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSPHEsphAtom::ComputeSPHEsphAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3)
    error->all(FLERR,"Number of arguments for compute sph/esph/atom command != 3");
  if (atom->esph_flag != 1)
    error->all(FLERR,"Compute sph/esph/atom command requires atom_style sph)");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  esphvector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPHEsphAtom::~ComputeSPHEsphAtom()
{
  memory->sfree(esphvector);
}

/* ---------------------------------------------------------------------- */

void ComputeSPHEsphAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"esphvector/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute esphvector/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeSPHEsphAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow evector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(esphvector);
    nmax = atom->nmax;
    esphvector = (double *) memory->smalloc(nmax*sizeof(double),"esphvector/atom:esphvector");
    vector_atom = esphvector;
  }

  double *esph = atom->esph;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              esphvector[i] = esph[i];
      }
      else {
              esphvector[i] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeSPHEsphAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
