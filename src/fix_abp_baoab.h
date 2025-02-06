/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   By Cory Hargus (cory.hargus@gmail.com).
   Adapted from fix_bd_baoab by Guang Shi, July 2016.
   See: B. Leimkuhler and C. Matthews, App. Math. Res. Exp. 2013(1) 2013.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Active Brownian Particles with 2nd order BAOAB integrator.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(abp/BAOAB,FixABPbaoab)

#else

#ifndef LMP_FIX_ABP_BAOAB_H
#define LMP_FIX_ABP_BAOAB_H

#include "fix.h"

namespace LAMMPS_NS {

class FixABPbaoab : public Fix {
 public:
  FixABPbaoab(class LAMMPS *, int, char **);
  virtual ~FixABPbaoab() {}
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void reset_dt();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

 protected:
  class RanMars *random;
  double dt;
  double Dx, Dtheta, gamma_x, omega, v0, seed;
  double xnoise_scale, anoise_scale;
  double xnoise_0, xnoise_1, anoise;
  bool translational_noise;

  double **rnum;
  int nvalues;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
