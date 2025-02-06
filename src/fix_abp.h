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
   By Cory Hargus (cory.hargus@gmail.com)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Chiral Active Brownian Particles with 1st order Euler–Maruyama integration.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(abp,FixABP)

#else

#ifndef LMP_FIX_ABP_H
#define LMP_FIX_ABP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixABP : public Fix {
 public:
  FixABP(class LAMMPS *, int, char **);
  virtual ~FixABP() {}
  int setmask();
  virtual void init();
  virtual void reset_dt();
  virtual void final_integrate();

 protected:
  class RanPark *random;
  double dt;
  double Dx, Dtheta, gamma_x, omega, v0, seed;
  double xnoise_scale, anoise_scale;
  double xnoise_0, xnoise_1, anoise;
  bool translational_noise;
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
