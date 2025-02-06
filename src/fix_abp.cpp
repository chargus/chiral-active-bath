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


#include "stdio.h"
#include "string.h"
#include "comm.h"
#include "fix_abp.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "math.h"
#include "random_park.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

// example command
// fix abp DX DTHETA GAMMA_X OMEGA V0 SEED
FixABP::FixABP(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"abp") != 0 && narg < 9)
    error->all(FLERR,"Illegal fix abp command");

  Dx = utils::numeric(FLERR,arg[3],false,lmp);
  Dtheta = utils::numeric(FLERR,arg[4],false,lmp);
  gamma_x = utils::numeric(FLERR,arg[5],false,lmp);
  omega = utils::numeric(FLERR,arg[6],false,lmp);
  v0 = utils::numeric(FLERR,arg[7],false,lmp);
  int seed = utils::inumeric(FLERR,arg[8],false,lmp);
  
  if (gamma_x==0)  {
    translational_noise=0;  // No translational noise
  } else {
    translational_noise=1;
  }

  // allocate the random number generator
  random = new RanPark(lmp,seed + comm->me); // use processor-unique seed
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixABP::setmask()
{
  int mask = 0;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixABP::init()
{
  dt = update->dt;
}

void FixABP::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixABP::final_integrate()
{
  double noise_0,noise_1,noise_2;
  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  // double *theta = atom->theta;  // Angle of director (using dipole moment)
  int flag, cols;
  int index = atom->find_custom("thetavec", flag, cols);
  double *theta = atom->dvector[index];
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  xnoise_scale = sqrt(2 * Dx * dt);
  anoise_scale = sqrt(2 * Dtheta * dt);
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      // First, update theta.
      anoise   = anoise_scale * random->gaussian();
      theta[i] += omega*dt + anoise;  // First, integrate the director angle
      if (theta[i] > MY_PI){  // Map theta back into [-pi, pi) if needed
        theta[i] -= MY_2PI;
      }
      if (theta[i] < -MY_PI){
        theta[i] += MY_2PI;
      }

      // Then update velocities. Exclude translational noise from the velocity
      // for the convenience of averaging, e.g. in a VACF calculation. Thus,
      // if used, it must be added later, e.g. as a kB*T*delta(t) contribution
      // to the VACF.
      v[i][0] = (f[i][0] / gamma_x) + v0 * cos(theta[i]); //+ xnoise_0;  For convenience, remove
      v[i][1] = (f[i][1] / gamma_x) + v0 * sin(theta[i]); //+ xnoise_1;  noise from velocities.

      // Update the positions, either with or without translational noise.
      if (translational_noise) {
        xnoise_0 = xnoise_scale * random->gaussian();
        xnoise_1 = xnoise_scale * random->gaussian();
        x[i][0] += (dt * v[i][0]) + xnoise_0;
        x[i][1] += (dt * v[i][1]) + xnoise_1;
      }
      else {
        x[i][0] += (dt * v[i][0]);
        x[i][1] += (dt * v[i][1]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

