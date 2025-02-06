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

#include <math.h>
#include "math_extra.h"
#include <stdio.h>
#include <string.h>
#include "fix_abp_baoab.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "memory.h"
#include "random_mars.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixABPbaoab::FixABPbaoab(LAMMPS *lmp, int narg, char **arg) :
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

  dynamic_group_allow = 1;
  time_integrate = 1;

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);

  rnum = NULL;

  nvalues = 3;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++){
    rnum[i][0] = random->gaussian();
    rnum[i][1] = random->gaussian();
    rnum[i][2] = random->gaussian();
}
}

/* ---------------------------------------------------------------------- */

int FixABPbaoab::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixABPbaoab::init()
{
  dt = update->dt;
  atom->add_custom("thetavec", 1, 0);
  xnoise_scale = sqrt(0.5 * Dx * dt);  // Factor of 0.5 for BAOAB
  anoise_scale = sqrt(0.5 * Dtheta * dt);  // Factor of 0.5 for BAOAB
  // dtv = update->dt;  // timestep
  // dtf = t_period * force->ftm2v;
  // gfactor = sqrt(0.5*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixABPbaoab::initial_integrate(int vflag)
{
  double dtm;
  double randf;
  double rx, ry, rtheta;

  // update v and x of atoms in group

  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  int flag, cols;
  int index = atom->find_custom("thetavec", flag, cols);
  double *theta = atom->dvector[index];
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  // in LJ units, t_target is given in kB*T/epsilon
  // Note, to avoid noisy correlation functions, we only add the
  // deterministic forces to the velocities. When computing the VACF later,
  // we need to add kB*T/gamma * delta(t), i.e. a delta function with
  // magnitude scaled by the bath diffusivity D = kB*T/gamma.
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (translational_noise) {
      rx = random->gaussian();
      ry = random->gaussian();
        // xnoise_0 = xnoise_scale * random->gaussian();
        // xnoise_1 = xnoise_scale * random->gaussian();
      }
      rtheta = random->gaussian();
      v[i][0] = (f[i][0] / gamma_x) + v0 * cos(theta[i]); //+ xnoise_0;  For convenience, remove
      v[i][1] = (f[i][1] / gamma_x) + v0 * sin(theta[i]); //+ xnoise_1;  noise from velocities.

      dtm = dt / mass[type[i]];
      randf = sqrt(mass[type[i]]) * xnoise_scale;
      if (translational_noise) {
        x[i][0] += (dtm * v[i][0]) + randf * (rx + rnum[i][0]);
        x[i][1] += (dtm * v[i][1]) + randf * (ry + rnum[i][1]);
      }
      else {
        x[i][0] += (dtm * v[i][0]);
        x[i][1] += (dtm * v[i][1]);
      }

      theta[i] += omega*dt + anoise_scale * (rtheta + rnum[i][2]);  // Integrate director angle

      if (theta[i] > MY_PI){  // Map theta back into [-pi, pi) if needed
        theta[i] -= MY_2PI;
      }
      if (theta[i] < -MY_PI){
        theta[i] += MY_2PI;
      }

      rnum[i][0] = rx;
      rnum[i][1] = ry;
      rnum[i][2] = rtheta;
  }
}

/* ---------------------------------------------------------------------- */

void FixABPbaoab::reset_dt()
{
  dt = update->dt;
}

/* ----------------------------------------------------------------------
   allocate atom-based array for rnum
------------------------------------------------------------------------- */

void FixABPbaoab::grow_arrays(int nmax)
{
  memory->grow(rnum,nmax,3,"fix_bd_baoab:rnum");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixABPbaoab::copy_arrays(int i, int j, int delflag)
{
  for (int m = 0; m < nvalues; m++)
    rnum[j][m] = rnum[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixABPbaoab::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nvalues; m++) buf[m] = rnum[i][m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixABPbaoab::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nvalues; m++) rnum[nlocal][m] = buf[m];
  return nvalues;
}
