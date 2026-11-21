/**
# Capillary breakup of a randomly corrugated liquid cylinder (axisymmetric)

Usage: ./ligament MAXLEVEL SEED [TOLERANCE]

Units: mean radius R = 1, liquid density 1, surface tension 1, so that the
capillary time sqrt(rho R^3/sigma) is 1. x is the axial direction, y is the
radial direction, the axis of symmetry is the bottom boundary.

Volume convention: with axi.h, dv() = Delta^2*y, i.e. the volume per radian.
The true volume is 2*pi*sum(f*dv()). All volumes written by this code are
true volumes (the factor 2*pi is included). */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"

#define NMODES 25
#define RMS 0.08
#define TEND 14.
#define MINLEVEL 6
#define W 2. // initial mean diameter

int MAXLEVEL = 10, seed = 1;
double phase[NMODES];
double vol0, dtmin = HUGE, t_pinch = -1.;
long cells_max = 0;
int mgp_max = 0;

double radius (double x)
{
  double a = RMS*sqrt(2./NMODES), eps = 0.;
  for (int n = 1; n <= NMODES; n++)
    eps += a*cos(2.*pi*n*x/L0 + phase[n-1]);
  return 1. + eps;
}

int main (int argc, char * argv[])
{
  if (argc > 1) MAXLEVEL = atoi (argv[1]);
  if (argc > 2) seed = atoi (argv[2]);
  if (argc > 3) TOLERANCE = atof (argv[3]); // Poisson tolerance, default 1e-3
  srand48 (seed);
  for (int n = 0; n < NMODES; n++)
    phase[n] = 2.*pi*drand48();

  L0 = 100.;
  periodic (right);
  rho1 = 1., rho2 = 1./830.;
  mu1 = 0.01, mu2 = 0.01/45.;
  f.sigma = 1.;
  init_grid (1 << MINLEVEL);
  run();
}

/**
True liquid volume (factor 2*pi included). */

double liquid_volume()
{
  double v = 0.;
  foreach (reduction(+:v))
    v += f[]*dv();
  return 2.*pi*v;
}

event init (t = 0)
{
  /**
  The maximum slope of the interface is about 0.25, so that a band of
  two (parent) cells on each side of the interface is enough to
  guarantee that all interfacial cells are at MAXLEVEL. */

  refine (level < MAXLEVEL && fabs(y - radius(x)) < 2.*Delta);
  fraction (f, radius(x) - y);
  vol0 = liquid_volume();
  fprintf (stderr, "# level %d seed %d Delta_min %g vol0 %.10g (pi*L0*(1+rms^2) = %.10g)\n",
	   MAXLEVEL, seed, L0/(1 << MAXLEVEL), vol0, pi*L0*(1. + sq(RMS)));
  fprintf (stderr, "# t i dt cells wall ndrops volume mgp.i mgu.i interfacial_cells interfacial_cells_below_maxlevel\n");
}

/**
With the default settings the wavelet error of f is computed with the
interface-reconstruction prolongation of the VOF tracer. This error is
small where the interface is nearly straight, so interfacial cells can
be coarsened below MAXLEVEL. Compiling with -DFULL_INTERFACE=1 adapts
on a copy of f which uses the default bilinear prolongation, which
keeps the whole interface at MAXLEVEL. */

event adapt (i++)
{
#if FULL_INTERFACE
  scalar ff[];
  foreach()
    ff[] = f[];
  adapt_wavelet ({ff,u}, (double[]){1e-3,1e-2,1e-2}, MAXLEVEL, MINLEVEL);
#else
  adapt_wavelet ({f,u}, (double[]){1e-3,1e-2,1e-2}, MAXLEVEL, MINLEVEL);
#endif
}

event stats (i++)
{
  if (grid->tn > cells_max) cells_max = grid->tn;
  if (i > 0 && dt < dtmin) dtmin = dt;
  if (mgp.i > mgp_max) mgp_max = mgp.i;
}

/**
Drops are counted every 0.1 time unit (to detect the first pinch-off).
The list of drops is written every time unit. Columns of the drops
files: index, true volume, axial centroid, D/W where D = (6V/pi)^(1/3).
The axial centroid is a circular mean because the domain is periodic
in x (tag() connects cells across the periodic boundary). */

event drops (t = 0; t <= TEND; t += 0.1)
{
  scalar m[];
  foreach()
    m[] = f[] > 1e-3;
  int n = tag (m);

  double v[n], c[n], s[n];
  for (int j = 0; j < n; j++)
    v[j] = c[j] = s[j] = 0.;
  foreach (serial)
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      c[j] += dv()*f[]*cos(2.*pi*x/L0);
      s[j] += dv()*f[]*sin(2.*pi*x/L0);
    }

  if (n > 1 && t_pinch < 0.)
    t_pinch = t;

  /**
  Number of interfacial cells, and how many of them are coarser than
  MAXLEVEL (this should be zero). */

  long ni = 0, nc = 0;
  foreach (reduction(+:ni) reduction(+:nc))
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      ni++;
      if (level < MAXLEVEL) nc++;
    }

  fprintf (stderr, "%g %d %g %ld %g %d %.10g %d %d %ld %ld\n",
	   t, i, dt, grid->tn, perf.t, n, liquid_volume(), mgp.i, mgu.i,
	   ni, nc);

  static int count = 0; // this event is called every 0.1 time unit
  if (count++ % 10 == 0) {
    char name[80];
    sprintf (name, "drops-%d-%d-%02d.dat", MAXLEVEL, seed, count/10);
    FILE * fp = fopen (name, "w");
    fprintf (fp, "# t = %g, %d drops\n# j volume xc D/W\n", t, n);
    for (int j = 0; j < n; j++) {
      double xc = atan2(s[j], c[j])*L0/(2.*pi);
      if (xc < 0.) xc += L0;
      fprintf (fp, "%d %.8g %.6g %.6g\n", j, 2.*pi*v[j], xc,
	       cbrt(12.*v[j])/W);
    }
    fclose (fp);
  }
}

event end (t = TEND)
{
  double vol = liquid_volume();
  fprintf (stderr, "# END level %d seed %d steps %d wall %g cells_mean %g "
	   "cells_max %ld updates_per_s %g dtmin %g mgp_max %d t_pinch %g "
	   "vol0 %.10g vol %.10g relerr %.3e\n",
	   MAXLEVEL, seed, i, perf.t, perf.tnc/(double)i, cells_max,
	   perf.speed, dtmin, mgp_max, t_pinch, vol0, vol, (vol - vol0)/vol0);
}
