
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <climits>

#include "meep.hpp"
#include "meep_internals.hpp"

/* This file contains a generic function for looping over all of the
   points in all of the chunks that intersect some given grid_volume.  This
   is used for everything from HDF5 output to applying source volumes to
   integrating energy and flux.  It's fairly tricky because of the
   parallelization, arbitrary chunk divisions, symmetries, and periodic
   boundary conditions, but at least all of the trickiness is in one
   place.  It is designed so that the inner loops over the actual grid
   points can be tight and fast (using the LOOP_OVER_IVECS macro).

   Many of the loops over chunks involve some sort of integration-like
   computation, and so we also perform the additional task of calculating
   the integration weights for each point -- mainly, this involves weighting
   the boundary points appropriately so that the sum approximates (via
   linear interpolation) a continuous integral over the supplied grid_volume. */

/****************************************************************************

                     Integration Weights

We want the integral from a to b, assuming linear interpolation of fn
(function values on grid points n).  Most interior points have weight
1, but the points just inside and outside the boundaries have
different weights.  Call the weights for the points just *outside* the
starting and ending boundaries s0 and e0, respectively, and weights
for the points just *inside* the boundaries s1 and e1.  Then we have
to handle the following cases:

1) a and b separated by at least 2 grid points, e.g.:

x  |    x       x       x  |    x
0  a    1       2       3  b    4

first segment: f(x) = f0 (1 - x) + f1 x
      -- \int_a^1 f(x) dx = f0 (1 - a)^2/2 + f1 (1 - a^2) / 2

last segment: f(x) = f3 (4 - x) + f4 (x - 3)
     -- \int_3^b f(x) dx = f3 [1 - (4-b)^2] / 2 + f4 (b - 3)^2 / 2

integral = f0 (1 - a)^2/2         <---- f0 s0
         + f1 (1 - a^2/2)         <---- f1 s1
         + f2
         + f3 (1 - (4-b)^2 / 2)   <---- f3 e1
         + f4 (b - 3)^2 / 2       <---- f4 e0

In terms of starting and ending weights:
         w0 = 1 - a
         w1 = b - 3

         s0 = w0^2 / 2
         s1 = 1 - (1 - w0)^2 / 2
         e0 = w1^2 / 2
         e1 = 1 - (1 - w1)^2 / 2

2) one grid point between a and b.

x  |    x  |    x
0  a    1  b    2

integral = f0 (1 - a)^2 / 2
         + f1 [(1 - a^2) + 1 - (2 - b)^2] / 2
         + f3 (b - 1)^2 / 2

s0 = w0^2 / 2
e0 = w1^2 / 2
s1 = e1 = 1 - (1 - w0)^2 / 2 - (1 - w1)^2 / 2

3) no grid points between a and b.

x  |      |    x
0  a      b    1

integral = f0 [ (1-a)^2 - (1-b)^2 ] / 2 + f1 [ b^2 - a^2 ] / 2
         = f0 [ w0^2 - (1-w1)^2 ]  / 2 + f1 [ w1^2 - (1-w0)^2 ] / 2

s0 = e1 = w0^2/2 - (1-w1)^2/2
e0 = s1 = w1^2/2 - (1-w0)^2/2

4) as (3), but a = b: interpolation, not integration:

 -- want: f0 * w0 + f1 * w1

s0 = w0
e0 = w1 = 1 - w0

                        --------------

          Integration Weights in Cylindrical Coordinates
                   FIXME: implement this below?

Ideally, we should have different weights for the R direction of
cylindrical coordinates, i.e. for integrating f(r) r dr, because again
we want to perfectly integrate any linear f(r).  Thus, the integration
weights will depend upon r.  Note, however, that we also have an r in
the dV, so we will have to divide the weights by this factor.

1) a and b separated by at least 2 grid points, e.g.:

x  |    x       x       x  |    x
i  a   i+1     i+2     i+3  b  i+4

(where r = i * inva).

linear interpolation in [i,i+1): f(x) = f_i (i+1 - x) + f_{i+1} (x-i)

   want: \int_a^b f(x) x dx

in terms of starting and ending weights:
         w0 = (i+1) - a
         w1 = b - (i+3)

integral = f_i [-w0^3 / 3 + (i+1) w0^2 / 2]                       <- s0 i
      + f_{j=i+1} [w0^3 / 3 - (j+1) w0^2 / 2 + j w0 + j/2 + 1/6]  <- s1 (i+1)
      + f_{j=i+2} j                                               <- 1 (i+2)
      + f_{j=i+3} [-w1^3 / 3 - (j-1) w1^2 / 2 + j w1 + j/2 - 1/6] <- e1 (i+3)
      + f_{j=i+4} [w1^3 / 3 + (j-1) w1^2 / 2]                     <- e0 (i+4)

      (thanks to Maple for doing the annoying algebra)
  (yes, I have tested that it correctly integrates linear f(r))

Note that the coefficients need to be divided by i, i+1, etcetera to
get s0, s1, etcetera; this gives an interior-point weight of 1 as
before.  For i->infinity, this should converge to the weights from
before.  Avoiding division by zero is more tricky, because the weight
at j=0 is not necessarily zero, due to the interpolation.  It might be
better to pre-include the dV in the weight for edge elements, with
appropriate logic in the IVEC_LOOP_WEIGHT macro.  Tricky.

The above is also not correct for integrals that cross x=0, because
it should really be the integral of f(x) |x|.  Even interior points
probably need special handling in that case.  For sanity, we would
just divide the integration region into positive and negative r and
integrate them separately somehow.  Grrr.

2) one grid point between a and b.

x  |    x   |    x
i  a   i+1  b   i+2

integral = f_i [-w0^3 / 3 + (i+1) w0^2 / 2]                 <- s0 i
      + f_{j=i+1} [w0^3 / 3 - (j+1) w0^2 / 2 + j w0 +
                  -w1^3 / 3 - (j-1) w1^2 / 2 + j w1]        <- {s1,e1} (i+1)
      + f_{j=i+2} [w1^3 / 3 + (j-1) w1^2 / 2]               <- e0 (i+2)

3) no grid points between a and b.

x  |      |    x
i  a      b   i+1

integral = f_i [-w0^3/3 + (i+1) w0^2/2
              + -w1^3/3 - (i-1) w1^2/2 + i w1 - i/2 - 1/6]  <- s0 i
   + f_{j=i+1} [ w0^3/3 - (j+1) w0^2/2 + j w0
              +  w1^3/3 + (j-1) w1^2/2 - j/2 + 1/6]         <- e0 (i+1)

4) as (3), but a = b: interpolation, not integration: same as above

****************************************************************************/

using namespace std;

namespace meep {

/***************************************************************/
/* get_field_components is a utility routine, designed to be   */
/* called by chunkloop functions, for fetching values of field */
/* components at grid points, accounting for the complications */
/* of symmetry and yee-grid averaging.                         */
/***************************************************************/
chunkloop_field_components::chunkloop_field_components(fields_chunk *fc, component cgrid,
                                                       std::complex<double> shift_phase,
                                                       const symmetry &S, int sn, int num_fields,
                                                       const component *components)
    : fc(fc), parent_components(num_fields), phases(num_fields), offsets(2 * num_fields),
      values(num_fields) {
  // for each requested component, get symmetry-parent component, yee-grid offsets, and phase shift
  for (int nc = 0; nc < num_fields; nc++) {
    parent_components[nc] = S.transform(components[nc], -sn);
    phases[nc] = shift_phase * S.phase_shift(parent_components[nc], sn);
    ptrdiff_t ofs1 = 0, ofs2 = 0;
    if (cgrid == Centered) fc->gv.yee2cent_offsets(parent_components[nc], ofs1, ofs2);
    offsets[2 * nc] = ofs1;
    offsets[2 * nc + 1] = ofs2;
  }
}

void chunkloop_field_components::update_values(ptrdiff_t idx) {
  for (size_t nc = 0; nc < values.size(); nc++) {
    // do appropriate averaging to get value of field component at grid point
    component cparent = parent_components[nc];
    ptrdiff_t ofs1 = offsets[2 * nc], ofs2 = offsets[2 * nc + 1];
    double favg[2] = {0.0, 0.0}; // real, imag parts
    for (int reim = 0; reim < 2; reim++) {
      const realnum *fgrid = fc->f[cparent][reim];
      if (!fgrid) continue;
      favg[reim] =
          0.25 * (fgrid[idx] + fgrid[idx + ofs1] + fgrid[idx + ofs2] + fgrid[idx + ofs1 + ofs2]);
    }
    values[nc] = phases[nc] * std::complex<double>(favg[0], favg[1]);
  }
}

/* The following two functions convert a vec to the nearest ivec
   in the dielectric (odd-coordinate) grid, either rounding down (floor)
   or up (ceil).  In the special case where a component of the vec is
   *exactly* on a component of the ivec, we add the corresponding
   component of equal_shift (which should be either -2, 0, or +2).
   (equal_shift is there to prevent us from counting edge points twice.) */

ivec fields::vec2diel_floor(const vec &pt, double a, const ivec &equal_shift) {
  ivec ipt(pt.dim);
  LOOP_OVER_DIRECTIONS(pt.dim, d) {
    ipt.set_direction(d, 1 + 2 * int(floor(pt.in_direction(d) * a - .5)));
    if (ipt.in_direction(d) == pt.in_direction(d))
      ipt.set_direction(d, ipt.in_direction(d) + equal_shift.in_direction(d));
  }
  return ipt;
}
ivec fields::vec2diel_ceil(const vec &pt, double a, const ivec &equal_shift) {
  ivec ipt(pt.dim);
  LOOP_OVER_DIRECTIONS(pt.dim, d) {
    ipt.set_direction(d, 1 + 2 * int(ceil(pt.in_direction(d) * a - .5)));
    if (ipt.in_direction(d) == pt.in_direction(d))
      ipt.set_direction(d, ipt.in_direction(d) + equal_shift.in_direction(d));
  }
  return ipt;
}

static inline int iabs(int i) { return (i < 0 ? -i : i); }

/* Integration weights at boundaries (c.f. long comment at top).   */
/* This code was formerly part of loop_in_chunks, now refactored   */
/* as a separate routine so we can call it from get_array_metadata.*/
void compute_boundary_weights(grid_volume gv, const volume &where, ivec &is, ivec &ie,
                              bool snap_empty_dimensions, vec &s0, vec &e0, vec &s1, vec &e1) {
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    double w0, w1;
    w0 = 1. - where.in_direction_min(d) * gv.a + 0.5 * is.in_direction(d);
    w1 = 1. + where.in_direction_max(d) * gv.a - 0.5 * ie.in_direction(d);
    if (ie.in_direction(d) >= is.in_direction(d) + 3 * 2) {
      s0.set_direction(d, w0 * w0 / 2);
      s1.set_direction(d, 1 - (1 - w0) * (1 - w0) / 2);
      e0.set_direction(d, w1 * w1 / 2);
      e1.set_direction(d, 1 - (1 - w1) * (1 - w1) / 2);
    }
    else if (ie.in_direction(d) == is.in_direction(d) + 2 * 2) {
      s0.set_direction(d, w0 * w0 / 2);
      s1.set_direction(d, 1 - (1 - w0) * (1 - w0) / 2 - (1 - w1) * (1 - w1) / 2);
      e0.set_direction(d, w1 * w1 / 2);
      e1.set_direction(d, s1.in_direction(d));
    }
    else if (where.in_direction_min(d) == where.in_direction_max(d)) {
      if (snap_empty_dimensions) {
        if (w0 > w1)
          ie.set_direction(d, is.in_direction(d));
        else
          is.set_direction(d, ie.in_direction(d));
        w0 = w1 = 1.0;
      }
      s0.set_direction(d, w0);
      s1.set_direction(d, w1);
      e0.set_direction(d, w1);
      e1.set_direction(d, w0);
    }
    else if (ie.in_direction(d) == is.in_direction(d) + 1 * 2) {
      s0.set_direction(d, w0 * w0 / 2 - (1 - w1) * (1 - w1) / 2);
      e0.set_direction(d, w1 * w1 / 2 - (1 - w0) * (1 - w0) / 2);
      s1.set_direction(d, e0.in_direction(d));
      e1.set_direction(d, s0.in_direction(d));
    }
    else
      meep::abort("bug: impossible(?) looping boundaries");
  }
}

/* Generic function for computing loops within the chunks, often
   integral-like things, over a grid_volume WHERE.  The job of this
   function is to call CHUNKLOOP() for each chunk that intersects
   WHERE, passing it the chunk, the range of integer coordinates to
   loop over, the integration weights for the boundary points, and the
   bloch phase shift, translational shift, and symmetry operation to
   transform the chunk to the actual integration location.  (N.B.
   we apply the symmetry first to the chunk, *then* the shift.)

   We also pass CHUNKLOOP() dV0 and dV1, such that the integration
   "grid_volume" dV is dV0 + dV1 * iloopR, where iloopR is the loop
   variable (starting from 0 at the starting integer coord and
   incrementing by 1) corresponding to the direction R.  Note that, in
   the LOOP_OVER_IVECS macro, iloopR corresponds to the loop variable
   loop_i2 in Dcyl (cylindrical coordinates).  In other coordinates,
   dV1 is 0.  Note also that by "grid_volume" dV we mean the integration
   unit corresponding to the dimensionality of WHERE (e.g. an area if
   WHERE is 2d, etc.)

   In particular, the loop's point coordinates are calculated on the
   Yee grid for component cgrid.  cgrid == Centered is a good choice
   if you want to work with a combination of multiple field
   components, because all of the field components can be interpolated
   onto this grid without communication between chunks.

   The integration weights are chosen to correspond to integrating the
   linear interpolation of the function values from these grid points.

   For a simple example of an chunkloop routine, see the
   tests/integrate.cpp file.

   The parameters USE_SYMMETRY (default = true) and SNAP_EMPTY_DIMENSIONS
   (default = false) are for use with not-quite-integration-like
   operations.  If use_symmetry is false, then we do *not* loop over
   all possible symmetry transformations of the chunks to see if they
   intersect WHERE; we only use chunks that, untransformed, already
   intersect the grid_volume.  If SNAP_EMPTY_DIMENSIONS is true, then for empty
   (min = max) dimensions of WHERE, instead of interpolating, we
   "snap" them to the nearest grid point.  */

void fields::loop_in_chunks(field_chunkloop chunkloop, void *chunkloop_data, const volume &where,
                            component cgrid, bool use_symmetry, bool snap_empty_dimensions) {
  if (coordinate_mismatch(gv.dim, cgrid))
    meep::abort("Invalid fields::loop_in_chunks grid type %s for dimensions %s\n",
                component_name(cgrid), dimension_name(gv.dim));
  if (where.dim != gv.dim)
    meep::abort("Invalid dimensions %d for WHERE in fields::loop_in_chunks", where.dim);

  if (cgrid == Permeability) cgrid = Centered;

  /* Find the corners (is and ie) of the smallest bounding box for
     wherec, on the grid of odd-coordinate ivecs (i.e. the
     "dielectric/centered grid") and then shift back to the yee grid for c. */
  vec yee_c(gv.yee_shift(Centered) - gv.yee_shift(cgrid));
  ivec iyee_c(gv.iyee_shift(Centered) - gv.iyee_shift(cgrid));
  volume wherec(where + yee_c);
  ivec is(vec2diel_floor(wherec.get_min_corner(), gv.a, zero_ivec(gv.dim)) - iyee_c);
  ivec ie(vec2diel_ceil(wherec.get_max_corner(), gv.a, zero_ivec(gv.dim)) - iyee_c);

  vec s0(gv.dim), e0(gv.dim), s1(gv.dim), e1(gv.dim);
  compute_boundary_weights(gv, where, is, ie, snap_empty_dimensions, s0, e0, s1, e1);

  int original_vol = 1;
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    if (boundaries[High][d] == Periodic && boundaries[Low][d] == Periodic) {
      original_vol *= (ie.in_direction(d) - is.in_direction(d)) / 2 + 1;
    }
    else {
      int vol_size_d = (min(user_volume.big_owned_corner(cgrid), ie).in_direction(d) -
                        max(user_volume.little_owned_corner(cgrid), is).in_direction(d)) /
                           2 +
                       1;
      original_vol *= (vol_size_d > 0 ? vol_size_d : 0);
    }
  }

  // loop over symmetry transformations of the chunks:
  int vol_sum = 0;
  for (int sn = 0; sn < (use_symmetry ? S.multiplicity() : 1); ++sn) {
    component cS = S.transform(cgrid, -sn);
    volume gvS = S.transform(gv.surroundings(), sn);
    vec L(gv.dim);
    ivec iL(gv.dim);

    // n.b. we can't just S.transform(lattice_vector,sn), 'cause of signs
    LOOP_OVER_DIRECTIONS(gv.dim, d) {
      direction dS = S.transform(d, -sn).d;
      L.set_direction(d, fabs(lattice_vector(dS).in_direction(dS)));
      iL.set_direction(d, iabs(ilattice_vector(dS).in_direction(dS)));
    }

    // figure out range of lattice shifts for which gvS intersects where:
    ivec min_ishift(gv.dim), max_ishift(gv.dim);
    LOOP_OVER_DIRECTIONS(gv.dim, d) {
      if (boundaries[High][S.transform(d, -sn).d] == Periodic) {
        min_ishift.set_direction(
            d,
            int(floor((where.in_direction_min(d) - gvS.in_direction_max(d)) / L.in_direction(d))));
        max_ishift.set_direction(d, int(ceil((where.in_direction_max(d) - gvS.in_direction_min(d)) /
                                             L.in_direction(d))));
      }
      else {
        min_ishift.set_direction(d, 0);
        max_ishift.set_direction(d, 0);
      }
    }

    // loop over lattice shifts
    ivec ishift(min_ishift);
    do {
      complex<double> ph = 1.0;
      vec shift(gv.dim, 0.0);
      ivec shifti(gv.dim, 0);
      LOOP_OVER_DIRECTIONS(gv.dim, d) {
        shift.set_direction(d, L.in_direction(d) * ishift.in_direction(d));
        shifti.set_direction(d, iL.in_direction(d) * ishift.in_direction(d));
        ph *= pow(eikna[d], ishift.in_direction(d));
      }

      for (int i = 0; i < num_chunks; ++i) {
        if (!chunks[i]->is_mine()) continue;
        grid_volume gvu(chunks[i]->gv);
        ivec _iscoS(S.transform(gvu.little_owned_corner(cS), sn));
        ivec _iecoS(S.transform(gvu.big_owned_corner(cS), sn));
        ivec iscoS(max(user_volume.little_owned_corner(cgrid), min(_iscoS, _iecoS))),
            iecoS(max(_iscoS, _iecoS)); // fix ordering due to to transform

        // With symmetry, the upper half of the original chunk is kept and includes one extra pixel.
        // When looped over all symmetries, pixels outside the lower boundary
        // "user_volume.little_owned_corner(cgrid)" is excluded. isym finds the lower boundary of
        // the upper half chunk. Then in each direction that chunk isn't the upper, checked by
        // ((S.transform(d, sn).d != d) != (S.transform(d, sn).flipped)), the end point iecoS will
        // shift to not go beyond isym
        ivec isym(gvu.dim, INT_MAX);
        LOOP_OVER_DIRECTIONS(gvu.dim, d) {
          int off_sym_shift = ((gv.iyee_shift(cgrid).in_direction(d) != 0)
                                   ? gv.iyee_shift(cgrid).in_direction(d) + 2
                                   : 2);
          isym.set_direction(d, S.i_symmetry_point.in_direction(d) - off_sym_shift);
        }

        ivec iscS(max(is - shifti, iscoS));
        ivec chunk_corner(gvu.little_owned_corner(cgrid));
        LOOP_OVER_DIRECTIONS(gv.dim, d) {
          if ((S.transform(d, sn).d != d) != (S.transform(d, sn).flipped))
            iecoS.set_direction(d, min(isym, iecoS).in_direction(d));
        }
        ivec iecS(min(ie - shifti, iecoS));

        if (iscS <= iecS) { // non-empty intersection
          // Determine weights at chunk looping boundaries:
          ivec isc(S.transform(iscS, -sn)), iec(S.transform(iecS, -sn));
          vec s0c(gv.dim, 1.0), s1c(gv.dim, 1.0), e0c(gv.dim, 1.0), e1c(gv.dim, 1.0);
          iscS += shifti;
          iecS += shifti;
          LOOP_OVER_DIRECTIONS(gv.dim, d) {
            direction dS = S.transform(d, sn).d;
            if (iscS.in_direction(dS) == is.in_direction(dS)) {
              s0c.set_direction(d, s0.in_direction(dS));
              s1c.set_direction(d, s1.in_direction(dS));
            }
            else if (iscS.in_direction(dS) == is.in_direction(dS) + 2) {
              s0c.set_direction(d, s1.in_direction(dS));
            }
            if (iecS.in_direction(dS) == ie.in_direction(dS)) {
              e0c.set_direction(d, e0.in_direction(dS));
              e1c.set_direction(d, e1.in_direction(dS));
            }
            else if (iecS.in_direction(dS) == ie.in_direction(dS) - 2) {
              e0c.set_direction(d, e1.in_direction(dS));
            }
            if (iecS.in_direction(dS) == iscS.in_direction(dS)) {
              double w = std::min(s0c.in_direction(d), e0c.in_direction(d));
              s0c.set_direction(d, w);
              e0c.set_direction(d, w);
              s1c.set_direction(d, w);
              e1c.set_direction(d, w);
            }
            else if (iecS.in_direction(dS) == iscS.in_direction(dS) + 1 * 2) {
              double w = std::min(s0c.in_direction(d), e1c.in_direction(d));
              s0c.set_direction(d, w);
              e1c.set_direction(d, w);
              w = std::min(s1c.in_direction(d), e0c.in_direction(d));
              s1c.set_direction(d, w);
              e0c.set_direction(d, w);
            }
            else if (iecS.in_direction(dS) == iscS.in_direction(dS) + 2 * 2) {
              double w = std::min(s1c.in_direction(d), e1c.in_direction(d));
              s1c.set_direction(d, w);
              e1c.set_direction(d, w);
            }

            // swap endpoints/weights if in wrong order due to S.transform
            if (isc.in_direction(d) > iec.in_direction(d)) {
              int iswap = isc.in_direction(d);
              isc.set_direction(d, iec.in_direction(d));
              iec.set_direction(d, iswap);
              double swap = s0c.in_direction(d);
              s0c.set_direction(d, e0c.in_direction(d));
              e0c.set_direction(d, swap);
              swap = s1c.in_direction(d);
              s1c.set_direction(d, e1c.in_direction(d));
              e1c.set_direction(d, swap);
            }
          }

          // Determine integration "volumes" dV0 and dV1;
          double dV0 = 1.0, dV1 = 0.0;
          LOOP_OVER_DIRECTIONS(gv.dim, d) {
            if (where.in_direction(d) > 0.0) dV0 *= gv.inva;
          }
          if (gv.dim == Dcyl) {
            dV1 = dV0 * 2 * pi * gv.inva;
            dV0 *= 2 * pi * fabs((S.transform(chunks[i]->gv[isc], sn) + shift).in_direction(R));
          }

          chunkloop(chunks[i], i, cS, isc, iec, s0c, s1c, e0c, e1c, dV0, dV1, shifti, ph, S, sn,
                    chunkloop_data);
          int loop_vol = 1;
          LOOP_OVER_DIRECTIONS(gv.dim, d) {
            loop_vol *= (iec.in_direction(d) - isc.in_direction(d)) / 2 + 1;
          }
          vol_sum += loop_vol;
        }
      }

      LOOP_OVER_DIRECTIONS(gv.dim, d) {
        if (ishift.in_direction(d) + 1 <= max_ishift.in_direction(d)) {
          ishift.set_direction(d, ishift.in_direction(d) + 1);
          break;
        }
        ishift.set_direction(d, min_ishift.in_direction(d));
      }
    } while (ishift != min_ishift);
  }
  int vol_sum_all = sum_to_all(vol_sum);
  if (use_symmetry && vol_sum_all != original_vol)
    master_printf("WARNING vol mismatch:, original_vol %i, looped vol_sum %i \n", original_vol,
                  vol_sum_all);
}

} // namespace meep
