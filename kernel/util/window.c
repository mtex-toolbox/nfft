/*
 * Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "api.h"

#if defined(DIRAC_DELTA)
  static const INT m2K_[] = {0};
#elif defined(GAUSSIAN)
  static const INT m2K_[] = {0, 1, 3, 6, 7, 9, 11, 13, 15, 17, 19, 21, 22, 23, 24};
#elif defined(B_SPLINE)
  static const INT m2K_[] = {0, 0, 4, 7, 10, 13, 15, 17, 19, 22, 24};
#elif defined(SINC_POWER)
  static const INT m2K_[] = {0, 0, 2, 5, 8, 11, 12, 14, 16, 18, 21, 23, 24, 24};
#else /* Kaiser-Bessel is the default. */
  static const INT m2K_[] = {1, 3, 7, 9, 14, 17, 20, 23, 24};
#endif

/**
 * Returns an appropriate value of the parameter K used with the PRE_LIN_PSI
 * flag for a given value of the cut-off parameter m.
 */
INT Y(m2K)(const INT m)
{
  int j = MIN(((int)(m)), ((int)((sizeof(m2K_) / sizeof(m2K_[0])) - 1)));
  return (INT)((1U << m2K_[j]) * (m + 2));
}

/**
 * Returns the default window cut off m for the selected window
 */
NFFT_INT Y(get_default_window_cut_off)()
{
  return (NFFT_INT)(WINDOW_HELP_ESTIMATE_m);
}

const char *Y(get_window_name)()
{
  return STRINGIZE(WINDOW_NAME);
}

/* dst[l] = sum_k c[k w + l] s^k, the values of the run side by side */
static inline void kb_horner(R *restrict dst, const R *restrict c, const INT D,
  const INT w, const R s)
{
  INT k, l;

  for (l = 0; l < w; l++)
    dst[l] = c[D * w + l];

  for (k = D - 1; k >= 0; k--)
    for (l = 0; l < w; l++)
      dst[l] = dst[l] * s + c[k * w + l];
}

void Y(kb_poly_run)(R *dst, const R *c, INT D, INT w, R s)
{
  /* a constant width for the usual cutoffs vectorises the run */
  switch (w)
  {
    case 6: kb_horner(dst, c, D, 6, s); break;
    case 8: kb_horner(dst, c, D, 8, s); break;
    case 10: kb_horner(dst, c, D, 10, s); break;
    case 12: kb_horner(dst, c, D, 12, s); break;
    case 14: kb_horner(dst, c, D, 14, s); break;
    case 16: kb_horner(dst, c, D, 16, s); break;
    case 18: kb_horner(dst, c, D, 18, s); break;
    default: kb_horner(dst, c, D, w, s);
  }
}

/* Interpolates the run at the D+1 Chebyshev points of the lowest degree D
 * whose polynomials match the exact run to 64 eps of its peak on a check grid.
 * The window is entire, so the degree stays small: 13 in double for m = 4..6. */
void Y(kb_poly_init)(R *poly, R b, R lg_tail, R peak_inv, INT m)
{
  const INT w = 2 * m + 2, nq = 4 * KB_POLY_DMAX;
  R *f = (R*) Y(malloc)((size_t)((KB_POLY_DMAX + 1) * w) * sizeof(R));
  R *T = (R*) Y(malloc)((size_t)((KB_POLY_DMAX + 1) * (KB_POLY_DMAX + 1)) * sizeof(R));
  R *ch = (R*) Y(malloc)((size_t)((KB_POLY_DMAX + 1) * w) * sizeof(R));
  R e[w], p[w];
  INT D, i, j, k, l, q;

  poly[0] = K(-1.0);

  for (D = 4; D <= KB_POLY_DMAX && poly[0] < K(0.0); D++)
  {
    const INT n = D + 1;
    R *c = poly + 1, err = K(0.0), peak = K(0.0);

    for (i = 0; i < n; i++)
      Y(kb_phi_run)(f + i * w, b, lg_tail, peak_inv, (R)m, m,
          (R)m + (COS(KPI * ((R)i + K(0.5)) / (R)n) + K(1.0)) / K(2.0));

    /* Chebyshev coefficients, then T[j k] the coefficient of s^k in T_j */
    for (j = 0; j < n; j++)
      for (l = 0; l < w; l++)
      {
        R s = K(0.0);
        for (i = 0; i < n; i++)
          s += f[i * w + l] * COS(KPI * (R)j * ((R)i + K(0.5)) / (R)n);
        ch[j * w + l] = (j ? K(2.0) : K(1.0)) * s / (R)n;
      }

    for (j = 0; j < n * n; j++)
      T[j] = K(0.0);
    T[0] = K(1.0);
    if (n > 1)
      T[n + 1] = K(1.0);
    for (j = 2; j < n; j++)
      for (k = 0; k < n; k++)
        T[j * n + k] = (k ? K(2.0) * T[(j - 1) * n + k - 1] : K(0.0)) - T[(j - 2) * n + k];

    for (k = 0; k < n; k++)
      for (l = 0; l < w; l++)
      {
        R s = K(0.0);
        for (j = k; j < n; j++)
          s += ch[j * w + l] * T[j * n + k];
        c[k * w + l] = s;
      }

    for (q = 0; q <= nq; q++)
    {
      const R t = (R)q / (R)nq;
      Y(kb_phi_run)(e, b, lg_tail, peak_inv, (R)m, m, (R)m + t);
      kb_horner(p, c, D, w, K(2.0) * t - K(1.0));
      for (l = 0; l < w; l++)
      {
        err = MAX(err, FABS(e[l] - p[l]));
        peak = MAX(peak, FABS(e[l]));
      }
    }

    if (err <= K(64.0) * EPSILON * peak)
      poly[0] = (R)D;
  }

  Y(free)(f);
  Y(free)(T);
  Y(free)(ch);
}
