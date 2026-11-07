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

/* Matlab gateway to the RNFFT, the NFFT of real-valued functions. The commands
 * are those of nfftmex; f is real, and f_hat holds the N_total coefficients of
 * the half spectrum in the order of the plan. */

#include "config.h"

#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include "nfft3.h"
#include "infft.h"
#include "imex.h"

#undef X
#define X(name) RNFFT(name)

#define PLANS_START 10 /* initial number of plans */
#define CMD_LEN_MAX 40 /* maximum length of command argument */

/* global flags */
#define RNFFT_MEX_FIRST_CALL (1U << 0)
static unsigned short gflags = RNFFT_MEX_FIRST_CALL;

static X(plan)** plans = NULL; /* plans */
static int plans_num_allocated = 0;
static char cmd[CMD_LEN_MAX];

#if defined(NFFT_SINGLE)
#define MX_CLASS mxSINGLE_CLASS
typedef float mxR;
#else
#define MX_CLASS mxDOUBLE_CLASS
typedef double mxR;
#endif

static inline void get_guru(const mxArray *prhs[], int d, int *N, int *M, int *n,
  int *m, unsigned int *f1, unsigned int *f2)
{
  int k;

  DM(if (mxGetNumberOfElements(prhs[1]) != (size_t)(2*d+5))
    mexErrMsgTxt("init_guru: Cell array must have 2*d+5 many entries: {d, N(1),...,N(d), M, n(1),...,n(d), m, nfft_flag, fftw_flags}");)

  for (k = 0; k < d; k++)
  {
    N[k] = mxGetScalar(mxGetCell(prhs[1], (unsigned int)(1+k)));
    DM(if (N[k] < 2 || (N[k]%2 != 0))
      mexErrMsgTxt("init_guru: input vector N must consist of even natural numbers");)
  }

  *M = mxGetScalar(mxGetCell(prhs[1], (unsigned int)(d+1)));
  DM(if (*M < 1)
    mexErrMsgTxt("init_guru: input argument M must be a natural number");)

  for (k = 0; k < d; k++)
    n[k] = mxGetScalar(mxGetCell(prhs[1], (unsigned int)(d+2+k)));

  *m = mxGetScalar(mxGetCell(prhs[1], (unsigned int)(2*d+2)));
  DM(if (*m < 1)
    mexErrMsgTxt("init_guru: input argument m must be a natural number");)

  *f1 = mxGetScalar(mxGetCell(prhs[1], (unsigned int)(2*d+3)));
  *f2 = mxGetScalar(mxGetCell(prhs[1], (unsigned int)(2*d+4)));
}

/* the plan named by the second argument */
static X(plan)* get_plan(const int nrhs, const mxArray *prhs[], const int nargs)
{
  int i;

  nfft_mex_check_nargs(nrhs, nargs, "rnfft: Wrong number of arguments.");
  i = nfft_mex_get_int(prhs[1], "rnfft: Input argument plan must be a scalar.");
  DM(if (i < 0 || i >= plans_num_allocated || plans[i] == 0)
    mexErrMsgTxt("Plan was not initialized or has already been finalized");)
  return plans[i];
}

static inline int mkplan(void)
{
  int i = 0;

  while (i < plans_num_allocated && plans[i] != 0)
    i++;
  if (i == plans_num_allocated)
  {
    int l;
    X(plan)** plans_old = plans;

    if (plans_num_allocated >= INT_MAX - PLANS_START - 1)
      mexErrMsgTxt("rnfft: Too many plans already allocated.");

    plans = Y(malloc)((plans_num_allocated+PLANS_START)*sizeof(X(plan)*));
    for (l = 0; l < plans_num_allocated; l++)
      plans[l] = plans_old[l];
    for (l = plans_num_allocated; l < plans_num_allocated+PLANS_START; l++)
      plans[l] = 0;
    if (plans_num_allocated > 0)
      Y(free)(plans_old);
    plans_num_allocated += PLANS_START;
  }
  plans[i] = Y(malloc)(sizeof(X(plan)));
  return i;
}

/* cleanup on mex function unload */
static void cleanup(void)
{
  int i;

  if (!(gflags & RNFFT_MEX_FIRST_CALL))
  {
    for (i = 0; i < plans_num_allocated; i++)
      if (plans[i])
      {
        X(finalize)(plans[i]);
        Y(free)(plans[i]);
        plans[i] = 0;
      }

    if (plans_num_allocated > 0)
    {
      Y(free)(plans);
      plans = NULL;
      plans_num_allocated = 0;
    }
    gflags |= RNFFT_MEX_FIRST_CALL;
  }
}

/* a real Matlab array, double or single, into the n entries of a */
static void get_real(const mxArray *p, R *a, const size_t n, const char *errmsg)
{
  size_t k;

  DM(if (mxGetNumberOfElements(p) != n || mxIsComplex(p))
    mexErrMsgTxt(errmsg);)
  if (mxIsDouble(p))
  {
    const double *v = mxGetPr(p);
    for (k = 0; k < n; k++)
      a[k] = v[k];
  }
  else if (mxIsSingle(p))
  {
    const float *v = (float*)mxGetData(p);
    for (k = 0; k < n; k++)
      a[k] = v[k];
  }
  else
    mexErrMsgTxt(errmsg);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (gflags & RNFFT_MEX_FIRST_CALL)
  {
    /* Force Matlab to load libfftw3. There is at least one version of Matlab
     * which otherwise crashes upon invocation of this mex function. */
    mexEvalString("fft([1,2,3,4]);");

    nfft_mex_install_mem_hooks();
    nfft_mex_lock_openmp_runtime();

    mexAtExit(cleanup);
    gflags &= ~RNFFT_MEX_FIRST_CALL;
  }

  /* command string */
  DM(if (nrhs == 0)
    mexErrMsgTxt("At least one input required.");)

  DM(if (!mxIsChar(prhs[0]))
    mexErrMsgTxt("First argument must be a string.");)

  if (mxGetString(prhs[0], cmd, CMD_LEN_MAX))
    mexErrMsgTxt("Could not get command string.");

  if (strcmp(cmd,"init_1d") == 0)
  {
    int i, N, M;

    nfft_mex_check_nargs(nrhs,3,"Wrong number of arguments for init_1d.");
    nfft_mex_get_nm(prhs,&N,&M);
    i = mkplan();
    X(init_1d)(plans[i],N,M);
    plhs[0] = mxCreateDoubleScalar((double)i);
  }
  else if (strcmp(cmd,"init_2d") == 0)
  {
    int i, N1, N2, M;

    nfft_mex_check_nargs(nrhs,4,"Wrong number of arguments for init_2d.");
    nfft_mex_get_n1n2m(prhs,&N1,&N2,&M);
    i = mkplan();
    X(init_2d)(plans[i],N1,N2,M);
    plhs[0] = mxCreateDoubleScalar((double)i);
  }
  else if (strcmp(cmd,"init_3d") == 0)
  {
    int i, N1, N2, N3, M;

    nfft_mex_check_nargs(nrhs,5,"Wrong number of arguments for init_3d.");
    nfft_mex_get_n1n2n3m(prhs,&N1,&N2,&N3,&M);
    i = mkplan();
    X(init_3d)(plans[i],N1,N2,N3,M);
    plhs[0] = mxCreateDoubleScalar((double)i);
  }
  else if (strcmp(cmd,"init") == 0)
  {
    int i, t, d, M;

    nfft_mex_check_nargs(nrhs,3,"Wrong number of arguments for init.");
    DM(if (!mxIsDouble(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) > 2)
      mexErrMsgTxt("init: input argument N must be a double array");)
    d = mxGetNumberOfElements(prhs[1]);
    if (d < 1)
      mexErrMsgTxt("init: input argument N must be a double array of length >= 1");
    {
      int N[d];
      const double *N_vec = mxGetPr(prhs[1]);

      for (t = 0; t < d; t++)
      {
        N[t] = (int)N_vec[t];
        DM(if (N[t] < 2 || (N[t]%2 != 0))
          mexErrMsgTxt("init: input vector N must consist of even natural numbers");)
      }
      M = nfft_mex_get_int(prhs[2],"init: input argument M must be a scalar.");
      DM(if (M < 1)
        mexErrMsgTxt("init: input argument M must be a natural number");)
      i = mkplan();
      X(init)(plans[i],d,N,M);
      plhs[0] = mxCreateDoubleScalar((double)i);
    }
  }
  else if (strcmp(cmd,"init_guru") == 0)
  {
    int d, i;

    nfft_mex_check_nargs(nrhs,2,"Wrong number of arguments for init_guru");
    DM(if (!mxIsCell(prhs[1]) || mxGetNumberOfElements(prhs[1]) < 3)
      mexErrMsgTxt("init_guru: second argument must be a cell array of length 2*d+5");)
    DM(if (!mxIsScalar(mxGetCell(prhs[1], 0)))
      mexErrMsgTxt("init_guru: cell array entry at index 0 (dimension d) must be scalar");)
    d = mxGetScalar(mxGetCell(prhs[1], 0));
    DM(if (d < 1)
      mexErrMsgTxt("init_guru: cell array entry at index 0 (dimension d) must be a natural number");)
    {
      int N[d], n[d], m, M;
      unsigned int f1, f2;

      get_guru(prhs,d,N,&M,n,&m,&f1,&f2);
      i = mkplan();
      X(init_guru)(plans[i],d,N,M,n,m,
        f1 | MALLOC_X | MALLOC_F | MALLOC_F_HAT | FFTW_INIT, f2);
      plhs[0] = mxCreateDoubleScalar((double)i);
    }
  }
  else if (strcmp(cmd,"precompute_psi") == 0)
    X(precompute_one_psi)(get_plan(nrhs,prhs,2));
  else if (strcmp(cmd,"trafo") == 0)
    X(trafo)(get_plan(nrhs,prhs,2));
  else if (strcmp(cmd,"adjoint") == 0)
    X(adjoint)(get_plan(nrhs,prhs,2));
  else if (strcmp(cmd,"trafo_direct") == 0)
    X(trafo_direct)(get_plan(nrhs,prhs,2));
  else if (strcmp(cmd,"adjoint_direct") == 0)
    X(adjoint_direct)(get_plan(nrhs,prhs,2));
  else if (strcmp(cmd,"finalize") == 0)
  {
    X(plan) *p = get_plan(nrhs,prhs,2);
    int i = nfft_mex_get_int(prhs[1],"rnfft: Input argument plan must be a scalar.");

    X(finalize)(p);
    Y(free)(p);
    plans[i] = 0;
  }
  else if (strcmp(cmd,"solver") == 0)
  {
    X(plan) *p = get_plan(nrhs,prhs,4);
    unsigned iterations = nfft_mex_get_int(prhs[2],"rnfft: Input argument iterations must be a scalar.");
    unsigned flags = nfft_mex_get_int(prhs[3],"rnfft: Input argument flags must be a scalar.");
    SOLVER(plan_c2r) ip;
    unsigned l;
    INT k;

    SOLVER(init_advanced_c2r)(&ip, (NFFT(mv_plan_c2r)*)p, flags ? flags : CGNR);
    for (k = 0; k < p->M_total; k++)
      ip.y[k] = p->f[k];
    /* the coefficients in the plan are the initial guess */
    for (k = 0; k < p->N_total; k++)
      ip.f_hat_iter[k] = p->f_hat[k];
    SOLVER(before_loop_c2r)(&ip);
    for (l = 0; l < iterations; l++)
      SOLVER(loop_one_step_c2r)(&ip);
    for (k = 0; k < p->N_total; k++)
      p->f_hat[k] = ip.f_hat_iter[k];
    SOLVER(finalize_c2r)(&ip);
  }
  else if (strcmp(cmd,"get_x") == 0)
  {
    X(plan) *p = get_plan(nrhs,prhs,2);
    const mwSize dims[2] = {p->d, p->M_total};
    INT k;
    mxR *x;

    plhs[0] = mxCreateNumericArray(2, dims, MX_CLASS, mxREAL);
    x = (mxR*)mxGetData(plhs[0]);
    for (k = 0; k < p->d * p->M_total; k++)
      x[k] = p->x[k];
  }
  else if (strcmp(cmd,"get_f") == 0)
  {
    X(plan) *p = get_plan(nrhs,prhs,2);
    const mwSize dims[2] = {p->M_total, 1};
    INT j;
    mxR *f;

    plhs[0] = mxCreateNumericArray(2, dims, MX_CLASS, mxREAL);
    f = (mxR*)mxGetData(plhs[0]);
    for (j = 0; j < p->M_total; j++)
      f[j] = p->f[j];
  }
  else if (strcmp(cmd,"get_f_hat") == 0)
  {
    X(plan) *p = get_plan(nrhs,prhs,2);
    const mwSize dims[2] = {p->N_total, 1};
    INT k;
    mxR *fr, *fi;

    plhs[0] = mxCreateNumericArray(2, dims, MX_CLASS, mxCOMPLEX);
    fr = (mxR*)mxGetData(plhs[0]);
    fi = (mxR*)mxGetImagData(plhs[0]);
    for (k = 0; k < p->N_total; k++)
    {
      fr[k] = CREAL(p->f_hat[k]);
      fi[k] = CIMAG(p->f_hat[k]);
    }
  }
  else if (strcmp(cmd,"set_x") == 0)
  {
    X(plan) *p = get_plan(nrhs,prhs,3);

    DM(if (mxGetNumberOfDimensions(prhs[2]) > 2 || mxGetM(prhs[2]) != (size_t)p->d
      || mxGetN(prhs[2]) != (size_t)p->M_total)
      mexErrMsgTxt("Input argument x must be a real d x M array.");)
    get_real(prhs[2], p->x, (size_t)(p->d * p->M_total),
      "Input argument x must be a real d x M array.");
  }
  else if (strcmp(cmd,"set_f") == 0)
  {
    X(plan) *p = get_plan(nrhs,prhs,3);

    get_real(prhs[2], p->f, (size_t)p->M_total,
      "Input argument f must be a real array of M entries.");
  }
  else if (strcmp(cmd,"set_f_hat") == 0)
  {
    X(plan) *p = get_plan(nrhs,prhs,3);
    const INT n = p->N_total;
    INT k;

    DM(if (mxGetM(prhs[2]) != (size_t)n || mxGetN(prhs[2]) != 1)
      mexErrMsgTxt("Input argument f_hat must be a column vector of N_total entries.");)
    if (mxIsDouble(prhs[2]))
    {
      const double *fr = mxGetPr(prhs[2]), *fi = mxGetPi(prhs[2]);
      for (k = 0; k < n; k++)
        p->f_hat[k] = fi ? fr[k] + II * fi[k] : fr[k];
    }
    else if (mxIsSingle(prhs[2]))
    {
      const float *fr = (float*)mxGetData(prhs[2]), *fi = (float*)mxGetImagData(prhs[2]);
      for (k = 0; k < n; k++)
        p->f_hat[k] = fi ? fr[k] + II * fi[k] : fr[k];
    }
    else
      mexErrMsgTxt("Input argument f_hat must be a double array");
  }
  else if (strcmp(cmd,"display") == 0)
  {
    X(plan) *p = get_plan(nrhs,prhs,2);

    mexPrintf("   pointer: %p\n",(void*)p);
    mexPrintf("         d: %d\n",(int)p->d);
    mexPrintf("   N_total: %d\n",(int)p->N_total);
    mexPrintf("   M_total: %d\n",(int)p->M_total);
    mexPrintf("         m: %d\n",(int)p->m);
    mexPrintf("     flags: %u\n",p->flags);
    mexPrintf("fftw_flags: %u\n",p->fftw_flags);
  }
  else if (strcmp(cmd,"get_num_threads") == 0)
    plhs[0] = mxCreateDoubleScalar((double)Y(get_num_threads)());
  else if (strcmp(cmd,"set_num_threads") == 0)
    Y(set_num_threads)(nfft_mex_set_num_threads_check(nrhs, prhs, (void **)plans,
      plans_num_allocated));
  else if (strcmp(cmd,"has_threads_enabled") == 0)
    plhs[0] = mxCreateDoubleScalar((double)Y(has_threads_enabled)());
  else if (strcmp(cmd,"get_default_window_cut_off_m") == 0)
    plhs[0] = mxCreateDoubleScalar((double)WINDOW_HELP_ESTIMATE_m);
  else if (strcmp(cmd,"get_epsilon") == 0)
    plhs[0] = mxCreateDoubleScalar((double)Y(float_property)(NFFT_EPSILON));
  else
    mexErrMsgTxt("rnfft: Unknown command.\n");
}
