#include "fftw3.h"
#include "stdlib.h"

#ifndef _LTFAT_MEX_FILE
#define _LTFAT_MEX_FILE

#define ISNARGINEQ 2
#define TYPEDEPARGS 0
#define SINGLEARGS

static fftw_plan* p_old = 0;

static void ifftrealAtExit(void)
{
  if(p_old!=0)
  {
     fftw_destroy_plan(*p_old);
     free(p_old);
  }
}

#endif // _LTFAT_MEX_FILE - INCLUDED ONCE

#define MEX_FILE __BASE_FILE__
#include "ltfat_mex_template_helper.h"

#if defined(LTFAT_SINGLE) || defined(LTFAT_DOUBLE)
#include "ltfat_types.h"
#include "config.h"
#include "fftw3.h"

// Calling convention:
//  comp_ifftreal(f,N);

void LTFAT_NAME(ltfatMexFnc)( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  #ifdef LTFAT_DOUBLE
  if(p_old==0)
  {
      mexAtExit(ifftrealAtExit);
  }
  #endif
  mwSize ii, L, W, L2;
  LTFAT_FFTW(plan) p;
  LTFAT_REAL *f, s;
  LTFAT_REAL *fin_r, *fin_i;


  L2 = (mwSize) mxGetM(prhs[0]);
  W  = (mwSize) mxGetN(prhs[0]);
  L  = (mwSize) mxGetScalar(prhs[1]);


  if(L/2+1!=L2)
    mexErrMsgTxt("Invalid output length");

  fin_r = (LTFAT_REAL*)mxGetPr(prhs[0]);
  fin_i = (LTFAT_REAL*)mxGetPi(prhs[0]);
  // Case when input is real
  if(!mxIsComplex(prhs[0]))
  {
    mxArray* tmpIn = ltfatCreateMatrix(L2, W, LTFAT_MX_CLASSID , mxCOMPLEX);
    LTFAT_REAL *fin_r_old = (LTFAT_REAL*)mxGetPr(prhs[0]);
    fin_r = (LTFAT_REAL*)mxGetPr(tmpIn);
    fin_i = (LTFAT_REAL*)mxGetPi(tmpIn);
    for(mwIndex jj=0;jj<L2*W;jj++)
    {
        fin_r[jj]= fin_r_old[jj];
        fin_i[jj]= (LTFAT_REAL )0.0;
    }
  }

  // Create output and get pointer
  plhs[0] = ltfatCreateMatrix(L, W, LTFAT_MX_CLASSID , mxREAL);
  f= (LTFAT_REAL*) mxGetPr(plhs[0]);

  // This section is not being compiled. It contains a segmentation
  // faults. The idea is to pass Matlab's split memory layout directly
  // to FFTW

  LTFAT_FFTW(iodim) dims[1], howmanydims[1];

  // Create plan. Copy data from cin to f.
  dims[0].n = L;
  dims[0].is = 1;
  dims[0].os = 1;

  howmanydims[0].n = W;
  howmanydims[0].is = L2;
  howmanydims[0].os = L;

  // The calling prototype

  // fftw_plan fftw_plan_guru_split_dft_c2r(
  //        int rank, const fftw_iodim *dims,
  //        int howmany_rank, const fftw_iodim *howmany_dims,
  //        double *ri, double *ii, double *out,
  //        unsigned flags);


  p = LTFAT_FFTW(plan_guru_split_dft_c2r)(1, dims,
				   1, howmanydims,
				   fin_r,fin_i,f, FFTW_OPTITYPE);

  if(p_old!=0)
  {
    fftw_destroy_plan(*p_old);
    free(p_old);
  }
  p_old = malloc(sizeof(p));
  memcpy(p_old,&p,sizeof(p));


  // Real IFFT.
  LTFAT_FFTW(execute)(p);

  //LTFAT_FFTW(destroy_plan)(p);

  // Scale, because FFTW's normalization is different.
  s  = (LTFAT_REAL) (1.0/((LTFAT_REAL)L));
  for (ii=0; ii<L*W; ii++)
    {
      f[ii] *=s;
    }

  return;
}
#endif


/*
#include "mex.h"
#include "config.h"
#include "fftw3.h"
#include "ltfat-mex-helper.h"

// Calling convention:
//  comp_ifftreal(f,N);

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )
{
  int ii, L, W, L2;
  fftw_plan p;
  double *f, s;

#ifndef FIXME_THIS_SECTION_CONTAINS_AN_ERROR

  L2 = mxGetM(prhs[0]);
  W  = mxGetN(prhs[0]);
  L  = (int)mxGetScalar(prhs[1]);

  // Create output and get pointer
  plhs[0] = mxCreateDoubleMatrix(L, W, mxREAL);
  f=mxGetPr(plhs[0]);

  // This section is not being compiled. It contains a segmentation
  // faults. The idea is to pass Matlab's split memory layout directly
  // to FFTW

  fftw_iodim dims[1], howmanydims[1];

  // Create plan. Copy data from cin to f.
  dims[0].n = L;
  dims[0].is = 1;c2r
  dims[0].os = 1;

  howmanydims[0].n = W;
  howmanydims[0].is = L2;
  howmanydims[0].os = L;

  // The calling prototype

  // fftw_plan fftw_plan_guru_split_dft_c2r(
  //        int rank, const fftw_iodim *dims,
  //        int howmany_rank, const fftw_iodim *howmany_dims,
  //        double *ri, double *ii, double *out,
  //        unsigned flags);


  p = fftw_plan_guru_split_dft_c2r(1, dims,
				   1, howmanydims,
				   mxGetPr(prhs[0]), mxGetPi(prhs[0]),
				   f, FFTW_ESTIMATE);

  // Real IFFT.
  fftw_execute(p);

#else

  // This section copies to a combined layout, and use the exact same
  // FFTW call as the Octave interface.

  fftw_complex *c_combined;
c2r
  L2 = mxGetM(prhs[0]);
  W  = mxGetN(prhs[0]);
  L  = (int)mxGetScalar(prhs[1]);

  // Create output and get pointer
  plhs[0] = mxCreateDoubleMatrix(L, W, mxREAL);
  f=mxGetPr(plhs[0]);

  c_combined = mxMalloc(L2*W*sizeof(fftw_complex));

  split2combined(L2*W, prhs[0], c_combined);

  // Create plan. Copy data from f to cout.
  p = fftw_plan_many_dft_c2r(1, &L, W,
			     c_combined, NULL,
			     1, L2,
			     f, NULL,
			     1, L,
			     FFTW_ESTIMATE);

  // Real IFFT.
  fftw_execute(p);

  mxFree(c_combined);

#endif

  fftw_destroy_plan(p);

  // Scale, because FFTW's normalization is different.
  s  = 1.0/L;
  for (ii=0; ii<L*W; ii++)
    {
      f[ii] *=s;
  }


  return;

}
*/
