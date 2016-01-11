#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define REALARGS
#define OCTFILENAME comp_dgtreal_long // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                     c=comp_dgtreal_long(f,gf,a,M);\n\
                     Yeah."


#include "ltfat_oct_template_helper.h"
// octave_idx_type is 32 or 64 bit signed integer
/*
  dgtreal_fb forwarders
*/

static inline void fwd_dgtreal_long(const double *f, const double *g,
                                  const octave_idx_type L,
							      const octave_idx_type W, const octave_idx_type a,
							      const octave_idx_type M, Complex *cout)
{
   dgtreal_long_d(f,g,L,W,a,M,reinterpret_cast<double (*)[2]>(cout));
}

static inline void fwd_dgtreal_long(const float *f, const float *g,
                                const octave_idx_type L,
							    const octave_idx_type W, const octave_idx_type a,
							    const octave_idx_type M, FloatComplex *cout)
{
   dgtreal_long_s(f,g,L,W,a,M,reinterpret_cast<float (*)[2]>(cout));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
     DEBUGINFO;
	 const octave_idx_type a = args(2).int_value();
     const octave_idx_type M = args(3).int_value();

	 MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));
	 MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(1));
  
     const octave_idx_type L = f.rows();
     const octave_idx_type W = f.columns();
     const octave_idx_type N = L/a;

     const octave_idx_type M2 = (M/2)+1;
	 
	 dim_vector dims_out(M2,N*W);  
     dims_out.chop_trailing_singletons();

     MArray<LTFAT_COMPLEX> cout(dims_out); 
     cout.fill(0);
	 
	 fwd_dgtreal_long(f.data(),g.data(),L,W,a,M,cout.fortran_vec());
	 
     return octave_value(cout);
}


/* OLD CODE
#include <octave/oct.h>
#include "config.h"
#include "ltfat.h"

DEFUN_DLD (comp_dgtreal_long, args, ,
  "This function calls the C-library\n\
  c=comp_dgtreal_long(f,gf,a,M);\n\
  Yeah.")
{

  const int a = args(2).int_value();
  const int M = args(3).int_value();

  const Matrix f = args(0).matrix_value();
  const Matrix g = args(1).matrix_value();
  
  const int L = f.rows();
  const int W = f.columns();
  const int N = L/a;

  const int M2 = (M/2)+1;
  
  ComplexMatrix cout(M2,N*W);  
  
  dgtreal_long((double*)f.data(),(double*)g.data(),
	    L, W, a, M,(ltfat_complex*)cout.data());
  
  return octave_value (cout);

}
*/