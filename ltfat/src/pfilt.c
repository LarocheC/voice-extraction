#include "config.h"
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include "ltfat.h"


LTFAT_EXTERN void
LTFAT_NAME(pfilt_fir_rr)(const LTFAT_REAL *f, const LTFAT_REAL *g,
	      const int L, const int gl,
	      const int W, const int a,
	      LTFAT_REAL *cout)
{
   /*  --------- initial declarations -------------- */

   int l, n, w;

   LTFAT_REAL *gw;

   LTFAT_REAL *gb;
   LTFAT_REAL fw;

   const LTFAT_REAL *fbd;
   
   /*  ----------- calculation of parameters and plans -------- */
   
   const int N=L/a;

   /* These are floor operations. */
   const int glh=gl/2;

   /* This is a ceil operation. */
   const int glh_d_a=(int)ceil((glh*1.0)/(a));
   
   gw   = (LTFAT_REAL*)ltfat_malloc(gl*sizeof(LTFAT_REAL));
   
   /*  ---------- main body ----------- */
  
   /* Do the fftshift of g to place the center in the middle and
    * conjugate it.
    */
   
   for (l=0;l<glh;l++)
   {
      gw[l]=g[l+(gl-glh)];
   }
   for (l=glh;l<gl;l++)
   {
      gw[l]=g[l-glh];
   }

   for (w=0;w<W;w++)
   {	 
      /*----- Handle the first boundary using periodic boundary conditions.*/
      for (n=0; n<glh_d_a; n++)
      {
	 gb=gw;
	 fbd=f+L-(glh-n*a)+L*w;
	 fw=0.0;
	 for (l=0;l<glh-n*a;l++)
	 {
	    fw +=fbd[l]*gb[l];
	 }
	 fbd=f-(glh-n*a)+L*w;
	 for (l=glh-n*a;l<gl;l++)
	 {
	    fw +=fbd[l]*gb[l];
	 }
	 cout[n+w*N]=fw;
      }	     
   
      /* ----- Handle the middle case. --------------------- */
      for (n=glh_d_a; n<(L-(gl+1)/2)/a+1; n++)
      {
	 gb=gw;
	 fbd=f+(n*a-glh+L*w);
	 fw=0;
	 for (l=0;l<gl;l++)
	 {
	    fw +=fbd[l]*gb[l];
	 }
	 cout[n+w*N]=fw;	 
      }
      
      /* Handle the last boundary using periodic boundary conditions. */   
      for (n=(L-(gl+1)/2)/a+1; n<N; n++)
      {
	 gb=gw;
	 fbd=f+(n*a-glh+L*w);
	 fw=0;
	 for (l=0;l<L-n*a+glh;l++)
	 {
	    fw +=fbd[l]*gb[l];
	 }
	 fbd=f-(L-n*a+glh)+L*w;
	 for (l=L-n*a+glh;l<gl;l++)
	 {	
	    fw +=fbd[l]*gb[l];
	 }
	 cout[n+w*N]=fw;
      }
   }

    /* -----------  Clean up ----------------- */   
   ltfat_free(gw);

}
