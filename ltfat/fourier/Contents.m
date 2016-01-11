% LTFAT - Basic Fourier and DCT analysis.
%
%  Peter L. Søndergaard, 2008 - 2013.
%
%  Support routines
%    FFTINDEX       -  Index of positive and negative frequencies.
%    MODCENT        -  Centered modulo operation.
%    FLOOR23        -  Previous number with only 2,3 factors
%    FLOOR235       -  Previous number with only 2,3,5 factors
%    CEIL23         -  Next number with only 2,3 factors
%    CEIL235        -  Next number with only 2,3,5 factors
%    NEXTFASTFFT    -  Next efficient FFT size (2,3,5,7).
%  
%  Basic Fourier analysis
%    DFT            -  Unitary discrete Fourier transform.
%    IDFT           -  Inverse of DFT.
%    FFTREAL        -  FFT for real valued signals.
%    IFFTREAL       -  Inverse of FFTREAL.
%    PLOTFFT        -  Plot FFT coefficients.
%    PLOTFFTREAL    -  Plot FFTREAL coefficients.
%
%  Simple operations on periodic functions
%    INVOLUTE       -  Involution.
%    PEVEN          -  Even part of periodic function.
%    PODD           -  Odd part of periodic function.
%    PCONV          -  Periodic convolution.
%    CONVOLVE       -  Fast, non-periodic convolution.
%    PXCORR         -  Periodic cross correlation.
%    ISEVENFUNCTION -  Test if function is even.
%    MIDDLEPAD      -  Cut or extend even function.
%
%  Functions
%    EXPWAVE        -  Complex exponential wave.
%    PCHIRP         -  Periodic chirp.
%    SHAH           -  Shah distribution.
%    PHEAVISIDE     -  Periodic Heaviside function.
%    PRECT          -  Periodic rectangle function.
%    PSINC          -  Periodic sinc function.
%
%  Window functions
%    PGAUSS         -  Periodic Gaussian.
%    PSECH          -  Periodic SECH.
%    PBSPLINE       -  Periodic B-splines.
%    FIRWIN         -  FIR windows (Hanning,Hamming,Blackman,...).
%    FIRKAISER      -  FIR Kaiser-Bessel window.
%    FIR2LONG       -  Extend FIR window to LONG window.
%    LONG2FIR       -  Cut LONG window to FIR window.
%
%  Filtering
%    FIRFILTER      -  Construct an FIR filter.
%    BLFILTER       -  Construct a band-limited filter.
%    WARPEDBLFILTER -  Warped, band-limited filter.
%    PFILT          -  Apply filter with periodic boundary conditions.
%    MAGRESP        -  Magnitude response plot.
%    TRANSFERFUNCTION - Computer the transfer function of a filter.
%
%  Hermite functions and fractional Fourier transforms
%    PHERM          -  Periodic Hermite functions.
%    HERMBASIS      -  Orthonormal basis of Hermite functions.    
%    DFRACFT        -  Discrete Fractional Fourier transform
%    FFRACFT        -  Fast Fractional Fourier transform
%
%  Approximation of continuous functions
%    FFTRESAMPLE    -  Fourier interpolation.
%    DCTRESAMPLE    -  Cosine interpolation.
%    PDERIV         -  Derivative of periodic function.
%
%  Cosine and Sine transforms.
%    DCTI           -  Discrete cosine transform type I
%    DCTII          -  Discrete cosine transform type II
%    DCTIII         -  Discrete cosine transform type III
%    DCTIV          -  Discrete cosine transform type IV
%    DSTI           -  Discrete sine transform type I
%    DSTII          -  Discrete sine transform type II
%    DSTIII         -  Discrete sine transform type III
%    DSTIV          -  Discrete sine transform type IV
%
%  For help, bug reports, suggestions etc. please send email to
%  ltfat-help@lists.sourceforge.net
%
%   Url: http://ltfat.sourceforge.net/doc/fourier/Contents.php

% Copyright (C) 2005-2013 Peter L. Søndergaard <soender@users.sourceforge.net>.
% This file is part of LTFAT version 1.4.2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


