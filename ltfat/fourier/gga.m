function c = gga(f,indvec,dim)
%GGA Generalized Goertzel algorithm
%   Usage:  y = gga(x,indvec)
%
%   Input parameters:
%         x      : Input data.
%         indvec : Indices to calculate. 
%
%   Output parameters:
%         c      : Coefficient vector.
%
%   c=GGA(f,indvec) computes the discrete-time fourier transform DTFT of
%   f at 'indices' contained in indvec, using the generalized second-order
%   Goertzel algorithm. Thanks to the generalization, the 'indices' can be 
%   non-integer valued in the range 0 to Ls-1, where Ls is the length of
%   the first non-singleton dimension of f. Index 0 corresponds to the
%   DC component and integers in indvec result in the classical DFT
%   coefficients. If indvec is empty or ommited, indvec is assumed to be
%   0:Ls-1.
%
%   c=GGA(f,indvec,dim) computes the DTFT samples along the dimension dim.
%
%   *Remark:**
%   Besides the generalization the algorithm is also shortened by one
%   iteration compared to the conventional Goertzel.
%
%   Examples:
%   ---------
%   
%   Calculating DTFT samples of interest:
% 
%     % Generate input signal
%     k = 0:2^10-1;
%     f = 5  sin(2*pi*k*0.05 + pi/4) + 2  sin(2*pi*k*0.1031 - pi/3);
%
%     % Non-integer indices of interest
%     kgga = 102.9:0.05:109.1;
%     % For the purposes of plot, remove the integer-valued elements
%     kgga = setdiff(kgga,k);
%
%     % This is equal to fft(f)
%     ck = gga(f,k);
%
%     fprintf('GGA to FFT error: %dn',norm(ck-fft(f)));
%
%     % DTFT samples just for non-integer indices
%     ckgga = gga(f,kgga);
%
%     % Plot modulus of coefficients
%     figure(f1);
%     hold on;
%     stem(k,abs(ck),'k');
%     stem(kgga,abs(ckgga),'r:');
%     limX = [102.9 109.1];
%     set(gca,'XLim',limX);
%     set(gca,'YLim',[0 1065]);
%
%     figure(f2);
%     hold on;
%     stem(k,angle(ck),'k');
%     stem(kgga,angle(ckgga),'r:');
%     set(gca,'XLim',limX);
%     set(gca,'YLim',[-pi pi]);
%
%   References:
%     P. Sysel and P. Rajmic. Goertzel algorithm generalized to non-integer
%     multiples of fundamental frequency. EURASIP Journal on Advances in
%     Signal Processing, 2012(1):56, 2012.
%     
%
%   Url: http://ltfat.sourceforge.net/doc/fourier/gga.php

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
       
% The original copyright goes to
% 2013 Pavel Rajmic, Brno University of Technology, Czech Rep.


%% Check the input arguments
if nargin < 1
    error('%s: Not enough input arguments.',upper(mfilename))
end

if isempty(f)
    error('%s: X must be a nonempty vector or a matrix.',upper(mfilename))
end

if nargin<3
  dim=[];  
end;

[f,~,Ls,~,dim,permutedsize,order]=assert_sigreshape_pre(f,[],dim,'GGA');

if nargin > 1 && ~isempty(indvec)
   if ~isreal(indvec) || ~isvector(indvec)
      error('%s: INDVEC must be a real vector.',upper(mfilename))
   end
else
   indvec = 0:Ls-1;
end

c = comp_gga(f,indvec);

permutedsize(1)=numel(indvec);

c=assert_sigreshape_post(c,dim,permutedsize,order);


