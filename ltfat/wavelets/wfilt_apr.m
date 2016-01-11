function [h,g,a,info] = wfilt_apr(N)
%WFILT_APR Almost Perfect Reconstruction Filter Bank for Non-redundant, Approximately Shift-Invariant, ComplexWavelet Transforms
%   Usage: [h,g,a] = wfilt_apr(N);
%
%   [h,g,a] = WFILT_APR(N) computes an almost perfect reconstruction
%   filter bank for non-redundant, approximately shift-invariant,
%   Complex Wavelet Transforms. Critically subsampled. 
%
%
%   References:
%     R. Hosseini and M. Vafadust. Almost perfect reconstruction filter bank
%     for non-redundant, approximately shift-invariant, complex wavelet
%     transforms. Journal of Wavelet Theory and Applications, 2(1):1-14,
%     2008.
%     
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/wfilt_apr.php

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

a= [3;3;3];
switch(N)
    case 1
% from the software package filters1.m
harr = [
-0.0195   (-0.0046 - 0.0018*i) (-0.0046 + 0.0018*i)
-0.0028   (0.0160 - 0.0205*i)  (0.0160 + 0.0205*i)
-0.0373   (-0.0050 - 0.0079*i) (-0.0050 + 0.0079*i)
-0.0271   (0.0530 - 0.0753*i)  (0.0530 + 0.0753*i)
 0.2735   (0.0504 + 0.0402*i)  (0.0504 - 0.0402*i)
 0.0127   (0.1180 - 0.1403*i)  (0.1180 + 0.1403*i)
 0.6495   (0.3612 + 0.3391*i)  (0.3612 - 0.3391*i) 
 0.0369   (-0.4554 + 0.4554*i)  (-0.4554 - 0.4554*i)
 0.6495   (-0.3391 - 0.3612*i)  (-0.3391 + 0.3612*i) 
 0.0127   (0.1403 - 0.1180*i)  (0.1403 + 0.1180*i) 
 0.2735   (-0.0402 - 0.0504*i) (-0.0402 + 0.0504*i)
-0.0271   (0.0753 - 0.0530*i)  (0.0753 + 0.0530*i)
-0.0373   (0.0079 + 0.0050*i)  (0.0079 - 0.0050*i)
-0.0028   (0.0205 - 0.0160*i)  (0.0205 + 0.0160*i)
-0.0195   (0.0018 + 0.0046*i)  (0.0018 - 0.0046*i)
];
    case 2
   harr = [
-0.0124  (-0.0002 + 0.0041*i)  (-0.0002 - 0.0041*i)
 0.0062  (-0.0095 + 0.0090*i)  (-0.0095 - 0.0090*i)
-0.0184  (-0.0014 + 0.0032*i)  (-0.0014 - 0.0032*i)
 0.0306  (-0.0280 + 0.0134*i)  (-0.0280 - 0.0134*i)
 0.0611  ( 0.0024 - 0.0351*i)  ( 0.0024 + 0.0351*i)
 0.0171  (-0.0304 + 0.0042*i)  (-0.0304 - 0.0042*i)
 0.0883  ( 0.0372 - 0.1169*i)  ( 0.0372 + 0.1169*i)
-0.2966  ( 0.0352 + 0.0177*i)  ( 0.0352 - 0.0177*i)
-0.0589  ( 0.0938 - 0.1716*i)  ( 0.0938 + 0.1716*i)
-0.6231  ( 0.3841 + 0.3107*i)  ( 0.3841 - 0.3107*i)
-0.1181  (-0.4446 + 0.4446*i)  (-0.4446 - 0.4446*i)
-0.6231  (-0.3107 - 0.3841*i)  (-0.3107 + 0.3841*i)
-0.0589  ( 0.1716 - 0.0938*i)  ( 0.1716 + 0.0938*i)
-0.2966  (-0.0177 - 0.0352*i)  (-0.0177 + 0.0352*i)
 0.0883  ( 0.1169 - 0.0372*i)  ( 0.1169 + 0.0372*i)
 0.0171  (-0.0042 + 0.0304*i)  (-0.0042 - 0.0304*i)
 0.0611  ( 0.0351 - 0.0024*i)  ( 0.0351 + 0.0024*i)
 0.0306  (-0.0134 + 0.0280*i)  (-0.0134 - 0.0280*i)
-0.0184  (-0.0032 + 0.0014*i)  (-0.0032 - 0.0014*i)
 0.0062  (-0.0090 + 0.0095*i)  (-0.0090 - 0.0095*i)
-0.0124  (-0.0041 + 0.0002*i)  (-0.0041 - 0.0002*i)  
        ];

    otherwise
        error('%s: No such Almost Perfect Reconstruction Filter Bank filter. ',upper(mfilename));
end;

h=mat2cell(harr.',[1,1,1],length(harr));
g = h;
info.istight = 1;

