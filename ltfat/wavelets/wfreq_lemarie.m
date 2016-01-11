function [H,G] = wfreq_lemarie(L)
%WFREQ_LEMARIE  Battle and Lemarie filters frequency resp. sampling
%   Usage: [H,G]=wfreq_lemarie(L)
%
%   Input parameters:
%         N     : Number of samples of the frequency response.
%
%   [H,G]=wfreq_lemaire(L) calculates L samples of the Battle and
%   Lemarie filters frequency responses.
%
%   References:
%     S. G. Mallat. A theory for multiresolution signal decomposition: The
%     wavelet representation. IEEE Trans. Pattern Anal. Mach. Intell.,
%     11(7):674-693, July 1989. [1]http ]
%     
%     Henvisninger
%     
%     1. http://dx.doi.org/10.1109/34.192463
%     
%
%
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/wfreq_lemarie.php

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

% Original copyright goes to:
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
% Author: Jose Martin Garcia
% e-mail: Uvi_Wave@tsc.uvigo.es


% frequency axis
w=[0:2*pi/L:2*pi*(1-1/L)];
w(1)=eps;
w(L/2+1)=w(L/2+1)+1e-15;

% calculation of frequency response of analysis lowpass filter 
num=0;den=0;
K=36;
for k=-K:K,
	num=1./((w+2*pi*k).^8)+num;
	den=1./((2*w+2*pi*k).^8)+den;
end
H = cell(2,1);
H{1}=sqrt(num./(2.^8*den));
H{1}(1)=1;

H{2} = fftshift(H{1});
G = cell(2,1);
G{1} = fliplr(H{1});
G{2} = fliplr(H{2});




