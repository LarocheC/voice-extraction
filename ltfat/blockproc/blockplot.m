function blockplot(p,F,c)
%BLOCKPLOT Plot block coefficients
%   Usage: blockplot(p,F,c);
%
%   Input parameters:
%         p     : JAVA object of the class net.sourceforge.ltfat.SpectFrame.
%         F     : Frame object.
%         c     : Block coefficients.
%
%   BLOCKPLOT(p,F,c) appends the block coefficients c to the running 
%   coefficient plot in p.
%
%   Url: http://ltfat.sourceforge.net/doc/blockproc/blockplot.php

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

if size(c,2)>1
   error('%s: Only one channel input is supported.',upper(mfilename));
end

ctf = framecoef2tf(F,c(:,1));

if strcmp(F.blockalg,'sliced')
   % DO the coefficient overlapping or cropping
end

ctf = cast(ctf,'single');
javaMethod('append',p,ctf);



