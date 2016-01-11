function [AF,BF]=wfbtbounds(wt,L)
%FWTBOUNDS Frame bounds of DWT
%   Usage: fcond=fwtbounds(wt,L);
%          [A,B]=fwtbounds(wt,L);
%
%   WFBTBOUNDS(wt,L) calculates the ratio B/A of the frame bounds
%   of the filterbank specified by wt for a system of length
%   L. The ratio is a measure of the stability of the system. 
%
%   [A,B]=WFBTBOUNDS(wt,L) returns the lower and upper frame bounds
%   explicitly. 
%
%   See WFBT for explanation of parameter wt.
%
%   See also: wfbt, filterbankbounds
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/wfbtbounds.php

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


if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

wt = wfbtinit({'strict',wt},'nat');

if L~=wfbtlength(L,wt)
    error(['%s: Specified length L is incompatible with the length of ' ...
           'the time shifts.'],upper(mfilename));
end;


for ii=1:numel(wt.nodes)
   a = wt.nodes{ii}.a;
   assert(all(a==a(1)),'%s: One of the basic wavelet filterbanks is not uniform.',upper(mfilename));
end

% Do the equivalent filterbank
[g,a] = wfbt2filterbank(wt);

maxGlen = max(cellfun(@(gEl) numel(gEl.h),g));
% maxa is also lcm(a)
maxa = max(a);

assert(L>=maxGlen,'%s: One of the filters is longer than L.',upper(mfilename));

gnew = {};

for m=1:numel(g)
   pk = maxa/a(m);
   
   for ii=0:pk-1
      gnew{end+1} = struct('h',g{m}.h,'offset',g{m}.offset-a(m)*ii);
   end
end

anew = maxa*ones(numel(gnew),1);


if nargout<2
   AF = filterbankbounds(gnew,anew,L);
elseif nargout == 2
   [AF, BF] = filterbankbounds(gnew,anew,L);
end
