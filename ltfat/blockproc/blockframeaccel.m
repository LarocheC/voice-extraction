function Fo = blockframeaccel(F, Lb, varargin)
%BLOCKFRAMEACCEL Precompute structures for block processing
%   Usage: F = blockframeaccel(F,Lb);
%
%   F=BLOCKFRAMEACCEL(F,Lb) have to be called for each frame object prior 
%   entering the main loop where BLOCKANA and BLOCKSYN are called.
%   The function work entirely like FRAMEACCEL but in addition, it prepares
%   structures for the processing of a consecutive stream of blocks.
%
%      'sliwin',sliwin   : Slicing window. sliwin have to be a window
%                            of length 2L. It is used in the slicing
%                            window approach.
%
%   Url: http://ltfat.sourceforge.net/doc/blockproc/blockframeaccel.php

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



definput.flags.blockalg = {'naive','sliced','segola'};
definput.keyvals.sliwin = [];
[flags,kv]=ltfatarghelper({},definput,varargin);

assert(~(~flags.do_sliced && ~isempty(kv.sliwin)),...
   '%s: Definig slicing window without setting the ''silced'' flag.',mfilename);

if flags.do_sliced 
   if isempty(kv.sliwin)
      kv.sliwin = 'hann';
   end

   if ~isnumeric(kv.sliwin)
      kv.sliwin = fftshift(firwin(kv.sliwin,2*Lb));
   else
      if numel(kv.sliwin)~=2*Lb
         error('%s: The slicing window length has to be 2*Lb=%i.',upper(mfilename),2*Lb);
      end
   end

   Fo = frameaccel(F,2*Lb);
   Fo.sliwin = kv.sliwin;
elseif flags.do_segola
   Fo = frameaccel(F,Lb);

   if ~isfield(Fo,'winLen')
      error(['%s: Segment overlap cannot be used with this analysis frame.,'...
             ' It does not have FIR windows.'],upper(mfilename));
   end
    
   switch(Fo.type) 
      case {'dgt','dgtreal'}
         Fo = frameaccel(F,Lb+Fo.winLen-1+Fo.a);
         assert(Fo.a <= Lb ,['%s: Time step is bigger than the',...
      ' block length.'],mfilename);
      case {'filterbank','filterbankreal','ufilterbank','ufilterbankreal'}
         Fo = frameaccel(F,Lb+Fo.winLen-1+max(Fo.a(:,1)));
         assert(all(Fo.a(:,2)==1), '%s: Fractional subsampling is not supported',upper(mfilename) );
         assert(max(Fo.a(:,1)) <= Lb ,['%s: Time step is bigger than the',...
      ' block length.'],upper(mfilename));
         Fo.lcma =  filterbanklength(1,F.a);
      case {'dwilt'}
         Fo = frameaccel(F,Lb+Fao.winLen-1+2*Fo.M);
         Fo.a = 2*Fo.M;
      case {'wmdct'}
         Fo = frameaccel(F,Lb+Fao.winLen-1+Fo.M);
         Fo.a = Fo.M;
   end
   


elseif flags.do_naive
   Fo = frameaccel(F,Lb);
end

Fo.blockalg = flags.blockalg;
