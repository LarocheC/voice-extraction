function demo_blockproc_basicloop(source,varargin)
%DEMO_BLOCKPROC_BASICLOOP Basic real-time audio manipulation
%   Usage: demo_blockproc_basicloop('gspi.wav')
%
%   For additional help call DEMO_BLOCKPROC_BASICLOOP without arguments.
%
%   The demo runs simple playback loop allowing to set gain in dB.
%
%
%   Url: http://ltfat.sourceforge.net/doc/demos/demo_blockproc_basicloop.php

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

if demo_blockproc_header(mfilename,nargin)
   return;
end

% Basic Control pannel (Java object)
p = blockpanel({
               {'GdB','Gain',-20,20,0,21},...
               });

bufLen = 1024;
% Setup blocktream
if isoctave
   block(source,varargin{:},'L',bufLen);
else
   block(source,varargin{:},'loadind',p,'L',bufLen);
end

flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
   gain = blockpanelget(p,'GdB');
   gain = 10^(gain/20);
   
   [f,flag] = blockread();
   blockplay(f*gain);
end
p.close();
