function demo_blockproc_slidingsgram(source,varargin)
%DEMO_BLOCKPROC_SLIDINGSGRAM Basic real-time rolling spectrogram visualization
%   Usage: demo_blockproc_slidingsgram('gspi.wav')
%
%   For additional help call DEMO_BLOCKPROC_SLIDINGSGRAM without arguments.
%
%   This demo shows a simple rolling spectrogram of whatever is specified in
%   source. 
%
%   Url: http://ltfat.sourceforge.net/doc/demos/demo_blockproc_slidingsgram.php

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

% Basic sepctrogram figure (Java object)
fobj = blockfigure();

%
bufLen = 1024;
% Setup blocktream
if isoctave
   fs=block(source,varargin{:},'L',bufLen);
else
   fs=block(source,varargin{:},'loadind',p,'L',bufLen);
end

% Using dgtreal with 20ms hann window, hop factor 80, 1000 channels.
% Redundancy factor 12.5
winLenms = 20;
a = 80;
M = 1000;

F = frame('dgtreal',{'hann',floor(fs*winLenms/1e3)},a,M);
F = blockframeaccel(F, bufLen,'segola');


flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
   % Obtain the global gain value
   gain = blockpanelget(p,'GdB');
   gain = 10^(gain/20);

   % Read the next block of samples
   [f,flag] = blockread();
   f=f*gain;
   
   % Do analysis using the specified frame. 
   c = blockana(F, f); 
   
   % Draw the first channel coefficients
   blockplot(fobj,F,c(:,1));
   
   % Enqueue to play
   blockplay(f);
end
p.close();
fobj.close();
