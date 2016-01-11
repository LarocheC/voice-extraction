function [f,valid] = blockread(L)
%BLOCKREAD Read one block from input
%   Usage: blockread(L)
%       
%   Input parameters:
%      L    : Number of samples.
%   Output parameters:
%      f     : Samples.
%      valid : Input data valid flag.
%
%   Function also control the playback, so it does not have to rely on
%   whether user called BLOCKPLAY.
% 
%   Block streaming uses several buffers to compensate for the processing
%   delay variation. 
%
%   Url: http://ltfat.sourceforge.net/doc/blockproc/blockread.php

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

persistent Lwav;
persistent clearStr;
persistent readTime;
persistent t1;
persistent t2;
%global delayLog;
%global delayLog2;

if nargin<1
   L = block_interface('getBufLen'); 
end

do_updateGUI = 0;
do_updateBAR = 0;

loadind = block_interface('getDispLoad');

if ischar(loadind)
   if strcmp('bar',loadind)
      do_updateBAR = 1;
   end
elseif isjava(loadind)
   do_updateGUI = 1;
else
   error('%s: Something went wrong. Should not ever get here.',upper(mfilename));
end

if do_updateBAR || do_updateGUI 
   if block_interface('getPageNo')>0
      procTime = toc(t1);
      res = 2;

      fs= playrec('getSampleRate');
      load = floor(100*(procTime+readTime)/(L/fs));
      
      if do_updateBAR
         msg = sprintf(['Load : |',repmat('*',1,ceil(min([load,100])/res)),repmat(' ',1,floor((100-min([load,100]))/res)),'| \n']);
         droppedStr = sprintf('Dropped samples: %i\n',playrec('getSkippedSampleCount'));
         fprintf([clearStr,msg,droppedStr]);
         clearStr = repmat(sprintf('\b'), 1, length(msg)+length(droppedStr));
      elseif do_updateGUI
         javaMethod('updateBar',loadind,load);
      else
         error('%s: Something went wrong. Should not ever get here.',upper(mfilename));
      end
   
      block_interface('setSkipped',playrec('getSkippedSampleCount'));
      if playrec('getSkippedSampleCount') > block_interface('getSkipped')
         block_interface('setSkipped',playrec('getSkippedSampleCount'));
      end

   else
      clearStr = '';
      %delayLog = [];
      %delayLog2 = [0];
      procTime = 0;
   end
   t2 = tic;
end

valid = 1;
source = block_interface('getSource');
pos = block_interface('getPos') +1; % convert to the Matlab indexing
block_interface('incPageNo');
pageNo = block_interface('getPageNo');
classid = block_interface('getClassId');

% Update sample counter
block_interface('setPos',pos+L-1); % convert back the Matlab indexing

%%%
%% REC, source is a mic/aux, no loopback
%
if strcmp(source,'rec')
   recChanList = block_interface('getRecChanList');
   
   if do_updateBAR || do_updateGUI
      readTime = toc(t2);
   end
   % Issue reading buffers up to max
   while block_interface('getEnqBufCount') <= block_interface('getBufCount')
      block_interface('pushPage', playrec('rec', L, recChanList));
   end
   pageList = block_interface('getPageList');
   % Block until the first page is loaded
   while(playrec('isFinished', pageList(1)) == 0)
   end
   % Read the data. Cast to the specified type
   f = cast(playrec('getRec',pageList(1)),classid);
   % Delete page
   playrec('delPage', pageList(1));
   % Throw away the page id
   block_interface('popPage');
   
%%%   
%% PLAYREC, source is a mic, loopback to an output
%
elseif strcmp(source,'playrec')
   recChanList = block_interface('getRecChanList');
   if pageNo<=1
      blockplay(zeros(L,numel(recChanList),classid));
   end
   
   % Enqueue already processed
   fhat = block_interface('getEnqueuedToPlay');
   if isempty(fhat)
      fhat = zeros(L,numel(recChanList),classid);
   end
   chanList = block_interface('getPlayChanList');

   % Copy input channel to all output chanels.
   fhat = repmat(fhat,1,numel(chanList));
   % Play and record
   block_interface('pushPage',playrec('playrec', fhat, chanList, -1, recChanList));
   
   if do_updateBAR || do_updateGUI
      readTime = toc(t2);
   end
   pageList = block_interface('getPageList');
   % Playback is block_interface('getBufCount') behind the input
   if block_interface('getPageNo') <= block_interface('getBufCount')
      f = zeros(L,numel(recChanList),classid);
   else
      % Block until the first page is loaded
      while(playrec('isFinished', pageList(1)) == 0)
      end
      % Read the data
      f = cast(playrec('getRec',pageList(1)),classid);
      playrec('delPage', pageList(1));
      % Throw away the page id
      block_interface('popPage');
   end

%%%   
%% PLAY: Source is a *.wav file
%
elseif isa(source,'function_handle')
   % Get play channel list (could be chached) 
   chanList = block_interface('getPlayChanList');
   % Get already processed (from blockplay)
   fhat = block_interface('getEnqueuedToPlay');

   % Create something if blockplay was not called
   if isempty(fhat)
      fhat = zeros(L,numel(chanList),classid);
   end

   % Broadcast single input channel to all output chanels.
   if size(fhat,2)==1
      fhat = repmat(fhat,1,numel(chanList));
   end

   % Number of wav samples (is chached, since it is disk read operation)
   Lwav = block_interface('getLs');


   % Determine valid samples
   endSample = min(pos + L - 1, Lwav(1));
   %f = cast(wavread(source,[pos, endSample]),block_interface('getClassId'));
   f = source(pos,endSample);
   % Pad with zeros if some samples are missing
   if (pos + L - 1) >= Lwav(1)
      ftmp = zeros(L,Lwav(2),classid);
      ftmp(1:size(f,1),:) = f;
      f = ftmp;
      % Rewind if loop option was set.
      if block_interface('isLoop')
         block_interface('setPos',0);
         % Throw away stored overlaps.
         if ~isempty(block_interface('getAnaOverlap'))
            block_interface('setAnaOverlap',[]);
         end
         if ~isempty(block_interface('getSynOverlap'))
            block_interface('setSynOverlap',[]);
         end
      else
         valid = 0;
      end
   end


   % playrec('play',... - enques fhat to be played
   % block_interface('pushPage', - stores page number in an inner FIFO
   % queue
   block_interface('pushPage', playrec('play', fhat, chanList));

   if do_updateBAR || do_updateGUI
      readTime = toc(t2);
   end
   % If enough buffers are enqued, block the execution here until the 
   % first one is finished.
   if block_interface('getEnqBufCount') > block_interface('getBufCount')
      pageId = block_interface('popPage');
      % "Aggresive" chceking if page was played.
      % Another (supposedly slower) option is:
      % playrec('block',pageId);
      while(playrec('isFinished', pageId) == 0), end;
   end
end

if pageNo<=1
   playrec('resetSkippedSampleCount');
end

if do_updateBAR || do_updateGUI
   t1=tic;
end


