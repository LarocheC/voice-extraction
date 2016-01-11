function [c,info]=wfbt(f,wt,varargin)
%WFBT   Wavelet FilterBank Tree
%   Usage:  c=wfbt(f,wt);
%           [c,info]=wfbt(...);
%
%   Input parameters:
%         f     : Input data.
%         wt    : Wavelet Filterbank tree definition.
%
%   Output parameters:
%         c    : Coefficients stored in a cell-array.
%         info : Transform parameters struct.
%
%   c=WFBT(f,wt) returns coefficients c obtained by applying a wavelet 
%   filterbank tree defined by wt to the input data f. In addition, the
%   function returns struct. info containing transform parameters. It can
%   be conviniently used for the inverse transform IWFBT e.g. 
%   fhat = iWFBT(c,info). It is also required by the PLOTWAVELETS function.
%
%   wt defines a tree shaped filterbank structure build from the 
%   elementary two (or more) channel wavelet filters. The tree can have any
%   shape and thus provide a flexible frequency covering. The outputs of the
%   tree leaves are stored in c.
%
%   The wt parameter can have two formats:
%
%   1) Cell array containing 3 elements {w,J,treetype}, where w is
%      the basic wavelet filterbank definition as in FWT function, J*
%      stands for the depth of the tree and the flag treetype defines 
%      the type of the tree to be used. Supported options are:
%
%      'dwt'  
%        DWT tree. Just the low-pass output is decomposed further.
%
%      'full'
%        Full decomposition tree. Each output is decomposed up to level J.
%
%   2) Structure returned by the WFBTINIT function and possibly
%      modified by WFBTPUT and WFBTREMOVE.
%
%   If f is row/column vector, the coefficient vectors c{jj} are columns.
%   
%   If f is a matrix, the transformation is by default applied to each of
%   W columns [Ls, W]=size(f).
%
%   In addition, the following flag groups are supported:
%
%   'per','zero','odd','even'
%     Type of the boundary handling.
%
%   'freq','nat'
%     Frequency or natural ordering of the coefficient subbands. The direct
%     usage of the wavelet tree ('nat' option) does not produce coefficient
%     subbans ordered according to the frequency. To achieve that, some 
%     filter shuffling has to be done ('freq' option).  
%
%   Please see the help on FWT for a description of the boundary condition flags.
%
%   Examples:
%   ---------
%   
%   A simple example of calling the WFBT function using the "full decomposition" wavelet tree:
% 
%     f = gspi;
%     J = 7;
%     [c,info] = wfbt(f,{'sym10',J,'full'});
%     plotwavelets(c,info,44100,'dynrange',90);
%
%   See also: iwfbt, wfbtinit
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/wfbt.php

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


if(nargin<2)
   error('%s: Too few input parameters.',upper(mfilename));  
end

definput.import = {'fwt','wfbtcommon'};
definput.keyvals.dim = [];
[flags,kv]=ltfatarghelper({},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.forder);
    
%% ----- step 1 : Verify f and determine its length -------
[f,Ls]=comp_sigreshape_pre(f,upper(mfilename),0);

% Determine next legal input data length.
L = wfbtlength(Ls,wt,flags.ext);

% Pad with zeros if the safe length L differ from the Ls.
if(Ls~=L)
   f=postpad(f,L); 
end

%% ----- step 3 : Run computation
wtPath = nodesBForder(wt);
rangeLoc = rangeInLocalOutputs(wtPath,wt);
rangeOut = rangeInOutputs(wtPath,wt); % very slow
c = comp_wfbt(f,wt.nodes(wtPath),rangeLoc,rangeOut,flags.ext);

%% ----- Optionally : Fill info struct ----
if nargout>1
   info.fname = 'wfbt';
   info.wt = wt;
   info.ext = flags.ext;
   info.Lc = cellfun(@(cEl) size(cEl,1),c);
   info.Ls = Ls;
   info.fOrder = flags.forder;
   info.isPacked = 0;
end



