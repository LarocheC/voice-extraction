function [c,info]=uwfbt(f,wt,varargin)
%UWFBT   Undecimated Wavelet FilterBank Tree
%   Usage:  c=uwfbt(f,wt);
%           [c,info]=uwfbt(...);
%
%   Input parameters:
%         f   : Input data.
%         wt  : Wavelet Filterbank tree
%
%   Output parameters:
%         c     : Coefficients stored in L xM matrix.
%
%   c=UWFBT(f,wt) computes redundant time (or shift) representation c 
%   of the input signal f using the filterbank tree definition in wt and
%   using the "a-trous" algorithm. Number of columns in c (*M*) is defined
%   by the total number of outputs of the leaf nodes of the tree.
%   In addition, the function returns struct. info containing the transform
%   parameters. It can be conviniently used for the inverse transform IUWFBT
%   e.g. fhat = iUWFBT(c,info). It is also required by the PLOTWAVELETS
%   function.
%
%   If f is a matrix, the transformation is applied to each of W columns
%   and the coefficients in c are stacked along the third dimension.
%
%   Please see help for WFBT description of possible formats of wt and
%   description of frequency and natural ordering of the coefficient subbands.
%
%   Examples:
%   ---------
%   
%   A simple example of calling the UWFBT function using the "full decomposition" wavelet tree:
% 
%     f = greasy;
%     J = 8;
%     [c,info] = uwfbt(f,{'sym10',J,'full'});
%     plotwavelets(c,info,16000,'dynrange',90);
%
%   See also: iuwfbt, wfbtinit
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/uwfbt.php

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

definput.import = {'wfbtcommon'};
flags=ltfatarghelper({},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.forder);

%% ----- step 1 : Verify f and determine its length -------
[f,Ls]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',upper(mfilename));  
end

%% ----- step 2 : Prepare input parameters
wtPath = nodesBForder(wt);
nodesUps = nodeFiltUps(wtPath,wt);
rangeLoc = rangeInLocalOutputs(wtPath,wt);
rangeOut = rangeInOutputs(wtPath,wt);
%% ----- step 3 : Run computation
c = comp_uwfbt(f,wt.nodes(wtPath),nodesUps,rangeLoc,rangeOut);

%% ----- Optional : Fill the info struct. -----
if nargout>1
   info.fname = 'uwfbt';
   info.wt = wt;
   info.fOrder = flags.forder;
   info.isPacked = 0;
end

