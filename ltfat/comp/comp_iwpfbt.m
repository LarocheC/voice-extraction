function f=comp_iwpfbt(c,wtNodes,pOutIdxs,chOutIdxs,Ls,ext)
%COMP_IWFBT Compute Inverse Wavelet Packet Filter-Bank Tree
%   Usage:  f=comp_iwpfbt(c,wtNodes,pOutIdxs,chOutIdxs,Ls,ext)
%
%   Input parameters:
%         c          : Coefficients stored in cell array.
%         wtNodes    : Filterbank tree nodes (elementary filterbans) in
%                      reverse BF order. Cell array of structures of length nodeNo.
%         pOutIdxs   : Idx of each node's parent. Array of length nodeNo.
%         chOutIdxs  : Idxs of each node children. Cell array of vectors of
%                      length nodeNo.
%         ext        : Type of the forward transform boundary handling.
%
%   Output parameters:
%         f          : Reconstructed data in L*W array.
%
%
%   Url: http://ltfat.sourceforge.net/doc/comp/comp_iwpfbt.php

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

% Do non-expansve transform if ext=='per'
doPer = strcmp(ext,'per');

% For each node in tree in the BF order...
 for jj=1:length(wtNodes)
    % Node filters to a cell array
    gCell = cellfun(@(gEl) gEl.h(:),wtNodes{jj}.g(:),'UniformOutput',0);
    % Node filters subs. factors
    a = wtNodes{jj}.a;
    % Node filters initial skips
    if(doPer)
       offset = cellfun(@(gEl) gEl.offset,wtNodes{jj}.g);
    else
       offset = -cellfun(@(gEl) numel(gEl),gCell) + 1 + (a - 1);
    end
    
    if(pOutIdxs(jj))
       % Run filterbank and add to the existing subband.
       ctmp = comp_ifilterbank_td(c(chOutIdxs{jj}),gCell,a,size(c{pOutIdxs(jj)},1),offset,ext);
       c{pOutIdxs(jj)} = (1/sqrt(2))*(c{pOutIdxs(jj)}+ctmp);
    else
       % We are at the root.
       f = comp_ifilterbank_td(c(chOutIdxs{jj}),gCell,a,Ls,offset,ext);
    end
 end
     
 
