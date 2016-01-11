function f = comp_ifilterbank(c,g,a,L)
%COMP_IFILTERBANK Compute inverse filterbank
%
%   Url: http://ltfat.sourceforge.net/doc/comp/comp_ifilterbank.php

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

W = size(c{1},2);
M = numel(g);
classname = assert_classname(c{1});

f=zeros(L,W,classname);

% This routine must handle the following cases
%
%   * Time-side or frequency-side filters (test for  isfield(g,'H'))
%
%   * Cell array or matrix input (test for iscell(c))
%
%   * Regular or fractional subsampling (test for info.isfractional)

l=(0:L-1).'/L;
for m=1:M
    conjG=conj(comp_transferfunction(g{m},L));
        
    % For Octave 3.6 compatibility
    conjG=cast(conjG,classname);
    
    % Handle fractional subsampling (this implies frequency side filters)
    if isfield(g{m},'H') && numel(g{m}.H)~=L
        N=size(c{m},1);
        Llarge=ceil(L/N)*N;
        amod=Llarge/N;
        
        for w=1:W                        
            % This repmat cannot be replaced by bsxfun
            innerstuff=middlepad(circshift(repmat(fft(c{m}(:,w)),amod,1),-g{m}.foff),L);
            f(:,w)=f(:,w)+ifft(circshift(bsxfun(@times,innerstuff,circshift(conjG,-g{m}.foff)),g{m}.foff));
        end;                
    else
        if iscell(c)
            for w=1:W
                % This repmat cannot be replaced by bsxfun
                f(:,w)=f(:,w)+ifft(repmat(fft(c{m}(:,w)),a(m),1).*conjG);
            end;
        else
            for w=1:W
                % This repmat cannot be replaced by bsxfun
                f(:,w)=f(:,w)+ifft(repmat(fft(c(:,m,w)),a(m),1).*conjG);
            end;            
        end;
    end;
end;
