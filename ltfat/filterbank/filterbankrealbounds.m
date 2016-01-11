function [AF,BF]=filterbankrealbounds(g,a,L);
%FILTERBANKREALBOUNDS  Frame bounds of filter bank for real signals only
%   Usage: fcond=filterbankrealbounds(g,a);
%          [A,B]=filterbankrealbounds(g,a);
%
%   FILTERBANKREALBOUNDS(g,a,L) calculates the ratio B/A of the frame
%   bounds of the filterbank specified by g and a for a system of length
%   L. The ratio is a measure of the stability of the system.  Use this
%   function on the common construction where the filters in g only covers
%   the positive frequencies.
%
%   [A,B]=FILTERBANKREALBOUNDS(g,a) returns the lower and upper frame
%   bounds explicitly.
%
%   See also: filterbank, filterbankdual
%
%   Url: http://ltfat.sourceforge.net/doc/filterbank/filterbankrealbounds.php

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
  
if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if L~=filterbanklength(L,a)
    error(['%s: Specified length L is incompatible with the length of ' ...
           'the time shifts.'],upper(mfilename));
end;

[g,info]=filterbankwin(g,a,L,'normal');
M=info.M;

AF=Inf;
BF=0;
  
if info.isuniform
  % Uniform filterbank, use polyphase representation
  a=a(1);
  
  N=L/a;

  % G1 is done this way just so that we can determine the data type.
  G1=comp_transferfunction(g{1},L);
  thisclass=assert_classname(G1);
  G=zeros(L,M,thisclass);
  G(:,1)=G1;
  for ii=2:M
    G(:,ii)=comp_transferfunction(g{ii},L);
  end;
  
  Ha=zeros(a,M,thisclass);
  Hb=zeros(a,M,thisclass);
  
  for w=0:N-1
    idx_a = mod(w-(0:a-1)*N,L)+1;
    idx_b = mod((0:a-1)*N-w,L)+1;
    Ha = G(idx_a,:);
    Hb = conj(G(idx_b,:));
    
    % A 'real' is needed here, because the matrices are known to be
    % Hermitian, but sometimes Matlab/Octave does not recognize this.  
    work=real(eig(real(Ha*Ha'+Hb*Hb')));
    
    AF=min(AF,min(work));
    BF=max(BF,max(work));
    
  end;
  
  AF=AF/a;
  BF=BF/a;
  
else
    if info.ispainless
        % Compute the diagonal of the frame operator.
        f=comp_filterbankresponse(g,info.a,L,1);
        
        AF=min(f);
        BF=max(f);
    else
        error(['%s: There is no fast method to find the frame bounds of ' ...
               'this filterbank as it is neither uniform nor painless. ' ...
               'Please see FRAMEBOUNDS for an iterative method that can ' ...
               'solve the problem.'],upper(mfilename));                        
    end;        
end;

if nargout<2
  % Avoid the potential warning about division by zero.
  if AF==0
    AF=Inf;
  else
    AF=BF/AF;
  end;
end;
  


