function [fhat, fola] = blocksyn(F, c , Lb, fola)
%BLOCKSYN Blockwise synthesis interface
%   Usage: blocksyn(F, c, Lb)
%
%   Input parameters:
%      Fs   : Synthesis frame object.
%      c    : Coefficients of a block.
%      fola : Explicitly defined overlap.
%   Output parameters:
%      fhat : Reconstructed block of signal.
%      fola : Stored overlap.
%
%   c=BLOCKSYN(F,c,Lb) reconstructs the signal block fhat from the coefficients c*
%   using the frame defined by F.
%
%   *Note:* To get perfect reconstruction, the synthesis frame F must
%   be a dual frame of the analysis frame used in BLOCKANA.
%
%   See also: block, blockana, framedual   
%
%   References:
%     N. Holighaus, M. Dörfler, G. A. Velasco, and T. Grill. A framework for
%     invertible, real-time constant-Q transforms. IEEE Transactions on
%     Audio, Speech and Language Processing, 21(4):775 -785, 2013.
%     
%     Z. Průša. Segmentwise Discrete Wavelet Transform. PhD thesis, Brno
%     University of Technology, Brno, 2012.
%     
%
%   Url: http://ltfat.sourceforge.net/doc/blockproc/blocksyn.php

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
    
    if ~isstruct(F)
        error('%s: First agument must be a frame definition structure.',upper(mfilename));
    end;
    
    if nargin<4
       fola = [];
    end
    
    % Next block index start (from a global point of view, starting with zero)
    nextSb = block_interface('getPos');
    % Block index start (from a global point of view, starting with zero)
    Sb = nextSb-Lb;
    
    if strcmp(F.blockalg,'naive')
       % Most general. Should work for anything.
       % Produces awful block artifacts when coefficients are altered.
       fhat = F.frsyn(c);
       fhat = fhat(1:Lb,:);
    elseif strcmp(F.blockalg,'sliced')
        % General processing
        % Equal block length assumtion
        % Reconstruct
        f = F.frsyn(c);
        % Result should not be longer than 2*Lb
        f = f(1:2*Lb,:);
        % Multiply by a slicing window
        f = bsxfun(@times,F.sliwin,f);
        % Load and add overlap (first half)
        ol = loadOverlap(Lb,size(c,2),fola);
        olLen = size(ol,1);
        f(1:olLen,:) = f(1:olLen,:) + ol;
        % Store overlap (second half)
        if nargout>1
           fola=storeOverlap(f,Lb);
        else
           storeOverlap(f,Lb);
        end
        % Return first half
        fhat = f(1:Lb,:);
    elseif strcmp(F.blockalg,'segola')
       Lw = F.winLen;
       switch(F.type)
         case 'fwt'
           % The SegDWT algorithm
           J = F.J;
           w = F.g;
           m = numel(w.g{1}.h);
           a = w.a(1);
           blocksize = a^J;
           r = (a^J-1)/(a-1)*(m-1);
           Lbrec = (floor(nextSb/blocksize) - floor(Sb/blocksize))*blocksize;
           rSb = (a^J-1)/(a-1)*(m-a) + mod(Sb,a^J);
           over = r - rSb;
           f = block_ifwt(c,w,J,Lbrec);
           ol = loadOverlap(r-mod(Sb, a^J),size(c,2),fola);
           olLen = size(ol,1);
           f(1:olLen-over,:) = f(1:olLen-over,:) + ol(1+over:end,:);
           f = [ol(1:over,:);f];
           if nargout>1
              fola=storeOverlap(f,r-mod(nextSb, a^J));
           else
              storeOverlap(f,r-mod(nextSb, a^J));
           end
           fhat = f(1:Lb,:);
          case {'dgt','dgtreal'}
           % Time step 
           a = F.a; 
           % Length of the left half of the window
           Lwl = floor(Lw/2);
           % Length of the right part of the window
           Lwr = floor(Lw/2);

           Sbonelmax =  ceil((Lw-1)/a)*a + a-1;
           Sbolen = ceil((Lw-1)/a)*a + mod(Sb,a);
           % Next block overlap length
           nextSbolen = ceil((Lw-1)/a)*a + mod(nextSb,a);
           Lext = Sbolen + Lb - mod(nextSb,a);
           Lextc = Sbolen + Lb - nextSbolen + Lwl;
           
           startc = ceil(Lwl/a)+1;
           endc = ceil((Lextc)/a);

           
           cc = F.coef2native(c,size(c));
           chat = zeros(size(cc,1),ceil(Lext/a),size(cc,3));
           chat(:,startc:endc,:) = cc;
           chat = F.native2coef(chat); 
           f = F.frsyn(chat);
           %startIdx = Sbcolen*a-Lwr;
           %f = f(1+startIdx:startIdx+Lext,:);
           over = Sbonelmax - Sbolen;
           
           ol = loadOverlap(Sbonelmax-mod(Sb,a),size(c,2),fola);
           olLen = size(ol,1);
           f(1:olLen-over,:) = f(1:olLen-over,:) + ol(1+over:end,:);
           f = [ol(1:over,:);f];
           if nargout>1
              fola=storeOverlap(f,Sbonelmax-mod(nextSb,a));
           else
              storeOverlap(f,Sbonelmax-mod(nextSb,a));
           end
           fhat = f(1:Lb,:);
             
          otherwise
           error('%s: Unsupported frame.',upper(mfilename));
        end

    else
       error('%s: Frame was not created with blockaccel.',upper(mfilename));
    end

end

function overlap = loadOverlap(L,chan,overlap)
%LOADOVERLAP Loads overlap
%
%
   if isempty(overlap)
      overlap = block_interface('getSynOverlap');
   end
   
   if isempty(overlap)
     overlap = zeros(L,chan,block_interface('getClassId'));
   end
   Lo = size(overlap,1);
   if nargin<1
      L = Lo;
   end
   if L>Lo
      error('%s: Required more samples than stored.',upper(mfilename));
   end
   overlap = overlap(end-L+1:end,:);
end

function overlap = storeOverlap(fext,L)
%STOREOVERLAP Stores overlap
%
%
    if L>size(fext,1)
        error('%s: Storing more samples than passed.',upper(mfilename));
    end
    overlap = fext(end-L+1:end,:);
    
    if nargout<1
       block_interface('setSynOverlap',overlap); 
    end
end % STOREOVERLAP

