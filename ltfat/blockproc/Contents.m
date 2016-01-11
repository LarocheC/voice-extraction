% LTFAT - Block processing
%
%  Zdenek Prusa, 2013.
%
%  Basic methods
%    BLOCK          - Create a new block
%    BLOCKDEVICES   - List available devices
%    BLOCKREAD      - Read samples from file/device
%    BLOCKANA       - Block analysis
%    BLOCKSYN       - Block synthesis
%    BLOCKPLAY      - Play block (sound output)
%    BLOCKDONE      - Destroy the block object
%
%  Helper functions
%    BLOCK_FWT      - FWT processing
%    BLOCK_IFWT     - IFWT processing
%
%  For help, bug reports, suggestions etc. please send email to
%  ltfat-help@lists.sourceforge.net
%
%   Url: http://ltfat.sourceforge.net/doc/blockproc/Contents.php

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

