function devs = blockdevices()
%BLOCKDEVICES Lists audio devices
%   Usage: devs = blockdevices();
%
%   BLOCKDEVICES lists the available audio input and output devices. The
%   ID can be used in the BLOCK function to specify which device should
%   be used.
%
%   See also: block
%
%   Url: http://ltfat.sourceforge.net/doc/blockproc/blockdevices.php

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

devs = playrec('getDevices');

fprintf('\nAvailable output devices:\n');

for k=1:length(devs)
    if(devs(k).outputChans)
        fprintf('ID =%2d: %s (%s) %d channels, default latency %d--%d ms\n', ...
            devs(k).deviceID, devs(k).name, ...
            devs(k).hostAPI, devs(k).outputChans,...
            floor(1000*devs(k).defaultLowOutputLatency),...
            floor(1000*devs(k).defaultHighOutputLatency));

    end
end

fprintf('\nAvailable input devices:\n');

for k=1:length(devs)
    if(devs(k).inputChans)
        fprintf('ID =%2d: %s (%s) %d channels,  default latency %d--%d ms\n', ...
            devs(k).deviceID, devs(k).name, ...
            devs(k).hostAPI, devs(k).inputChans,...
            floor(1000*devs(k).defaultLowInputLatency),...
            floor(1000*devs(k).defaultHighInputLatency));

    end
end
