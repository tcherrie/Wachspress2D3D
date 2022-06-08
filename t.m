% Copyright (C) 2022 Théodore CHERRIERE (theodore.cherriere@ens-paris-saclay.fr)
% 
% This program is free software: you can redistribute it and/or modify  it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License  along with 
% this program.  If not, see <https://www.gnu.org/licenses/ <https://www.gnu.org/licenses/>>.

function tM=t(M)
%tM=t(M)
%__________________________________________________________________________
% Adapt pagetranspose depending on the Matlab version

vv = ver('MATLAB');
vvNum = str2double(vv.Release(3:end-2));
vvId = vv.Release(end-1:end-1);

vvOld=false;

if vvNum<2020
    vvOld=True;
elseif vvNum==2020 && vvId < 'b'
    vvOld=True;
end

if vvOld
    tM=permute(M,[2 1 3 4 5 6 7 8 9]); % if Matlab version < R2020B
else
    tM=pagetranspose(M);
end
end