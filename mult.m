%Copyright (C) 2022 Théodore CHERRIERE (theodore.cherriere@ens-paris-saclay.fr)

%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.

%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.

%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <https://www.gnu.org/licenses/>.

function result=mult(varargin)
% result=mult(A1,A2,...,AN)
%__________________________________________________________________________
% Compute the matrice product A(:,:,i,j..)*B(:,:,i,j..)*...*N(:,:,i,j..) 
% efficiently with pagemtimes (R2020B or later) or mmx for older releases
% (https://github.com/yuvaltassa/mmx)
%__________________________________________________________________________
% EXAMPLE :
%   
% A = rand(3,2,1e6);
% B = rand(2,3,1e6);
% C = rand(3,5,1e6);
%
% tic
% r1=zeros(3,5,1e6);
% for i=1:size(A,3)
%   r1(:,:,i)=A(:,:,i)*B(:,:,i)*C(:,:,i);
% end
% toc
%
% tic
% r2=mult(A,B,C);
% toc
% max(abs(r1-r2),[],'all')

vv = ver('MATLAB');
vvNum = str2double(vv.Release(3:end-2));
vvId = vv.Release(end-1:end-1);

vvOld=false;

if vvNum<2020
    vvOld=True;
elseif vvNum==2020 && vvId < 'b'
    vvOld=True;
end

result=varargin{1};

if vvOld
     for i=2:length(varargin)
        result=mmx('mult',resultat,varargin{i}); 
     end
    % check https://github.com/yuvaltassa/mmx to get MMX
else
    for i=2:length(varargin)
        result=pagemtimes(result,varargin{i});
    end
end

end

    