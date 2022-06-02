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

if version

    result=varargin{1};
    for i=2:length(varargin)
        result=pagemtimes(result,varargin{i}); % for new Matlab versions
        %resultat=mmx('mult',resultat,varargin{i}); % else (check
        % https://github.com/yuvaltassa/mmx)
    end
end

    