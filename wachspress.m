function [w,dwdx]=wachspress(x,domain)
% [w,dwdx]=wachspress(x,domain)
% INPUTS :
%   - x : the vector of the cartesian coordinates
%   - domain : interpolation domain
% OUTPUTS :
%   - w : Wachspress shape functions
%   - dwdx : gradients of w

[n,dim]=size(x);
domain.vertices=domain.vertices-sum(domain.vertices,1)/length(domain.vertices);
if dim==0
    w=ones(1,1,n);
    dwdx=zeros(0,1,n);
elseif dim==1
    x=shiftdim(x,-2);
    w=[-x;x]+0.5; % assume that 0 is the barycenter
    dwdx=repmat([-1;1],[1,1,n]);
elseif dim==2
    [w,dwdx]= wachspress2d_opt(x,domain);
elseif dim==3
    [w,dwdx]=wachspress3d_opt(x,domain);
else
   error("incompatible density field") 
end


end