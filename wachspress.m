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

function [w,dwdx]=wachspress(x,domain)
% [w,dwdx]=wachspress(x,domain)
% INPUTS :
%   - x : the vector of the cartesian coordinates
%   - domain : interpolation domain
% OUTPUTS :
%   - w : Wachspress shape functions
%   - dwdx : gradients of w

% Use wachspress2d_opt and wachspress3d_opt, optimized from :
% M. Floater, A. Gillette, and N. Sukumar,
% “Gradient bounds for Wachspress coordinates on polytopes,”
% SIAM J. Numer. Anal., vol. 52, no. 1, pp. 515–532, 2014,%
% doi: 10.1137/130925712

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

function [phi,dphi] = wachspress2d_opt(x,domain)
%% Evaluate Wachspress basis functions and their gradients in a convex polygon
%% Inputs:
% v    : [x1 y1; x2 y2; ...; xn yn], the n vertices of the polygon in ccw
% x    : [x(1) x(2)], the point at which the basis functions are computed
%% Outputs:
% phi  : output basis functions = [phi_1; ...; phi_n]
% dphi : output gradient of basis functions = [dphi_1; ...; dphi_n]

% based on the Matlab code from :
% M. Floater, A. Gillette, and N. Sukumar,
% “Gradient bounds for Wachspress coordinates on polytopes,”
% SIAM J. Numer. Anal., vol. 52, no. 1, pp. 515–532, 2014,%
% doi: 10.1137/130925712.

% vectorized by T. Cherriere (2021)

v=domain.vertices;
np=size(x,1);
n = size(v,1);

un = reshape(getNormals(v).',2,1,[]);
h = mult(reshape(v.',1,2,[]) - reshape(x.',1,2,1,np),un);
p=un./h;
i=1:n; im1=mod(i-2,n) + 1;
w=p(1,:,im1,:).*p(2,:,i,:)-p(1,:,i,:).*p(2,:,im1,:);
wsum=sum(w,3);
R= permute(p(:,:,im1,:)+p(:,:,i,:),[3,1,4,2]);
phi = t(shiftdim(w./wsum,1));
phiR = mult(t(phi),R);
dphi = phi.*(R-phiR);

end

function un = getNormals(v)
% Function to compute the outward unit normal to each edge
n = size(v,1);un = zeros(n,2);
ind1=mod(1:n,n)+1;
ind2=1:n;
d = v(ind1,:) - v(ind2,:);
un(ind2,:) = [d(ind2,2) -d(ind2,1)]./vecnorm(d,2,2);
end

function [phi, dphi] = wachspress3d_opt(x,domain)
% Evaluate Wachspress basis functions and their gradients
% in a convex polyhedron
%% Inputs:
% polyhedra. ...
% v    : [x1 y1 z1; x2 y2 z2; ...; xn yn zn], the n vertices of the polyhedron
% g    : cell array: [i1 i2 ... i_{k1}]; . . . ; [i1 i2 ... i_{kn}],
%        which are the n neighborhood graphs, each in some counter-clockwise order
%        as seen from outside the polyhedron
% un   : [x1 y1 z1; x2 y2 z2; ...; xm ym zm], unit normal to each facet
% x    : [x(1) x(2) x(3)], the point at which the basis functions are computed
% Outputs:% phi  : basis functions = [phi_1; ...; phi_n]
% dphi : gradient of basis functions = [dphi_1; ...; dphi_n]

% based on the Matlab code from :
% M. Floater, A. Gillette, and N. Sukumar, 
% “Gradient bounds for Wachspress coordinates on polytopes,”
% SIAM J. Numer. Anal., vol. 52, no. 1, pp. 515–532, 2014,
% doi: 10.1137/130925712.

% vectorized by T. Cherriere (2021)

np=size(x,1);
v=domain.vertices; un=domain.normals; g=domain.vertices2facets;

n = size(v,1);w = zeros(n,1,np); R = zeros(n,3,np);

x=reshape(x.',1,3,1,np);
un=reshape(un.',1,3,[]);
v=reshape(v.',1,3,[]);

for i = 1:n
    f = g{i};
    k = length(f);
    h = mult(v(:,:,i) - x,t(un(:,:,f)));
    p = permute(un(:,:,f)./ h,[1,2,4,3]);
    j = 1:k-2;
    wloc=permute(p(1,1,:,j).*(p(1,2,:,j+1).*p(1,3,:,k)-p(1,2,:,k).*p(1,3,:,j+1))+...
         p(1,1,:,j+1).*(p(1,3,:,j).*p(1,2,:,k)-p(1,3,:,k).*p(1,2,:,j))+...
        p(1,1,:,k).*(p(1,2,:,j).*p(1,3,:,j+1)-p(1,2,:,j+1).*p(1,3,:,j)),[4,2,3,1]);
    Rloc=permute(p(:,:,:,j)+p(:,:,:,j+1)+p(:,:,:,k),[4,2,3,1]); 
    w(i,:,:) = sum(wloc,1);
    R(i,:,:) = mult(t(wloc),Rloc)./ w(i,:,:);
end

wsum = sum(w,1);
phi  = w./wsum;
phiR = mult(t(phi),R);
dphi= phi .* (R - phiR);

end
