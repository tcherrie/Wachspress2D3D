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