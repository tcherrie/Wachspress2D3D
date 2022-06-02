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