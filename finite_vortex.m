%% FUNCTION FINITE VORTEX

% this function computes the induced velocity at point P with given P_1 and
% P_2, which create a finite vortex line with this points at its extremes.

function V = finite_vortex(X1,X2,Xp)

r1=Xp-X1;
r2=Xp-X2;

n_r1=norm(r1);
n_r2=norm(r2);

V =cross(1/(4*pi)*(n_r1+n_r2)/(n_r1*n_r2*(n_r1*n_r2+dot(r1,r2)))*r1 ,r2);

end