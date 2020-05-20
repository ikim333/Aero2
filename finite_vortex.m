%%FUNCTION FINITE VORTEX
% aquesta funció calcula la velocitat induïda en el punt P donats P_1 i
% P_2, els quals formen una linia de vòrtex finita amb aquests punts en els
% seus extrems

function V = finite_vortex(X1,X2,Xp)

r1=Xp-X1;
r2=Xp-X2;

n_r1=norm(r1);
n_r2=norm(r2);

V =cross(1/(4*pi)*(n_r1+n_r2)/(n_r1*n_r2*(n_r1*n_r2+dot(r1,r2)))*r1 ,r2);

end