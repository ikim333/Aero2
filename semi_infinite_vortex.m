%% FUNCTION SEMI INFINITE VORTEX

% this function computes the velocity at point P if point P_1 tends to
% infinite (TEMA 4, DIAPO 9)

function V = semi_infinite_vortex(Ur,X2,Xp)

    r2 = Xp-X2;

    n_r2 = norm(r2);

    u_r2 = r2/n_r2;
 
    V = cross(1/(4*pi)*((1-dot(Ur,u_r2))/(norm(cross(Ur,r2)))^2)*Ur , r2); 

end

