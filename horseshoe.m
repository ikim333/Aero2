%% FUNCTION HORSESHOE

% this function computes the induced velocity at a point as the function of
% 2 consecutive points, it is formed by a finite vortex and two infinte
% vortex (TEMA 4, DIAPO10)

function V = horseshoe(Ur,X1,X2,Xp)

    V_left = semi_infinite_vortex(Ur,X1,Xp);
    V_top = finite_vortex(X1,X2,Xp);
    V_right = semi_infinite_vortex(Ur,X2,Xp);

    V = V_left + V_top - V_right;

end