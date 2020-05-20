%% FUNCTION HORSESHOE
% aquesta funci� calcula la velocitat indu�da en un punt en funci� de dos
% punts (que seran consecutius), formada per un v�rtex finit i dos infinits
% (TEMA4, DIAPO 10)

function V = horseshoe(Ur,X1,X2,Xp)

    V_left = semi_infinite_vortex(Ur,X1,Xp);
    V_top = finite_vortex(X1,X2,Xp);
    V_right = semi_infinite_vortex(Ur,X2,Xp);

    V = V_left + V_top - V_right;

end