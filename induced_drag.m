%% FUNCTION INDUCED DRAG
% this function computes the induced drag (TEMA 4, DIAPO 21)

function [alpha_ind, CD_i] = induced_drag (N, Cl, cl_0, cl_alpha, alpha, twist, alpha_ind, c, dy, s_wing)


    CD_i = 0; %valor on comença el sumatori

    for i=1:N
        alpha_ind(i) = (Cl(i)-cl_0)/cl_alpha - alpha - twist(i);
        CD_i = CD_i - Cl(i)*alpha_ind(i)*c(i)*dy(i)/s_wing;
    end    
    
end