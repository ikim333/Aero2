function [CD, CD_i] = CD_method_2 (N, Cl, cl_0, cl_alpha, alpha, twist, c, s_wing, dy, b)

    alpha_ind = zeros(N,1);
    [alpha_ind, CD_i] = induced_drag (N, Cl, cl_0, cl_alpha, alpha, twist, alpha_ind, c, dy, s_wing);

    Cdv = 0.0063 - 0.0033.*Cl + 0.0068.*Cl.^2;

    CDv = Cdv'*dy/b;

    CD = CDv + CD_i;
    