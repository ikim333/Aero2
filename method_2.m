function [Cl, CL, CD] = method_2 (a_ij, b_i, c, dy, b, N, cl_0, cl_alpha, alpha, twist, s_wing)

    %% circulation
    gamma = a_ij\b_i;
    gamma_nd = gamma./c; %non-dimensional gamma

    %% lift coefficient and total lift
    Cl=zeros(N,1);

    Cl = 2*gamma./c;

    CL = Cl'*dy/b;

    %% drag coefficient and total drag
    alpha_ind = zeros(N,1);
    [alpha_ind, CD_i] = induced_drag (N, Cl, cl_0, cl_alpha, alpha, twist, alpha_ind, c, dy, s_wing);

    Cdv = 0.0063 - 0.0033.*Cl + 0.0068.*Cl.^2;

    CDv = Cdv'*dy/b;

    CD = CDv + CD_i;
    
end