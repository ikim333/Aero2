function [Cl, CL] = CL_method_2 (a_ij, b_i, c, dy, b, N, cl_0, cl_alpha, alpha, twist, s_wing)

    %% circulation
    gamma = a_ij\b_i;
    gamma_nd = gamma./c; %non-dimensional gamma

    %% lift coefficient and total lift
    Cl=zeros(N,1);

    Cl = 2*gamma./c;

    CL = Cl'*dy/b;
    
end