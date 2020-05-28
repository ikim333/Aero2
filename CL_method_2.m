function [Cl, CL] = CL_method_2 (a_ij, b_i, c, dy, S, N)

    %% circulation
    gamma = a_ij\b_i;

    %% lift coefficient and total lift
    Cl=zeros(N,1);

    Cl = 2*gamma./c;

    CL = Cl'*dy/S;
    
end