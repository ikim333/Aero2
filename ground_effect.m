function [a_ij_ge] = ground_effect (N, X, Ur, Xp, c, cl_alpha, K_v, a_ij)

a_ij_ge = zeros(N,N);

    for i=1:N
        for j=1:N

        xim_j = X(j,:);    
        xim_j1 = X(j+1,:);  
        xim_j(3) = -X(j,3);
        xim_j1(3) = -X(j+1,3);

        V_prima = horseshoe(Ur,xim_j,xim_j1,Xp(i,:));

        a_prima = -0.5*c(i)*cl_alpha*dot(V_prima,K_v);
        
        a_ij_ge(i,j) = a_ij(i,j) - a_prima;
       
        end
    end
    
    
    
end