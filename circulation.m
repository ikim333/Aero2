%% FUNCTION CIRCULATION

% this function computes the circulation at each point of the wing
% throughout method 2 with the influent coefficient a and the independent
% term b, calculated next (TEMA 4, DIAPO 17 & 23)

function [gamma, gamma_nd] = circulation (N, c, U_inf, cl_0, cl_alpha, alpha, twist, Ur, X, Xp, K_v)

a = zeros(N,N);
b = zeros(N,1);


    for i=1:N
        
        b(i) = 0.5*c(i)*U_inf*(cl_0+cl_alpha*(alpha+twist(i)));
        
        for j=1:N
            
            if i==j
                V_l = semi_infinite_vortex(Ur,X(j,:),Xp(i,:));
                V_r = semi_infinite_vortex(Ur,X(j+1,:),Xp(i,:));
                V = V_l - V_r;
                a(i,j) = -0.5*c(i)*cl_alpha*dot(V,K_v)+1;
            
            else
                V = horseshoe(Ur,X(j,:),X(j+1,:),Xp(i,:));
                a(i,j) = -0.5*c(i)*cl_alpha*dot(V,K_v);
            
            end
        end
    end
    
    gamma = a\b;
    gamma_nd = gamma./c; %non-dimensional gamma

end