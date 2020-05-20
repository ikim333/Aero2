%% FUNCTION CIRCULATION
% aquetsa funció calcula la circulació en cada punt de l'ala estudiat
% mitjançant el mètode 2 a partir del coeficient d'influència a i el terme
% independent b, calculats a continuació (TEMA 4, DIAPOS 17 i 23)

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