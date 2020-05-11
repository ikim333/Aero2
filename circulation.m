function circulation


    for i=1:N
        
        b(i) = 0.5*c(i)*U_inf*(cl_0+cl_alpha*(alpha+twist(i)));
        
        for j=1:N_div
            
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

end