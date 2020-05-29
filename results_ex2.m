%% CL, CMLE, CP, XCP and airfoil data computation

function [res1, res2, res3, cp, ctrl_points, X] = results_ex2(r,alpha,N)
    %% Initial vectors definition
    
    Xc=linspace(1,0,N+1);   % Chord definition with n+1 nodes
    X=zeros(2*N+1,2);       % Coordinates from the points over the profile surface; X : global axis points  
    X(1,1)=1;
    X(2*N+1,1)=1;

    zc=zeros(N+1,1);
    nc=zeros(N+1,2);  
    zt=zeros(N+1,1);
    
    load('data_ex2.mat');   % We load the input data from the main program

    % Incident Flow
    u=[uinf*cosd(alpha(r)),uinf*sind(alpha(r))]; % Non disturbed flow velocity in both components

    %% Section 1
    for i=2:N
            zc(i)=profile_coord_ex2(Xc(i),p,m);     % Calculates the mean curvatrue line from a NACA profile.       
            nc(i,1)=-profile_deriv_ex2(Xc(i),p,m);  % Calculates the perpendicular to the previously calculated mean curvature line.     
            nc(i,2)=1;
            nc(i,:)=nc(i,:)/norm(nc(i,:));          % Normal unit vector to calculate geometry.
            zt(i,:)=thickness_ex2(Xc(i),t);         % Calculates the profile thickness with a given empiric formula.
            X(i,1)=Xc(i)-zt(i)*nc(i,1);             % Profile boundaries
            X(i,2)=zc(i)-zt(i)*nc(i,2);
            X(2*N-i+2,1)=Xc(i)+zt(i)*nc(i,1);       % X Panel geometry
            X(2*N-i+2,2)=zc(i)+zt(i)*nc(i,2);       % Z Panel geometry
    end
    
    %% Section 2
    
    ctrl_points=zeros(2*N,2);   % Control points coordinates in globals
    pan_len=zeros(2*N,1);       % Panels lenght
    dx=zeros(2*N,2);            % Panels dx and dy
    pan_sin=zeros(2*N,1);       % Panels sinus
    pan_cos=zeros(2*N,1);       % Panels cosinus
    pan_norm=zeros(2*N,2);      % Perpendicular vector to the panel
    
    for i=1:(2*N)
        ctrl_points(i,:)=(X(i+1,:)+X(i,:))/2;   % Control points calculation in the middle of the panel
        
        dx(i,:)=X(i+1,:)-X(i,:);                % Coordinates difference between first and second point from a panel
        
        pan_len(i)=sqrt(dx(i,1)^2+dx(i,2)^2);   % Panels lenght calculation
        
        pan_sin(i)=-dx(i,2)/pan_len(i);         % Panel sinus calculation
        pan_cos(i)=dx(i,1)/pan_len(i);          % Panel cosinus calculation
        
        pan_norm(i,:)=[pan_sin(i),pan_cos(i)];  % Sinus and Cosinus computation to the panel normal vector
    end
    
    ctrl_pan_local=zeros(2*N,2*N,2);            % Position in local coordinates from each panel from each control point
    
    for i=1:(2*N)
        for j=1:(2*N)
            ctrl_pan_local(i,j,:)=coord_change_ex2(ctrl_points(i,:), X(j,:),pan_sin(j),pan_cos(j));
        end
    end
    
    %% Section 3
    
    r1_2=zeros(2*N);
    r2_2=zeros(2*N);
    theta1=zeros(2*N);
    theta2=zeros(2*N);
    
    for i=1:(2*N)
        for j=1:(2*N)
            
            r1_2(i,j)=ctrl_pan_local(i,j,1)^2+ctrl_pan_local(i,j,2)^2;
            r2_2(i,j)=(ctrl_pan_local(i,j,1)-pan_len(j))^2+ctrl_pan_local(i,j,2)^2;
            
            theta1(i,j)=atan2(ctrl_pan_local(i,j,2),ctrl_pan_local(i,j,1));
            theta2(i,j)=atan2(ctrl_pan_local(i,j,2),(ctrl_pan_local(i,j,1)-pan_len(j)));
        end
    end
    
    vind_local=zeros(2*N,2*N,2);
    
    for i=1:(2*N)
        for j=1:(2*N)
            if j==i
                if  ctrl_pan_local(i,j,2) < 0
                    vind_local(i,j,1)=-0.5*gamma;
                    vind_local(i,j,2)=0;
                else
                    vind_local(i,j,1)=0.5*gamma;
                    vind_local(i,j,2)=0;
                end
            else
                vind_local(i,j,:)=[(1/(2*pi))*(theta2(i,j)-theta1(i,j)), (1/(4*pi))*log(r2_2(i,j)/r1_2(i,j))];  %càlcul de la velocitat induïda en x
            end
        end
    end
    
    %% Section 4
    vind_global=zeros(2*N,2*N,2);
    
    for i=1:(2*N)
        for j=1:(2*N)
            vind_global(i,j,:)=[vind_local(i,j,1)*pan_cos(j)+vind_local(i,j,2)*pan_sin(j),-vind_local(i,j,1)*pan_sin(j)+vind_local(i,j,2)*pan_cos(j)];  %càlcul de la velocitat induïda en x      
        end
    end
    
    A=zeros(2*N);
    
    for i=1:(2*N)
        if i==N/2
            A(i,:)=0;
            A(i,1)=1;
            A(i,2*N)=1;
        else
            for j=1:2*N
                if i==j
                    A(i,j)=-0.5;
                else
                    A(i,j)=vind_global(i,j,1)*pan_cos(i)-vind_global(i,j,2)*pan_sin(i);
                end        
            end
        end
    end
    
    B=zeros(2*N,1);
    
    for i=1:(2*N)
        if i==N/2
            B(i)=0;
        else
            B(i)=-u(1)*pan_cos(i)+u(2)*pan_sin(i);
        end
    end
    gamma=A\B;
   
    %% Lift Coefficient Calculation
    cl=zeros(2*N,1);
    for i =1:2*N
       cl(i)=gamma(i)*pan_len(i);
    end
    
    Cl(r)=(2/(uinf*c))*sum(cl);
    
    %% Pressure Coefficient Calculation
    pan_tg_vect(:,:)=[pan_cos(:),-pan_sin(:)]; % Tangent vector to the panel
    v_tg=zeros(2*N,1);
    aij_gammaj_prod=zeros(2*N,1);

    for i=1:2*N
        for j=1:2*N
            if j==i
                aij_gammaj_prod(j)=0;             
            else
                aij_gammaj_prod(j)=A(i,j)*gamma(j); % Here we calculate the product from the influence coefficients 
            end                                     % at the panel gamma in order to obtain the tangencial velocity.

        end
            v_tg(i)=dot(u,pan_tg_vect(i,:))+sum(aij_gammaj_prod)+0.5*gamma(i);
    end

    cp(:)=1-(v_tg(:)./uinf).^2;
    
    %% Momentum Coefficient at Leading Edge calculation
    
    cm0=zeros(2*N,1);
    for i=1:(2*N)
        cmLE(i)=cp(i)*dot(ctrl_points(i,:),dx(i,:));
    end
    
    CmLE(r)=sum(cmLE);
    xcp = abs(CmLE/Cl);
    
    res1=Cl(r);
    res2=CmLE(r);
    res3=xcp;
    
    
end

