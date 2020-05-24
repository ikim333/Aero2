function [res1, res2, cp, ctrl_points, X] = functionex2(r,alpha,N)
    Xc=linspace(1,0,N+1); %considerar n+1, és la x que defineix la corda
    X=zeros(2*N+1,2);         %coordenades dels punts sobre la superficie del perfil
                              %la x son els punts en eixos globals, els reals
    
    X(1,1)=1;
    X(2*N+1,1)=1;

    zc=zeros(N+1,1);
    nc=zeros(N+1,2); %vector normal
    zt=zeros(N+1,1);
    
    load('data_ex2.mat');

    % Incident Flow
    u=[uinf*cosd(alpha(r)),uinf*sind(alpha(r))]; %velocitat del flux no pertorbat descomposat en eles seves components

    for i=2:N
            zc(i)=profile_coord_ex2(Xc(i),p,m);        %funcio de l'exercici 1
            nc(i,1)=-profile_deriv_ex2(Xc(i),p,m);      %funcio de l'exercici 1
            nc(i,2)=1;
            nc(i,:)=nc(i,:)/norm(nc(i,:));  %Vector normal unitari per calcular la geometria
            zt(i,:)=thickness_ex2(Xc(i),t);
            X(i,1)=Xc(i)-zt(i)*nc(i,1);%punt 1, el 0,5 es perque no es l'espessor total, sino nomes la meitat superior
            X(i,2)=zc(i)-zt(i)*nc(i,2);
            X(2*N-i+2,1)=Xc(i)+zt(i)*nc(i,1);   %Geometria del panell x
            X(2*N-i+2,2)=zc(i)+zt(i)*nc(i,2);   %Geometria del panell z
    end

    ctrl_points=zeros(2*N,2);   %coordenades dels punts de control en coordenades globals
    pan_long=zeros(2*N,1);  %Longitud dels panells
    dx=zeros(2*N,2);        %dx i dy dels panells
    pan_sin=zeros(2*N,1);  %Sinus que forma el panell
    pan_cos=zeros(2*N,1);  %Cosinus que forma el panell
    pan_norm=zeros(2*N,2); %Vector normal al panell en el punt de control
    
    
    for i=1:(2*N)
        ctrl_points(i,:)=(X(i+1,:)+X(i,:))/2;     %càlcul dels punts de control al punt mig dels panells
        
        dx(i,:)=X(i+1,:)-X(i,:);              %Càlcul de la diferencia entre cordenades entre el primer punt i el segon d'un panell
        
        pan_long(i)=sqrt(dx(i,1)^2+dx(i,2)^2);%Càlcul de la longitud dels panells
        
        pan_sin(i)=-dx(i,2)/pan_long(i);     %Càlcul del sinus del panell
        pan_cos(i)=dx(i,1)/pan_long(i);      %Càlcul del cosinus del panell
        
        pan_norm(i,:)=[pan_sin(i),pan_cos(i)];   %Introducció dels valors del vector normal
    end
    
    ctrl_pan_local=zeros(2*N,2*N,2); %posició en coordenades locals des de cada panell de cada punt de control
    
    for i=1:(2*N)
        for j=1:(2*N)
            ctrl_pan_local(i,j,:)=coord_change_ex2(ctrl_points(i,:), X(j,:),pan_sin(j),pan_cos(j));
        end
    end
    
    r1_2=zeros(2*N);
    r2_2=zeros(2*N);
    theta1=zeros(2*N);
    theta2=zeros(2*N);
    
    for i=1:(2*N)
        for j=1:(2*N)
            
            r1_2(i,j)=ctrl_pan_local(i,j,1)^2+ctrl_pan_local(i,j,2)^2;
            r2_2(i,j)=(ctrl_pan_local(i,j,1)-pan_long(j))^2+ctrl_pan_local(i,j,2)^2;
            
            theta1(i,j)=atan2(ctrl_pan_local(i,j,2),ctrl_pan_local(i,j,1));
            theta2(i,j)=atan2(ctrl_pan_local(i,j,2),(ctrl_pan_local(i,j,1)-pan_long(j)));
        end
    end
    
    vind_local=zeros(2*N,2*N,2);
    
    for i=1:(2*N)
        for j=1:(2*N)
            if j==i
                if  ctrl_pan_local(i,j,2) < 0
                    vind_local(i,j,1)=-0.5*gam;
                    vind_local(i,j,2)=0;
                else
                    vind_local(i,j,1)=0.5*gam;
                    vind_local(i,j,2)=0;
                end
            else
                vind_local(i,j,:)=[(1/(2*pi))*(theta2(i,j)-theta1(i,j)), (1/(4*pi))*log(r2_2(i,j)/r1_2(i,j))];  %càlcul de la velocitat induïda en x
            end
        end
    end
    
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
    gam=A\B;
   
    cl=zeros(2*N,1);
    for i =1:2*N
       cl(i)=gam(i)*pan_long(i);
    end
    %Cl(r)=(2/(uinf*c))*sum(cl); %Calcul del coef.Sustentacio
    Cl(r)=(2/(uinf*c))*sum(cl);
    
    pan_tg_vect(:,:)=[pan_cos(:),-pan_sin(:)]; %Vector tangencial al panell
    v_tg=zeros(2*N,1);
    aij_gamj_prod=zeros(2*N,1);

    for i=1:2*N
        for j=1:2*N
            if j==i
                aij_gamj_prod(j)=0;             
            else
                aij_gamj_prod(j)=A(i,j)*gam(j); %Calcul del producte dels coeficients d'influencia ger gamma del panell
            end                                     %necessari per calcular la velocitat tangencial

        end
            v_tg(i)=dot(u,pan_tg_vect(i,:))+sum(aij_gamj_prod)+0.5*gam(i);
    end

    cp(:)=1-(v_tg(:)./uinf).^2; %Coeficient de pressió
    
    cm0=zeros(2*N,1);
    for i=1:(2*N)
        cmLE(i)=cp(i)*dot(ctrl_points(i,:),dx(i,:));
    end
    
    CmLE(r)=sum(cmLE);
    
    res1=Cl(r);
    res2=CmLE(r);  
end

