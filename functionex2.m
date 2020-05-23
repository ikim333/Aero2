function [res1, res2, res3] = functionex2(r)

    % Incident Flow
    alpha=3*(r-2);  %angle d'atac en graus
    Uinf=1;         %modul de la velocitat del flux no pertorbat
    u=[Uinf*cosd(alpha),Uinf*sind(alpha)]; %velocitat del flux no pertorbat descomposat en eles seves components

    for i=2:N
            Zc(i)=2_profile_coord(Xc(i),p,m);        %funcio de l'exercici 1
            Nc(i,1)=-z_der(Xc(i),p,m);      %funcio de l'exercici 1
            Nc(i,2)=1;
            Nc(i,:)=Nc(i,:)/norm(Nc(i,:));  %Vector normal unitari per calcular la geometria
            Zt(i,:)=thickness(Xc(i),t);
            X(i,1)=Xc(i)-Zt(i)*Nc(i,1);%punt 1, el 0,5 es perque no es l'espessor total, sino nomes la meitat superior
            X(i,2)=Zc(i)-Zt(i)*Nc(i,2);
            X(2*N-i+2,1)=Xc(i)+Zt(i)*Nc(i,1);   %Geometria del panell x
            X(2*N-i+2,2)=Zc(i)+Zt(i)*Nc(i,2);   %Geometria del panell z
    end
    
%     hold on
%     plot(X(:,1),X(:,2))
%     grid on
%     xlabel('x/c')
%     ylabel('y/c')
%     Title=['NACA: ',int2str(NACA)];
%     title(Title)
%     axis([0 1 -0.2 0.2])
    
    control=zeros(2*N,2);   %coordenades dels punts de control en coordenades globals
    l_panell=zeros(2*N,1);  %Longitud dels panells
    dx=zeros(2*N,2);        %dx i dy dels panells
    s_panell=zeros(2*N,1);  %Sinus que forma el panell
    c_panell=zeros(2*N,1);  %Cosinus que forma el panell
    n_control=zeros(2*N,2); %Vector normal al panell en el punt de control
    
    
    for i=1:(2*N)
        control(i,:)=(X(i+1,:)+X(i,:))/2;     %càlcul dels punts de control al punt mig dels panells
        
        dx(i,:)=X(i+1,:)-X(i,:);              %Càlcul de la diferencia entre cordenades entre el primer punt i el segon d'un panell
        
        l_panell(i)=sqrt(dx(i,1)^2+dx(i,2)^2);%Càlcul de la longitud dels panells
        
        s_panell(i)=-dx(i,2)/l_panell(i);     %Càlcul del sinus del panell
        c_panell(i)=dx(i,1)/l_panell(i);      %Càlcul del cosinus del panell
        
        n_control(i,:)=[s_panell(i),c_panell(i)];   %Introducció dels valors del vector normal
    end
    
    %scatter(control(:,1),control(:,2),'.')  %posa a la gràfica els punts de control
    
    control_local_panells=zeros(2*N,2*N,2); %posició en coordenades locals des de cada panell de cada punt de control
    
    for i=1:(2*N)
        for j=1:(2*N)
            control_local_panells(i,j,:)=canvi_coord(control(i,:), X(j,:),s_panell(j),c_panell(j));
        end
    end
    
    r1_2=zeros(2*N);
    r2_2=zeros(2*N);
    theta1=zeros(2*N);
    theta2=zeros(2*N);
    
    for i=1:(2*N)
        for j=1:(2*N)
            
            r1_2(i,j)=control_local_panells(i,j,1)^2+control_local_panells(i,j,2)^2;
            r2_2(i,j)=(control_local_panells(i,j,1)-l_panell(j))^2+control_local_panells(i,j,2)^2;
            
            theta1(i,j)=atan2(control_local_panells(i,j,2),control_local_panells(i,j,1));
            theta2(i,j)=atan2(control_local_panells(i,j,2),(control_local_panells(i,j,1)-l_panell(j)));
        end
    end
    
    Vind_loc=zeros(2*N,2*N,2);
    
    for i=1:(2*N)
        for j=1:(2*N)
            if j==i
                if  control_local_panells(i,j,2) < 0
                    Vind_loc(i,j,1)=-0.5*gam;
                    Vind_loc(i,j,2)=0;
                else
                    Vind_loc(i,j,1)=0.5*gam;
                    Vind_loc(i,j,2)=0;
                end
            else
                Vind_loc(i,j,:)=[(1/(2*pi))*(theta2(i,j)-theta1(i,j)), (1/(4*pi))*log(r2_2(i,j)/r1_2(i,j))];  %càlcul de la velocitat induïda en x
            end
        end
    end
    
    Vind_global=zeros(2*N,2*N,2);
    
    for i=1:(2*N)
        for j=1:(2*N)
            Vind_global(i,j,:)=[Vind_loc(i,j,1)*c_panell(j)+Vind_loc(i,j,2)*s_panell(j),-Vind_loc(i,j,1)*s_panell(j)+Vind_loc(i,j,2)*c_panell(j)];  %càlcul de la velocitat induïda en x      
        end
    end

    a=zeros(2*N);
    
    for i=1:(2*N)
        if i==N/2
            a(i,:)=0;
            a(i,1)=1;
            a(i,2*N)=1;
        else
            for j=1:2*N
                if i==j
                    a(i,j)=-0.5;
                else
                    a(i,j)=Vind_global(i,j,1)*c_panell(i)-Vind_global(i,j,2)*s_panell(i);
                end        
            end
        end
    end
    
    b=zeros(2*N,1);
    
    for i=1:(2*N)
        if i==N/2
            b(i)=0;
        else
            b(i)=-u(1)*c_panell(i)+u(2)*s_panell(i);
        end
    end
    gam=a\b;
   
    cl=zeros(2*N,1);
    for i =1:2*N
       cl(i)=gam(i)*l_panell(i);
    end
    %Cl(r)=(2/(Uinf*c))*sum(cl); %Calcul del coef.Sustentacio
    Cl(r)=(2/(Uinf*c))*sum(cl);
    
    tangencial(:,:)=[c_panell(:),-s_panell(:)]; %Vector tangencial al panell
    v_tangencial=zeros(2*N,1);
    producte_aij_gamj=zeros(2*N,1);

    for i=1:2*N
        for j=1:2*N
            if j==i
                producte_aij_gamj(j)=0;             
            else
                producte_aij_gamj(j)=a(i,j)*gam(j); %Calcul del producte dels coeficients d'influencia ger gamma del panell
            end                                     %necessari per calcular la velocitat tangencial

        end
            v_tangencial(i)=dot(u,tangencial(i,:))+sum(producte_aij_gamj)+0.5*gam(i);
    end

    cp(:)=1-(v_tangencial(:)./Uinf).^2; %Coeficient de pressió
    
    cm0=zeros(2*N,1);
    for i=1:(2*N)
        cmLE(i)=cp(i)*dot(control(i,:),dx(i,:));
    end
    
    CmLE(r)=sum(cmLE);

    q(r)=alpha; %angle d'atac de cada vegada que corre el programa
    
    %Matriu amb tots els resultats: Fila 1: alpha Fila 2:Cl Fila 3:Cm0
    resultats(1,r)=q(r);
    resultats(2,r)=Cl(r);
    resultats(3,r)=CmLE(r);  
end

