%% ASSIGNMENT OF AERODYNAMICS - 29/05/2020
%AUTHORS: Miquel Badia, Daniel Longaron, Arnau Reyes

%This program computes the lift and drag distribution along a wing using
%the lifting line method

%% EXERCISE 2

% Students:
%     Miquel Badia
%     Arnau Reyes
%     Daniel Longaron

clc
clear all

%% DATA INPUT

%Perfil
NACA=2415;
m=(floor(NACA/1000))/100;        %diferencia maxima de la linia mitja respecte la corda
p=(floor((NACA/100))-m*1000)/10; %posició de la alçada maxima de la linia mitja
t=(NACA-m*100000-p*1000)/100;    %thickness
c=1;   %Corda

%N=50;   %numero divisions, el número real de panells será 2*N
gam=1; %suposem un valor inicial de gamma
rho=1;     %rho del flux
uinf=1;         %modul de la velocitat del flux no pertorbat

alpha = [-3 0 3 6 9];
                            
CL_def = zeros(1,5);
CM_LE_def = zeros(1,5);

N = 60;

save('data_ex2.mat', 'NACA', 'm', 'p', 't', 'c', 'gam', 'rho', 'uinf', 'alpha'); 

for r=1:5
        [res1, res2, cp, ctrl_points,X] = functionex2(r,alpha,N);
        CL_def(r) = res1;
        CM_LE_def(r) = res2;
end

figure;
grid on;
xlabel('\alpha');

yyaxis left;
plot(alpha(:),CL_def(:));
ylabel('C_l');
axis([-3 9 -0.2 1.2])

yyaxis right;
plot(alpha(:),CM_LE_def(:));
ylabel('C_m_L_E');

%% CP calc at fixed alpha
r = 3; % alpha = 3º
[res1, res2, cp, ctrl_points,X] = functionex2(r,alpha,N);

figure;
grid on;
xlabel('x/c');

yyaxis left;
plot(X(:,1),X(:,2));
ylabel('y/c');
axis([0 1 -0.2 0.2])
    
yyaxis right;
plot(ctrl_points(:,1),-cp(:));
ylabel('-C_p');
title("Distribució Cp");

%% Panels convergence
Nvect = [16 32 64 128 256 512 1024];

% CL and CM_LE matrix for different panels and different alpha values

CL_mat = zeros(5);
CM_LE_mat = zeros(5);

for k=1:5
    N = Nvect(k);
    for r=1:5
        [res1, res2, cp, ctrl_points,X] = functionex2(r,alpha,N);
        CL_mat(k,r) = res1;
        CM_LE_mat(k,r) = res2;
    end
end


% Error calculation at fixed Alpha 
s = size(Nvect, 2);
r = 3;
CL_vectorplot = zeros(1,s);
CMLE_vectorplot = zeros(1,s);
error_CL = zeros(1,s-1);
error_CMLE = zeros(1,s-1);

for k=1:s
    N = Nvect(k);
        [res1, res2, cp, ctrl_points,X] = functionex2(r,alpha,N);
        CL_vectorplot(k) = res1;
        CMLE_vectorplot(k) = res2;
        if k >=2
            error_CL(k-1) = abs(CL_vectorplot(k)-CL_vectorplot(k-1));
            error_CMLE(k-1) =  abs(CL_vectorplot(k)-CL_vectorplot(k-1));
        end
end

Nplot = zeros(1,s-1);
for i = 1:(s-1)
   Nplot(i) = Nvect(i); 
end

%error calcul cl
figure;
grid on;
xlabel('N (panels)');

yyaxis left;
plot(Nplot,CL_vectorplot(1:end-1));
ylabel('CL');
%axis([0 1 -0.2 0.2])
    
yyaxis right;
plot(Nplot,error_CL);
ylabel('Error CL');
title("Error CL & CL vs number of panels");

%error calcul cmle
figure;
grid on;
xlabel('N (panels)');

yyaxis left;
plot(Nplot,CMLE_vectorplot(1:end-1));
ylabel('CMLE');
%axis([0 1 -0.2 0.2])
    
yyaxis right;
plot(Nplot,error_CMLE);
ylabel('Error CMLE');
title("Error CMLE & CMLE vs number of panels");