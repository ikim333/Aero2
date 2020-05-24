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

CL_mat = zeros(5);
CM_LE_mat = zeros(5);

Nvect=[16 32 64 128 256];

save('data_ex2.mat', 'NACA', 'm', 'p', 't', 'c', 'gam', 'rho', 'uinf', 'alpha'); 

for k=1:5
    N = Nvect(k);
    for r=1:5
        [res1, res2] = functionex2(r,alpha,N);
        CL_mat(k,r) = res1;
        CM_LE_mat(k,r) = res2;
    end
    
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