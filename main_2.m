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
c=1;                             %Corda

N=50;   %numero divisions, el número real de panells será 2*N
gam=1; %suposem un valor inicial de gamma
rho=1;     %rho del flux

    
Xc=linspace(1,0,N+1); %considerar n+1, és la x que defineix la corda
X=zeros(2*N+1,2);         %coordenades dels punts sobre la superficie del perfil
                              %la x son els punts en eixos globals, els reals

Zc=zeros(N+1,1);
Nc=zeros(N+1,2); %vector normal
Zt=zeros(N+1,1);

X(1,1)=1;
X(2*N+1,1)=1;

for r=1:5
    [res1, res2, res3] = functionex2(r);
end
resultats

figure;
grid on;
xlabel('\alpha');

yyaxis left;
plot(q(:),Cl(:));
ylabel('C_l');
axis([-3 9 -0.2 1.2])

yyaxis right;
plot(q(:),CmLE(:));
ylabel('C_m_L_E');