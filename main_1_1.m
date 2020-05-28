%% ASSIGNMENT OF AERODYNAMICS - 29/05/2020
%AUTHORS: Miquel Badia, Daniel Longaron, Arnau Reyes

%% EXERCISE 1.1.

clc
clear all

%% DATA INPUT

% Wing:
twist_t = -3 *(pi/180); %twist angle, rad
c_r = 1.2; %chord at root
c_t = 0.8; %chord at tip
b = 6; 
b_a = 1.8; 
cl_alpha = 6.3;
cl_0 = 0.2;
s_wing = 2*(0.5*b*(c_r+c_t)/2); %geometry

%boundary conditions
alpha = 5 * pi/180; %angle of attack
U_inf = 1; %air speed

Ur = [-1,0,0]; %vector that points U_inf_v's direction
K_v = [0,0,1];


%% PREPROCESSING

N = 64; %number of divisions
twist = zeros(N,1); %twist definition in each section
X = zeros(N+1,3); %3 dimension case
Xp = zeros(N,3); %coordinates of each section
dy = zeros(N,1); %thickness of each section
c = zeros(N,1); %chord of each section



%to define each section, the cosine distribution will define each section
%coordinates:

[X, Xp, c, twist, dy] = coordinates (N, b, twist_t, c_t, c_r);


%% PART 1
% lifting line method 2

% circulation
[a_ij, b_i] = coefficients (N, c, U_inf, cl_0, cl_alpha, alpha, twist, Ur, X, Xp, K_v);

[Cl, CL] = CL_method_2 (a_ij, b_i, c, dy, s_wing, N);
[CD, CD_i, CDv] = CD_method_2 (N, Cl, cl_0, cl_alpha, alpha, twist, c, s_wing, dy, b);


%% RESULTS

disp('Results por part 1.1')
CL
CD_i
CDv
CD

res = 'Plot results for exercise 1.1? (1 for YES, 0 for NO):  ';
res = input(res);

if res==1

figure(1);
plot(Xp(:,2),Cl(:));
grid on;
grid minor;
title('Lift coefficient distribution');
xlabel('y');
ylabel('Cl');

else
end














