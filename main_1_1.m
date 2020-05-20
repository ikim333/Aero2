%% ASSIGNMENT OF AERODYNAMICS - 29/05/2020
%AUTHORS: Miquel Badia, Daniel Longaron, Arnau Reyes

%This program computes the lift and drag distribution along a wing using
%the lifting line method

%% EXERCISE 1

% Students:
%     Miquel Badia
%     Arnau Reyes
%     Daniel Longaron

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

N = 164; %number of divisions
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
[gamma, gamma_nd] = circulation (N, c, U_inf, cl_0, cl_alpha, alpha, twist, Ur, X, Xp, K_v);

%lift coefficient and total lift
Cl=zeros(N,1);

Cl = 2*gamma./c;

CL = Cl'*dy/b;

% drag coefficient and total drag
alpha_ind = zeros(N,1);
[alpha_ind, CD_i] = induced_drag (N, Cl, cl_0, cl_alpha, alpha, twist, alpha_ind, c, dy, s_wing);

Cdv = 0.0063 - 0.0033.*Cl + 0.0068.*Cl.^2;

CDv = Cdv'*dy/b;

CD = CDv + CD_i

%% REPRESENTATION OF RESULTS

figure(1);
plot(Xp(:,2),Cl(:));
grid on;
grid minor;
xlabel('y');
ylabel('Cl');














