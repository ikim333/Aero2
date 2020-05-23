%% ASSIGNMENT OF AERODYNAMICS - 29/05/2020
%AUTHORS: Miquel Badia, Daniel Longaron, Arnau Reyes

%% EXERCISE 1.2.

clc
clear all

%% DATA INPUT

main_1_1;

h0 = 0.75*c_r; %height over a flat ground
X(:,3) = h0;
Xp(:,3) = h0;

%% METHOD 2 & GROUND EFFECT

[a_ij_ge] = ground_effect (N, X, Ur, Xp, c, cl_alpha, K_v, a_ij);

[Cl_ge, CL_ge, CD_ge] = method_2 (a_ij_ge, b_i, c, dy, b, N, cl_0, cl_alpha, alpha, twist, s_wing);

CL_ge
CD_ge

%% REPRESENTATION OF RESULTS

res = 'Plot results for exercise 1.2? (1 for YES, 0 for NO):  ';
res = input(res);

if res==1

figure(2);
plot(Xp(:,2),Cl_ge(:));
grid on;
grid minor;
xlabel('y');
ylabel('Cl_ge');

else
end




