%% ASSIGNMENT OF AERODYNAMICS - 29/05/2020
%AUTHORS: Miquel Badia, Daniel Longaron, Arnau Reyes

%% EXERCISE 1.4.

clc
clear all

%% DATA INPUT

% exercise 1.1 is computed to make changes on it and develop 1.4
main_1_1;

% alpha_ind will be needed, so the function is re-called:
alpha_ind = zeros(N,1);
[alpha_ind, CD_i] = induced_drag (N, Cl, cl_0, cl_alpha, alpha, twist, alpha_ind, c, dy, s_wing);

alpha = 5*pi/180;
delta_left = -5*pi/180;
delta_right = 5*pi/180;

% those sections affected by the aileron will have a new Cl:
Cl_ail = Cl;
for i=1:N
    %Left aileron
    if Xp(i,2) <= -b/2 + b_a
        Cl_ail(i) = 0.2 + 6.3*(alpha+alpha_ind(i)) + 4.1*delta_left;
    %Right aileron
    elseif Xp(i,2) >= b/2 - b_a
        Cl_ail(i) = 0.2 + 6.3*(alpha+alpha_ind(i)) + 4.1*delta_right;
    else
    end
    
end

CL_ail = Cl_ail'*dy/b;


%Compute roll moment

Roll_moment = 0.5.*dy.*c.*Cl_ail.*Xp(:,2);

CR=sum(Roll_moment)/(b*s_wing);

% drag with ailerons computed next
[CD_ail, CD_i_ail] = CD_method_2 (N, Cl_ail, cl_0, cl_alpha, alpha, twist, c, s_wing, dy, b);

%% RESULTS

disp('Results por part 1.4')
CL_ail
CD_i_ail
CD_ail
CR

res = 'Plot results for exercise 1.4? (1 for YES, 0 for NO):  ';
res = input(res);

if res==1

figure(3);
plot(Xp(:,2),Cl_ail(:));
grid on;
grid minor;
title('Lift coefficient distribution');
xlabel('y');
ylabel('C_l');

figure(4);
plot(Xp(:,2),Cr(:));
grid on;
grid minor;
title('Rolling moment coefficient distribution');
xlabel('y');
ylabel('C_r');

else
end