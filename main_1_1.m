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

%boundary conditions
alpha = 5 * pi/180; %angle of attack
U_inf = 1; %air speed


%% PREPROCESSING

N = 8; %number of divisions
twist = zeros(N,1); %twist definition in each section
X = zeros(N+1,3); %3 dimension case
Xp = zeros(N,3); %coordinates of each section
dy = zeros(N,1); %thickness of each section
c = zeros(N,1); %chord of each section

%to define each section, the cosine distribution will define each section
%coordinates:

[X, Xp, c, twist] = coordinates (N, b, twist_t, c_t, c_r);


%% PART 1
% lifting line method 2

% parameters to calculate circulation gamma:
















