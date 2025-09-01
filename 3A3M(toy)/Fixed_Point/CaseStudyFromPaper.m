close all;
clear;
clc
addpath('/Users/kiguli/Documents/My_papers/YALMIP-master')
addpath('/Users/kiguli/Documents/mosek/9.3/toolbox/r2015a')
%% Original linear system S_1
% x1dot = A1x1 + B1u1 + Bdw
% y1 = C1x1
A1 = [zeros(2) eye(2) zeros(2); zeros(2) zeros(2) eye(2); zeros(2) zeros(2) zeros(2)]; % n1xn1
B1 = [zeros(2);zeros(2);eye(2)]; % n1xp1
C1 = [eye(2) zeros(2) zeros(2)]; % mxn1
%% Abstract System S_2 (no disturbance)
% x2dot = A2x2 + B2u2
% y2 = C2x2
A2 = zeros(2); % n2xn2
B2 = eye(2); % n2xp2
C2 = eye(2); % mxn2
%assume n2 <= n1
%% Approximate Simulation Variables
K = -[52*eye(2) 52.3*eye(2) 13*eye(2)];
P = [eye(2) zeros(2) zeros(2)].';
Q = zeros(2);
lambda =1.1;

%% Calculate M and K using semi-definite programming
M = C1*C1.'; %positive definite symmetric matrix
Mbar = inv(M);
% solve using:
%Mbar = inv(M);
%Kbar = K*Mbar;
% [Mbar Mbar*C1.'; C1*Mbar Identity] >= 0
% Mbar*A1.' + A1*Mbar + Kbar.'*B1.'+B1*Kbar <= -2*lambda*Mbar
% Identity is mxm

Mat = [Mbar Mbar*C1.'; C1*Mbar eye(3,3)];
Kbar = sdpvar(1,30);
F = [Mat >= 0, Mbar*A1.' + A1*Mbar + Kbar.'*B1.'+B1*Kbar <= -2*lambda*Mbar];
optimize(F);

%% Simulation Function
% V(x1,x2) = norm(C1x1(0) - C2x2(0))
V = norm(y1-y2);
gamma_1 = norm(sqrt(M)*Bd)*norm(d)/lambda;
gamma_2 = norm(sqrt(M)*(B1*R - P*B2))*norm(u2)/lambda;
% as long as:
% gamma_1 
% + gamma_2
% < V(x1,x2)

%% Bounded disturbance
Bd = [-0.2 -0.2 0 0 0 0].';
dist = 1;

