clear; clc; close all;

K = [1000, -1000; -1000, 2000]; M = [2, 0; 0, 3];

syms wn m k
tot = K - (wn)*M;

eq1 = det(tot) == 0
nat_freq = sqrt(vpasolve(eq1,wn))
nat_fref_f = nat_freq*(1/(2*pi))
period = (1./nat_fref_f)

u12 = (1000-(nat_freq(1)^2)*2)/1000
u11 = 1000*u12/(1000-(nat_freq(1)^2)*2)
u1 = [u11; u12];

u22 = (1000-(nat_freq(2)^2)*2)/1000
u21 = 1000*u22/(1000-(nat_freq(2)^2)*2)
u2 = [u21; u22];

M1 = u1'*M*m*u1;
M2 = u2'*M*m*u2;

K1 = u1'*K*k*u1;
K2 = u2'*K*k*u2;

U1_scaled = u1/(sqrt(10/3));
U2_scaled = u2/(sqrt(5));

% RHS = [0; 0];
% LHS = K - (wn^2)*M
% LHS = double(subs(LHS,wn,nat_freq(1)))
% U = LHS\RHS