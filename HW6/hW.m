clear all; clc; close all

A = 10^6*[96, -24; -24, 8];
B = [5; -.625];
coefs = A\B

x = 0:.005:.5;
L = .5;
v = coefs(1)*(3*(x/L).^2 - 2*(x/L).^3) + L*coefs(2)*((x/L).^3 - (x/L).^2);
figure(); plot(x,v)

v_prime = coefs(1)*(6*x/(L^2) - 6*(x.^2)/(L^3)) + coefs(2)*L*(3*(x.^2)/(L^3) - 2*x/(L^2));
figure(); plot(x,v_prime)

EI = 10^6;
v_double_prime = coefs(1)*(6/(L^2) - 12*(x)/(L^3)) + coefs(2)*(6*(x)/(L^2) - 2/(L));
M = EI*v_double_prime;
figure(); plot(x,M)