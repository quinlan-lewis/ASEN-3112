function [wn, f, T, modes, ortho_M, Mg, Kg] = MDOF_calc(M, K)

syms lambda u1 u2
tot = K - lambda*M;
u = [1; u2];
eq1 = det(tot) == 0

%solve for the natural frequencies using vpasolve
wn = double(sqrt(vpasolve(eq1,lambda)));
f = wn/(2*pi);
T = 1/f;

%solving for modes based on natural frequencies and setting one value of
%mode to 1
eq2 = tot*u;
eq_solve1 = subs(eq2(1),lambda,(wn(1)^2)) == 0;
u12 = vpasolve(eq_solve1,u2);
eq_solve2 = subs(eq2(1),lambda,(wn(2)^2)) == 0;
u22 = vpasolve(eq_solve2,u2);

mode1 = double([1; u12]);
mode2 = double([1; u22]);
mode1 = mode1./(max(abs(mode1)));
mode2 = mode2./(max(abs(mode2)));
modes = double([mode1, mode2]);

%orthonormalize to mass
M1 = mode1'*M*mode1;
M2 = mode2'*M*mode2;
mode1_prime = (1/sqrt(M1))*mode1;
mode2_prime = (1/sqrt(M2))*mode2;
ortho_M = double([mode1_prime, mode2_prime]);

%calculating Mg and Kg
Mg = double(ortho_M'*M*ortho_M);
Kg = double(ortho_M'*K*ortho_M);

end

