clear; clc; close all;

fileList = dir(fullfile('OTW','S8-DAQ-OTW-G15G16.txt'));

N       = length(fileList); % NUMBER OF FILES
data    = cell(N,1); % INIT
L = 10; %length [in]
t = 1/16; %thickness [in]
G = 3.75*10^6; %[psi]
beta = 1/3;
r_e = 3/8; %[in]
r_i = r_e - t;
r_avg = (r_e + r_i)/2;

for i = 1:N
    data{i} = readtable(fullfile('OTW',fileList(i).name));
    data{i}.Properties.VariableNames = {'t_s','rot_def','shear_strain','Torque_inlbf','Axial_in'};
end

shear_strain_theory = t*(pi/180*data{:,:}.rot_def)/L;
shear_strain_exp = pi/180*(data{:,:}.shear_strain);
Torque = -1*data{:,:}.Torque_inlbf;

figure()
hold on; grid on;
scatter(shear_strain_exp,Torque)
scatter(shear_strain_theory,Torque)
title('Torque vs. Shear Strain')
ylabel('Negative Torque [inlbf]'); xlabel('Shear Strain [rad]'); ylim([-25 5]); xlim([-.002 .006])

%GJ for experimental and theoretical
[p_exp,s1] = polyfit(shear_strain_exp,Torque,1);
[p_theory,s2] = polyfit(shear_strain_theory,Torque,1);

GJ_exp = p_exp(1)*t
GJ_theory = p_theory(1)*t

J_beta = beta*2*pi*r_avg*(t^3);
GJ_pure_theory = J_beta*G

[tmp_exp,delta1] = polyval(p_exp,shear_strain_exp,s1);
plot(shear_strain_exp,tmp_exp,'LineWidth',2);
[tmp_theory,delta2] = polyval(p_theory,shear_strain_theory,s2);
plot(shear_strain_theory,tmp_theory,'k','LineWidth',2);
legend('Experimental Data','Open Thin Wall Theory','Least Square Fit for Experimental Data','Least Square Fit for Theory','location','South')

%Uncertainty Analysis
un_exp = ((0.3343 - -20.72)/(0.0003924 - -0.001363))*r_avg;
un_theory = ((0.3343 - -20.72)/(0.005285 - 0.005146))*r_avg;

%Uncertainty for Torque
% error_exp = zeros(length(shear_strain_exp),1);
% error_theory = zeros(length(shear_strain_exp),1);
% for i = 1:length(shear_strain_exp)    
%     error_exp(i) = (Torque(i) - tmp_exp(i))/tmp_exp(i);
%     error_theory(i) = (Torque(i) - tmp_theory(i))/tmp_theory(i);
% end
% error_exp_avg = mean(error_exp)
% error_theory_avg = mean(error_theory)