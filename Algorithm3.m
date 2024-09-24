clc;
clear;
close all;

% alpha_f -----> [0,30,70,85]

sigma_max = deg2rad(75);

d_i = 4;
alpha_f = deg2rad(0);   %Variable
theta_i = deg2rad(67.8);
alpha_i = sigma_max + theta_i;
% d_fdes = d_i;
d_fdes = 5;
sigma_i = sigma_max;
del_sigma_f = 0.0001;
epsilon = 0.005;
sigma_f1 = sigma_max;
theta_f = sigma_f1-alpha_f;

sigma_f2 = sigma_f(d_i,d_fdes,sigma_max,alpha_i,theta_i,alpha_f,del_sigma_f,epsilon,sigma_i,sigma_f1,theta_f);
sigma_f2 = rad2deg(sigma_f2);

% Calculation of sigma_f
function sigma_ff = sigma_f(d_i,d_fdes,sigma_max,alpha_i,theta_i,alpha_f,del_sigma_f,epsilon,sigma_i,sigma_f1,theta_f)


N = (alpha_i - alpha_f) / (theta_i - alpha_f + sigma_f1);
d_f_ = d_f(d_i, theta_i, sigma_i, sigma_f1, N, theta_f);

if abs(d_f_ - d_fdes) < epsilon
    sigma_ff = sigma_f1;
else
    sigma_f1 = sigma_f1 - del_sigma_f;
    sigma_ff = sigma_f(d_i,d_fdes,sigma_max,alpha_i,theta_i,alpha_f,del_sigma_f,epsilon,sigma_i,sigma_f1,theta_f);
end
end


% Calculation of d_f
function d_f = d_f(d_i, theta_i, sigma_i, sigma_f, N, theta_f)
    if sigma_i == sigma_f
        d_f = (d_i / cos(theta_i)) * (exp((theta_f - theta_i) * cot(sigma_i))) * sin(sigma_f);
    else
        d_f = (d_i / cos(theta_i)) * ((sin(sigma_f)^(N / (N - 1))) / (sin(sigma_i)^(1 / (N - 1))));
    end
end
