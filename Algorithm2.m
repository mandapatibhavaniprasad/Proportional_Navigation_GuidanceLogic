clc;
clear;
close all;

sigma_max = deg2rad(75);

d_i = 4;
alpha_f = deg2rad(0);
theta_i = deg2rad(67.8);
% d_fdes = d_i;
d_fdes = 5;
sigma_f = sigma_max;
del_alpha = 0.0001;
epsilon = 0.005;

alpha_i_max = theta_i+sigma_max;
alpha_i1 = alpha_i_max;
alpha_i2 = alpha_i(d_i,d_fdes,sigma_max,theta_i,alpha_f,del_alpha,epsilon,sigma_f,alpha_i_max,alpha_i1);
alpha_i2 = rad2deg(alpha_i2);

% Algorithm 2
function alpha_ii = alpha_i(d_i,d_fdes,sigma_max,theta_i,alpha_f,del_alpha,epsilon,sigma_f,alpha_i_max,alpha_i1)
sigma_i = alpha_i1 - theta_i;
theta_f = sigma_f - alpha_f;
N = (alpha_i1 - alpha_f) / (theta_i - alpha_f + sigma_f);
d_f_ = d_f(d_i, theta_i, sigma_i, sigma_f, N, theta_f);
if abs(d_f_ - d_fdes) < epsilon
    alpha_ii = alpha_i1;
else
    alpha_i1 = alpha_i1 - del_alpha;

    alpha_ii = alpha_i(d_i,d_fdes,sigma_max,theta_i,alpha_f,del_alpha,epsilon,sigma_f,alpha_i_max,alpha_i1);


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
