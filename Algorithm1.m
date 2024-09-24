clc;
clear;
close all;

sigma_max = deg2rad(75);
d1_maxspec = 8;
d_i = 4;

alpha_f = deg2rad(0);
alpha_i = deg2rad(91);
del_alpha = 0.001;
epsilon = 0.01;
sigma_f = sigma_max;

theta_i1 = theta_i(d_i, d1_maxspec, sigma_max, alpha_f, del_alpha, epsilon, alpha_i,sigma_f);
theta_i1 = rad2deg(theta_i1);

% Algorithm 1
function theta = theta_i(d_i, d1_maxspec, sigma_max, alpha_f, del_alpha, epsilon, alpha_i,sigma_f)
    
    theta = alpha_i - sigma_max;

    sigma_i = alpha_i - theta;
    
    theta_f = sigma_f - alpha_f;
    N = (alpha_i - alpha_f) / (theta - alpha_f + sigma_f);
    
    d1_max_ = d1_max(d_i, theta, sigma_i, sigma_f, N, alpha_i);
    d_f_ = d_f(d_i, theta, sigma_i, sigma_f, N, theta_f);

    if (d1_maxspec - d1_max_) < epsilon
        
        return;
    else
        alpha_i = alpha_i + del_alpha;
        theta = theta_i(d_i, d1_maxspec, sigma_max, alpha_f, del_alpha, epsilon, alpha_i,sigma_f);
    end
end




% Calculation of d1_max
function d1_max = d1_max(d_i, theta_i, sigma_i, sigma_f, N, alpha_i)
    if sigma_i == sigma_f
        
        d1_max = (d_i / cos(theta_i)) * (exp((deg2rad(90) - sigma_i - theta_i) * cot(sigma_i))) * cos(deg2rad(90) - sigma_i);
    else
        d1_max = (d_i / cos(theta_i)) * ((sin((deg2rad(90) * (N - 1) + (alpha_i - N * theta_i)) / N) / sin(sigma_i))^(1 / (N - 1)))*(cos((deg2rad(90) - (alpha_i - N*(theta_i)))/N));
        
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
