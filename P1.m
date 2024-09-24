clc;clear;close all;


V = 10;
theta_i = deg2rad(70);
alpha_i = deg2rad(0);
alpha_f = deg2rad(0);
sigma_i  = deg2rad(30:15:80);
sigma_f = deg2rad(30:15:80);
d_i = 4;
R_i = d_i/cos(theta_i);

N1 = N(theta_i, alpha_f, sigma_i, sigma_f(end));
N2 = N(theta_i, alpha_f, sigma_i(end), sigma_f);

figure(1)
stem(rad2deg(sigma_i),N1,"filled",LineStyle=":",LineWidth=6.0,Color="k")
xlim([20,85])
ylim([0,1.2])
ax = gca;
ax.LineWidth = 4;
set(gca, 'FontSize', 30);
xlabel("\textbf{$\sigma_i$ , deg}", 'Interpreter', 'latex',"FontSize",30,"FontAngle","italic")
ylabel("\textbf{N}", 'Interpreter', 'latex',"FontSize",30,"FontAngle","italic")
grid on
figure(2)
stem(rad2deg(sigma_f),N2,"filled",LineStyle=":",LineWidth=6.0,Color="k")
ax = gca;
ax.LineWidth = 4;
xlim([20,85])
ylim([0,1.7])
grid on
set(gca, 'FontSize', 30);
xlabel("\textbf{$\sigma_f$ , deg}","FontSize",30, 'Interpreter', 'latex',"FontAngle","italic")
ylabel("\textbf{N}","FontSize",30, 'Interpreter', 'latex',"FontAngle","italic")


% Navigation Gain Calculation
function N = N(theta_i, alpha_f, sigma_i, sigma_f)
N = (theta_i - alpha_f + sigma_i)./(theta_i - alpha_f + sigma_f);
end
