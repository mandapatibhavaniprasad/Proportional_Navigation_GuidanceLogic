clc;clear;close all;
global sigma_f sigma_i theta_i alpha_f

V = 2;
theta_i = deg2rad(70);

alpha_f = deg2rad(0);

sigma_i_arr = deg2rad(30:15:80);
sigma_f_arr = deg2rad(30:15:80);

d_i = 4;
% R_i = d_i/cos(theta_i);
R_i = 10.6;

for i = 1:length(sigma_i_arr)
    sigma_i = sigma_i_arr(i);
    sigma_f = sigma_f_arr(end);
    tRange = linspace(0,50,3000);
    x0 = -R_i*cos(theta_i);
    y0 = - R_i*sin(theta_i)+4;

    N_ = N(theta_i, alpha_f, sigma_i, sigma_f);

    options = odeset("Events",@Eventfunc,"RelTol",1e-10,"AbsTol",1e-10);

    [tSol,YSol_PN] = ode45(@(t, y) PN(t,y,sigma_i,theta_i,V,N_), tRange,[R_i, theta_i,x0,y0],options);
    R = YSol_PN(:,1);
    theta = YSol_PN(:,2);
    N_ = N(theta_i, alpha_f, sigma_i, sigma_f);
    % disp(N_)
    sigma = (sigma_i+theta_i - N_*theta_i) - (1-N_)*theta;
    figure(1)
    plot(tSol,R,LineWidth=6.0)
    hold on
    legend('$\sigma_i = 30$ deg', '$\sigma_i = 45$ deg', '$\sigma_i = 60$ deg', '$\sigma_i = 75$ deg', 'Interpreter', 'latex')
    h_legend = legend;
    set(h_legend, 'FontSize',30);
    grid on
    xlim([-2,12])
    ylim([-2,12])
    ax = gca;
    ax.LineWidth = 4;
    set(gca, 'FontSize', 30);
    xlabel("\textbf{t, s}", 'Interpreter', 'latex', "FontSize",30,"FontAngle","italic")
    ylabel("\textbf{R, m}", 'Interpreter', 'latex', "FontSize",30,"FontAngle","italic")
    d1 = R.*cos(theta);
    % d1_ = d1(theta, d_i, theta_i, sigma_i, sigma, N_, sigma_f);
    
    figure(2)
    plot(rad2deg(theta),d1,LineWidth=6.0)
    set(gca, 'FontSize', 30);
    legend('$\sigma_i = 30$ deg', '$\sigma_i = 45$ deg', '$\sigma_i = 60$ deg', '$\sigma_i = 75$ deg', 'Interpreter', 'latex',"Location","southeast")
    h_legend = legend;
    set(h_legend, 'FontSize',30);
    xlim([-90 90])
    ylim([0,10])
    ax = gca;
    ax.LineWidth = 4;
    grid on
    xlabel("\textbf{$\theta$, deg}", 'Interpreter', 'latex',"FontSize",30,"FontAngle","italic")
    ylabel("\textbf{$d_1$, m}", 'Interpreter', 'latex',"FontSize",30,"FontAngle","italic")
   
    
    hold on

end

for j = 1:length(sigma_i_arr)
    sigma_i = sigma_i_arr(end);
    sigma_f = sigma_f_arr(j);
    tRange = linspace(0,50,3000);
    N_ = N(theta_i, alpha_f, sigma_i, sigma_f);
    x0 = -d_i*cos(sigma_i+theta_i);
    y0 = - d_i*sin(sigma_i+theta_i)+4;
    options = odeset("Events",@Eventfunc,"RelTol",1e-10,"AbsTol",1e-10);

    [tSol,YSol_PN] = ode45(@(t, y) PN(t,y,sigma_i,theta_i,V,N_), tRange,[R_i, theta_i,x0,y0],options);
    R = YSol_PN(:,1);
    theta = YSol_PN(:,2);
    
    % disp(N_)
    sigma = (sigma_i+theta_i - N_*theta_i) - (1-N_)*theta;
    figure(3)
    plot(tSol,R,LineWidth=6.0)
    hold on
    legend('$\sigma_f = 30$ deg', '$\sigma_f = 45$ deg', '$\sigma_f = 60$ deg', '$\sigma_f = 75$ deg', 'Interpreter', 'latex')
    h_legend = legend;
    set(gca, 'FontSize', 30);
    grid on
    set(h_legend, 'FontSize',30);
    ax = gca;
    ax.LineWidth = 4;
    xlabel("\textbf{t, s}", 'Interpreter', 'latex', "FontSize",30,"FontAngle","italic")
    ylabel("\textbf{R, m}", 'Interpreter', 'latex', "FontSize",30,"FontAngle","italic")
    xlim([-2,12])
    ylim([-2,12])
    
    
    d1 = R.*cos(theta);
    % d1_ = d1(theta, d_i, theta_i, sigma_i, sigma, N_, sigma_f);

    figure(4)
    plot(rad2deg(theta),d1,LineWidth=6.0)
    legend('$\sigma_f = 30$ deg', '$\sigma_f = 45$ deg', '$\sigma_f = 60$ deg', '$\sigma_f = 75$ deg', 'Interpreter', 'latex',"Location","southeast")
    h_legend = legend;
    set(h_legend, 'FontSize',30);
    grid on
    set(gca, 'FontSize', 30);
    xlim([-90 90])
    ylim([0,10])
    ax = gca;
    ax.LineWidth = 4;
    xlabel("\textbf{$\theta$, deg}", 'Interpreter', 'latex',"FontSize",30,"FontAngle","italic")
    ylabel("\textbf{$d_1$, m}", 'Interpreter', 'latex',"FontSize",30,"FontAngle","italic")
    hold on
end



% Navigation Gain Calculation
function N = N(theta_i, alpha_f, sigma_i, sigma_f)
N = (theta_i - alpha_f + sigma_i)./(theta_i - alpha_f + sigma_f);
end

% LOS Range calculation
function dYdt = PN(~, y,sigma_i,theta_i,V,N)
R = y(1);
theta = y(2);
% x_UAV = y(3);
% y_UAV = y(4);
sigma = (sigma_i+theta_i - N*theta_i) - (1-N)*theta;
x_UAV_dot = V*cos(sigma+theta);
y_UAV_dot = V*sin(sigma+theta);
R_dot = -V*cos(sigma);
theta_dot = -V*sin(sigma)/R;
dYdt = [R_dot;theta_dot;x_UAV_dot;y_UAV_dot];
end


function [val, isterminal, direction] = Eventfunc(~,y)

global sigma_f sigma_i theta_i alpha_f
N_ = N(theta_i, alpha_f, sigma_i, sigma_f);
theta = y(2);
sigma = (sigma_i+theta_i - N_*theta_i) - (1-N_)*theta;
val = sigma+theta - alpha_f;
isterminal = 1;
direction = -1;

end
