clc;clear;close all;
global sigma_i theta_i sigma_f alpha_f 

V = 2;
theta_i = deg2rad(70);

alpha_f = deg2rad(0);

sigma_i_arr = deg2rad(30:15:80);
sigma_f_arr = deg2rad(30:15:80);

d_i = 4;
% R_i = d_i/cos(theta_i);
R_i = 10.6;
R_prev = 0;
figure(1)
plot(0,4,"k*",MarkerSize=9)
figure(2)
plot(0,4,"k*",MarkerSize=9)
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
    
    X = YSol_PN(:,3);
    Y = YSol_PN(:,4);

    figure(1)
    
    hold on
    plot(X,Y,LineWidth=6.0)
    set(gca, 'FontSize', 30);
    legend('beacon', '$\sigma_i = 30$ deg', '$\sigma_i = 45$ deg', '$\sigma_i = 60$ deg', '$\sigma_i = 75$ deg', 'Interpreter', 'latex',"Location","southeast")
    h_legend = legend;
    set(h_legend, 'FontSize',30);
    xlabel("\textbf{X, m}", 'Interpreter', 'latex',"FontSize",30)
    ylabel("\textbf{Y, m}", 'Interpreter', 'latex',"FontSize",30)
    xlim([-17,7.5])
    ylim([-8,11])
    ax = gca;
    ax.LineWidth = 4;
    grid on
  
    
end

for j = 1:length(sigma_i_arr)
    sigma_i = sigma_i_arr(end);
    sigma_f = sigma_f_arr(j);
    tRange = linspace(0,50,3000);
    N_ = N(theta_i, alpha_f, sigma_i, sigma_f);
    x0 = -R_i*cos(theta_i);
    y0 = - R_i*sin(theta_i)+4;
    options = odeset("Events",@Eventfunc,"RelTol",1e-10,"AbsTol",1e-10);

    [tSol,YSol_PN] = ode45(@(t, y) PN(t,y,sigma_i,theta_i,V,N_), tRange,[R_i, theta_i,x0,y0],options);
    R = YSol_PN(:,1);
    theta = YSol_PN(:,2);
    X = YSol_PN(:,3);
    Y = YSol_PN(:,4);

    figure(2)

    hold on
    plot(X,Y,LineWidth=6.0)
    legend('beacon', '$\sigma_f = 30$ deg', '$\sigma_f = 45$ deg', '$\sigma_f = 60$ deg', '$\sigma_f = 75$ deg', 'Interpreter', 'latex',"Location","southeast")
    h_legend = legend;
    set(h_legend, 'FontSize',30);
    set(gca, 'FontSize', 30);

    xlabel("\textbf{X, m}", 'Interpreter', 'latex',"FontSize",30)
    ylabel("\textbf{Y, m}", 'Interpreter', 'latex',"FontSize",30)
    xlim([-17,7.5])
    ylim([-8,11])
    ax = gca;
    ax.LineWidth = 4;
    grid on
    
    
end


% Navigation Gain Calculation
function N = N(theta_i, alpha_f, sigma_i, sigma_f)
N = (theta_i - alpha_f + sigma_i)./(theta_i - alpha_f + sigma_f);
end

% LOS Range calculation
function dYdt = PN(~, y,sigma_i,theta_i,V,N)
R = y(1);
theta = y(2);
x_UAV = y(3);
y_UAV = y(4);
sigma = (sigma_i+theta_i - N*theta_i) - (1-N)*theta;
x_UAV_dot = V*cos(sigma+theta);
y_UAV_dot = V*sin(sigma+theta);
R_dot = -V*cos(sigma);
theta_dot = -V*sin(sigma)/R;
dYdt = [R_dot;theta_dot;x_UAV_dot;y_UAV_dot];
end


function [val, isterminal, direction] = Eventfunc(~,y)
global sigma_i theta_i sigma_f alpha_f 
N_ =  N(theta_i, alpha_f, sigma_i, sigma_f);
theta = y(2);
sigma = (sigma_i+theta_i - N_*theta_i) - (1-N_)*theta;

val = sigma+theta;
isterminal = 1;
direction = -1;
end

