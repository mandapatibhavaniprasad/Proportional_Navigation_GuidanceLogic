clc;
clear;
close all;
global sigma_i theta_i sigma_f alpha_f d_i

V = 2;
theta_i = deg2rad(67.8);
sigma_i = deg2rad(75);
alpha_i = sigma_i + theta_i;


alpha_f = deg2rad(0);

sigma_f_arr = deg2rad([40.8746,53.3365,63.9075,73.0061]);

d_i = 4;
% R_i = d_i/cos(theta_i);
R_i = 10.6;
figure(1)
plot(0,4,"k*",MarkerSize=20)
for i = 1:length(sigma_f_arr)
    sigma_f = sigma_f_arr(i);


    tRange = linspace(0,30,3000);
    x0 = -R_i*cos(theta_i);
    y0 = - R_i*sin(theta_i)+4;

    N_ = N(theta_i, alpha_f, sigma_i, sigma_f);

    options = odeset("Events",@Eventfunc,"RelTol",1e-10,"AbsTol",1e-10);

    [tSol,YSol_PN] = ode45(@(t, y) PN(t,y,sigma_i,theta_i,V,N_), tRange,[R_i, theta_i,x0,y0],options);
    R = YSol_PN(:,1);
    theta = YSol_PN(:,2);
    X = YSol_PN(:,3);
    Y = YSol_PN(:,4);
    N_ = N(theta_i, alpha_f, sigma_i, sigma_f);
    sigma = (sigma_i+theta_i - N_*theta_i) - (1-N_)*theta;
    alpha = sigma+theta;
    theta_dot = -V.*sin(sigma)./R;
    alpha_dot = N_.*theta_dot;

    figure(1)
    hold on
    plot(X,Y,LineWidth=6.0)
    ax = gca;
    ax.LineWidth = 4;
    set(gca, 'FontSize', 30);
    xlim([-16,10])
    ylim([-8,10])
    legend({'beacon', '$d_f = 2$ m', '$d_f = 3$ m', '$d_f = 4$ m', '$d_f = 5$ m'}, 'Interpreter', 'latex', 'Location', 'southeast');
    h_legend = legend;
    set(h_legend, 'FontSize',30);
    xlabel("\textbf{X, m}", 'Interpreter', 'latex',"FontSize",30)
    ylabel("\textbf{Y, m}", 'Interpreter', 'latex',"FontSize",30)
    
    grid on
    hold on

end

hold off

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
global sigma_i theta_i sigma_f alpha_f

N_ =  N(theta_i, alpha_f, sigma_i, sigma_f);

theta = y(2);
sigma = (sigma_i+theta_i - N_*theta_i) - (1-N_)*theta;

val = sigma+theta - alpha_f;
isterminal = 1;
direction = -1;
end

