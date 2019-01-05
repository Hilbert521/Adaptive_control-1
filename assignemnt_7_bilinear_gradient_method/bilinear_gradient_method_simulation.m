clear all;
close all;
clc;

%%%% assignment 7, task 4.13. the equations used are on page 224 in I & S

%% Model: m*y'' + beta*y' + k*y = u
m = 20;
beta = 0.1;
k = 5;


%% state space model x' = Ax + Bu
%x = [y; y dot]
A = [  0      1;
     -k/m -beta/m];
B = [0; 1/m];

%% filter Lambda = s^2 + lambda1*s + lambda0
lambda1 = 2.5;
lambda0 = 2;

[A_f0, B_f0, C_f0, D_f0] = tf2ss([0, 0, 1], [1, lambda1, lambda0]);
[A_f1, B_f1, C_f1, D_f1] = tf2ss([0, -1, 0],[1, lambda1, lambda0]);
[A_f2, B_f2, C_f2, D_f2] = tf2ss([-1, 0, 0],[1, lambda1, lambda0]);

%% parameters
Gamma = [15 0;
         0  0.5];
gamma = 4;

%% simulation
time = 600; %seconds
h = 0.01;   %time step
N = time/h; %number of steps

%allocate memory
t = zeros(1, N);                %1 row, N columns
x = zeros(2, N);
u = zeros(1, N);
theta = zeros(2, N);
epsilon = zeros(1, N);
phi = zeros(2, N);
x_z = zeros(2, N);
x_z1 = zeros(2, N);
x_phi1 = zeros(2, N);
x_phi2 = zeros(2, N);
phi = zeros(2, N);

theta = zeros(2, N);
z_hat = zeros(1,N);
rho = zeros(1,N);


Matrix = ones(2, N);
theta_star = zeros(2,N);
theta_star(1, :) = Matrix(1, :)*m;
theta_star(2, :) = Matrix(2, :)*beta;


% initial values
x(:,1) = [0; 0];
theta(:, 1) = [1; 1];
rho(1) = 0.1;


for n=1:(N-1)
    t(n+1) = n*h;               
    
    u = 10*sin(0.5*t(n)) + 7*cos(t(n));
    %u = 10;
    u(n) = u;
    
    %simulate the system with Euler integration       
    x_dot = A*x(:, n) + B*u(n);    
    x(:, n+1) = x(:, n) + h*x_dot;
    y(:, n) = x(1, n);
    
    %find z1 with euler integration
    x_z1_dot = A_f0*x_z1(:, n) + B_f0*u(:,n);
    x_z1(:, n+1) = x_z1(:, n) + h*x_z1_dot;
    z1(:, n) = C_f0*x_z1(:, n) + D_f0*u(:,n);    
    
    %find z with euler integration
    x_z_dot = A_f0*x_z(:, n) + B_f0*y(:,n);
    x_z(:, n+1) = x_z(:, n) + h*x_z_dot;
    z(:, n) = C_f0*x_z(:, n) + D_f0*y(:,n);
    
    %find phi with euler integration
    x_phi1_dot = A_f2*x_phi1(:, n) + B_f2*y(n);
    x_phi1(:, n+1) = x_phi1(:, n) + h*x_phi1_dot;
    phi1(n) = C_f2*x_phi1(:, n) + D_f2*y(n);
    
    x_phi2_dot = A_f1*x_phi2(:,n) + B_f1*y(n);
    x_phi2(:,n+1) = x_phi2(:, n) + h*x_phi2_dot;
    phi2(n) = C_f1*x_phi2(:, n) + D_f1*y(n);
    
    phi(:, n) = [phi1(n); phi2(n)];
    
    %find z_hat
    z_hat(n) = rho(n)*(theta(:,n)'*phi(:,n) + z1(n));
    
    
    %find epsilon
    n_squared(n) = phi(:, n)'*phi(:,n) + z1(:, n)^2;
    m_squared(n) = 1 + n_squared(n);
    epsilon(n) = (z(n)-z_hat(n))/m_squared(n);
   % epsilon(n) = (z(:,n)-z_hat(:,n));
    
    %find rho
    rho_dot = gamma*epsilon(n)*(theta(:,n)'*phi(:,n) + z1(n));
    rho(n+1) = rho(n) + h*rho_dot;
    
   
    
    %theta
    theta_dot = Gamma*epsilon(n)*phi(:,n);
    theta(:,n+1) = theta(:, n) + h*theta_dot; 
    
end

%% plot


figure(3);
grid on;
title('Plot of theta and \rho');
hold on;
plot(t, theta);
plot(t, rho);
legend('m', '\beta', '\rho');






