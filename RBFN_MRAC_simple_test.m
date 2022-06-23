close all
clear all
clc

%Code for testing a Model Reference Adaptive Control (MRAC) scheme to
%control the pitch of a plane under wing-rock dynamics, accounting for
%modeling errors by using a Gaussian Radial Basis Function Network (RBFN).

%Wing-rock dynamics modelling error variance
omega_n = 0.01;

%State measurement standard deviation
std_meas = 0.1;

dt = 0.05;                                  %Sampling period of 0.05 seconds
T = 0:dt:60;                                %Running time of a minute

Y = zeros(2, length(T));                    %State vector
Y_meas = zeros(2, length(T));               %Measured state vector
Xrm = zeros(2, length(T));                  %Reference model state vector
%x0 = [4; 5];
x0 = [3.5; 5];                              %Initial state
Y(:,1) = x0;
Y_meas(:,1) = x0;
Xrm(:,1) = x0;
delta_true_vec = zeros(1, length(T));       %Vector of true modeling error values
delta_RBFN_vec = zeros(1, length(T));       %Vector of estimated modeling error values

%MRAC settings

%Gains for nu_pd feedback (use MRAC_Ideal_Tuner.m for tuning)
%K1 = 20;
%K2 = 20;

%K1 = 2;
%K2 = 2;

K1 = 10;
K2 = 10;

%Paramaters and function for MRAC RBFN weight update

A = [0 1; -K1 -K2];
B = [0; 1];

Q = eye(2);
P = lyap(A',Q);
Gamma_W = 20;       %Learning rate (TUNE AS NEEDED)
max_w = 10;         %Maximum individual weight magnitude

%Projection operator saturates the values of the updated weights depedning
%on the maximum magnitude stated by max_w
proj = @(e,Theta,Wp, dT) max(min(Wp - (dT*Gamma_W*e'*P*B)*Theta, max_w),-max_w); 

%Number of centers of the RBFN
nRBFN = 121;

%Setting RBFN centers evenly over a [-2 2] x [-2 2] state grid
%Spreading the centers over a greater grid may increase uncertainty
%estimation performance, although such tests haven't been made for the
%project

cij = zeros(2,nRBFN);
kk = 1;
for ii = -5:5
    for jj = -5:5
        cij(:,kk) = [0.4*ii; 0.4*jj];
        kk = kk+1;
    end
end
W_RBFN = ones(nRBFN,1);         %RBFN weights
mu = 1;                         %Gaussian RBF spread
RBFN = @(x) exp(-sum((x-cij).^2)/(2*mu^2)); %Gaussian RBFN function

%Reference model parameters
%A second order reference model is proposed.
omega_rm = 1;   %Natural frequency
eps_rm = 0.5;   %Damping ratio

%Reference model state matrix
A_rm = [0 1; -omega_rm^2 -2*eps_rm*omega_rm];

%Acceleration disturbances meant to create peaks in the reference states so
%that they don't stay at zero after settling from the initial state.
%Use Reference_Model_Tuner.m for testing different deltaXdd values to obtain differen
%reference states
deltaXdd = zeros(size(T));

%for x0 = [4;5];
%deltaXdd = deltaXdd + 3.75*tripuls(T-15,2.5);
%deltaXdd = deltaXdd - 3.75*tripuls(T-30,2.5);
%deltaXdd = deltaXdd + 3.75*tripuls(T-45,2.5);

%for x0 = [3.5; 5];
%deltaXdd = deltaXdd +6*tripuls(T-15,4);
%deltaXdd = deltaXdd - 6*tripuls(T-30,4);
%deltaXdd = deltaXdd + 6*tripuls(T-45,4);

deltaXdd = deltaXdd + 10*tripuls(T-15,10);
deltaXdd = deltaXdd - 10*tripuls(T-30,10);
deltaXdd = deltaXdd + 10*tripuls(T-45,10);

%nu pseudocontrol input functions 
nu_rm =@(x_rm) [0 1]*A_rm*x_rm;     %nu feedforward term (relies on reference model)
nu_pd = @(e) [K1 K2]*e;             %nu feedback term (relies on error feedback)
nu_ad_RBFN =@(x,Wp) Wp'*(RBFN(x))'; %nu adaptive term (relies on RBFN modeling error estimation)

for ii = 1:length(T)-1
    %ii
    %Reference model update
    [~, Yrm] = ode45(@(t, x) A_rm*x + [0; deltaXdd(ii)], [0 dt], Xrm(:,ii));
    Xrm(:,ii+1) = Yrm(end,:)';
    %Error between actual output and reference states
    error = Xrm(:,ii+1) - Y_meas(:,ii);
    %Updating RBFN weights
    theta = RBFN(Y_meas(:,ii));
    W_RBFN = proj(error, theta', W_RBFN, dt);
    
    %Obtaining pseudo control input nu
    nu_rm_i = nu_rm(Xrm(:,ii+1));
    nu_pd_i = nu_pd(error);
    nu_ad_RBFN_i = nu_ad_RBFN(Y_meas(:,ii), W_RBFN);
    
    %Obtaining estimated modeling error via Gaussian RBFN
    delta_RBFN_vec(1, ii+1) = nu_ad_RBFN_i;
    
    nu = nu_rm_i + nu_pd_i - nu_ad_RBFN_i;
    
    %Simulating wing-rock dynamics, given u = nu
    %by obtaining Stochastic Differential Equation result
    [res, delta_true] = wrd_sim(Y(:,ii), nu, omega_n, dt);
    Y(:,ii+1) = res;                                %Actual state measurement
    Y_meas(:,ii+1) = res + std_meas*randn(2, 1);    %Noisy state measurement
    delta_true_vec(1, ii+1) = delta_true;           %Actual modeling error
end

MSE1 = (1/length(T))*((Y(1,:)-Xrm(1,:))*(Y(1,:)-Xrm(1,:))');
MSE2 = (1/length(T))*((Y(2,:)-Xrm(2,:))*(Y(2,:)-Xrm(2,:))');
MSE3 = (1/length(T))*((delta_RBFN_vec-delta_true_vec)*(delta_RBFN_vec-delta_true_vec)');

fprintf("MSE for X1: %d \n", MSE1)          %MSE for angle
fprintf("MSE for X2: %d \n", MSE2)          %MSE for angular velocity
fprintf("MSE for delta X: %d \n", MSE3)     %MSE for uncertainty estimation

%Ploting state tracking and uncertainty estimation results.
%Uncertainty estimation becomes better if the spread of the RBFN covers the
%desired reference states and the learning rate Gamma_W is increased.

figure(1)
hold on
plot(T, Y(1,:), 'b', 'Linewidth', 2)
plot(T, Xrm(1, :), 'r--', 'Linewidth', 2)
xlabel('Time (sec)', 'interpreter', 'latex')
ylabel('$\theta$ (rad)', 'interpreter', 'latex')
legend({'MRAC Output','Reference Model'},'interpreter', 'latex')
hold off

figure(2)
hold on
plot(T, Y(2,:), 'b', 'Linewidth', 2)
plot(T, Xrm(2, :), 'r--', 'Linewidth', 2)
xlabel('Time (sec)', 'interpreter', 'latex')
ylabel('$\dot{\theta}$ (rad/s)', 'interpreter', 'latex')
legend({'MRAC Output','Reference Model'},'interpreter', 'latex')
hold off

figure(3)
hold on
plot(T, delta_RBFN_vec, 'b', 'Linewidth', 2)
plot(T, delta_true_vec, 'r--', 'Linewidth', 2)
xlabel('Time (sec)', 'interpreter', 'latex')
ylabel('$\Delta (x)$', 'interpreter', 'latex')
legend({'Estimated $\Delta (x)$','True $\Delta (x)$'},'interpreter', 'latex')
hold off

%Wing-rock dynamics one-step simulation, using the Euler Maruyama method
%for Stochastic Differential Equations (SDEs)
function [res, dx] = wrd_sim(x, nu, omega_n, dt)

    %Wing-rock dynamics modelling error mean
    mean_x = 0.8 + 0.2314*x(1) + 0.6918*x(2) - 0.6245*abs(x(1))*x(2) + ...
        0.0095*abs(x(2))*x(2) + 0.0214*x(1)^3;
    mu = [x(2); nu + mean_x];
    sigma = [0; sqrt(omega_n)];
    dW = randn(1,1)*sqrt(dt);
    res = x + mu*dt + sigma*dW;
    dx = mean_x + sqrt(omega_n)*dW/sqrt(dt);

end