close all
clear all
clc

%Code for tuning Model Reference Adaptive Controller (MRAC) Gains (K),
%without taking into account modelling errors (ideal case)
%Thus, Wing-Rock Dynamics are represented by a double integrator plant

dt = 0.05;                      %Sampling period of 0.05 seconds
T = 0:dt:60;                    %Running time of a minute


Y = zeros(2, length(T));        %State vector
Xrm = zeros(2, length(T));      %Reference model state vector
x0 = [3.5;5];                   %Initial state
Y(:,1) = x0;
Xrm(:,1) = x0;

%MRAC settings

%Gains for nu_pd feedback (MODIFY AS NEEDED FOR DESIRED IDEAL PERFORMANCE)

%K1 = 20;
%K2 = 20;

K1 = 10;
K2 = 10;

%Reference model parameters
%A second order reference model is proposed.
omega_rm = 1;   %Natural frequency
eps_rm = 0.5;   %Damping ratio

%Reference model state matrix
A_rm = [0 1; -omega_rm^2 -2*eps_rm*omega_rm];

%Acceleration disturbances meant to create peaks in the reference states so
%that they don't stay at zero after settling from the initial state.
%Use rm_test.m for testing different deltaXdd values to obtain differen
%reference states
deltaXdd = zeros(size(T));

%for x0 = [4;5];
% deltaXdd = deltaXdd + 3.75*tripuls(T-15,2.5);
% deltaXdd = deltaXdd - 3.75*tripuls(T-30,2.5);
% deltaXdd = deltaXdd + 3.75*tripuls(T-45,2.5);

%for x0 = [3.5; 5];
deltaXdd = deltaXdd + 10*tripuls(T-15,10);
deltaXdd = deltaXdd - 10*tripuls(T-30,10);
deltaXdd = deltaXdd + 10*tripuls(T-45,10);

%nu pseudocontrol input functions 
nu_rm =@(x_rm) [0 1]*A_rm*x_rm; %nu feedforward term (relies on reference model)
nu_pd = @(e) [K1 K2]*e;         %nu feedback term (relies on error feedback) (KEY TERM IN THIS TUNING PHASE)

%Wing-rock dynamics without modelling error (double integrator with control input guiding second integration term)
A_real = [0 1; 0 0];
B_real = [0;1];

for ii = 1:length(T)-1
    %ii
    %Reference model update
    [~, Yrm] = ode45(@(t, x) A_rm*x + [0; deltaXdd(ii)], [0 dt], Xrm(:,ii));
    Xrm(:,ii+1) = Yrm(end,:)';
    %Error between actual output and reference states
    error = Xrm(:,ii+1) - Y(:,ii);
    %Obtaining pseudo control input nu
    nu_rm_i = nu_rm(Xrm(:,ii+1));
    nu_pd_i = nu_pd(error);
    nu = nu_rm_i + nu_pd_i;
    %Simulating wing-rock dynamics, given u = nu
    [~, Yreal] = ode45(@(t, x) A_real*x + B_real*nu, [0 dt], Y(:,ii));
    Y(:,ii+1) = Yreal(end,:)';
end

%Visualization of ideal performance. Increasing K1 and K2 values leads to
%better (ideal) reference tracking but leads to greater input value
%requirements

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