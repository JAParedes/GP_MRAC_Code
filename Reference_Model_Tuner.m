close all
clear all
clc

%Code for tuning reference model response for Model Reference Adaptive
%Control (MRAC) reference tracking.

%A second order reference model is proposed.
omega_rm = 1;       %Natural frequency
eps_rm = 0.5;       %Damping ratio
A_rm = [0 1; -omega_rm^2 -2*eps_rm*omega_rm];   %Reference model state matrix

%Plant initial state
%x0 = [4;5];
x0 = [3.5;5];

dt = 0.05;                  %Sampling period of 0.05 seconds
T = 0:dt:60;                %Running time of a minute
X = zeros(2, length(T));    %State vector
X(:,1) = x0;

%Acceleration disturbances meant to create peaks in the reference states so
%that they don't stay at zero after settling from the initial state.
deltaXdd = zeros(size(T));

%Tripuls is used to create triangular pulses in the desired acceleration at
%given timestamps.
%This allows the reference states to change while keeping the reference
%dynamics consistent

%for x0 = [4;5];
% deltaXdd = deltaXdd + 3.75*tripuls(T-15,2.5);
% deltaXdd = deltaXdd - 3.75*tripuls(T-30,2.5);
% deltaXdd = deltaXdd + 3.75*tripuls(T-45,2.5);

%for x0 = [3.5; 5];
deltaXdd = deltaXdd + 10*tripuls(T-15,4);
deltaXdd = deltaXdd - 10*tripuls(T-30,4);
deltaXdd = deltaXdd + 10*tripuls(T-45,4);

%Simulation of reference model
for ii = 1:length(T)-1

[~, Yrm] = ode45(@(t, x) A_rm*x + [0; deltaXdd(ii)], [0 dt], X(:,ii));
X(:,ii+1) = Yrm(end,:)';

end

%Visualization of desired states.
%deltaXdd can be changed to modify the resulting peaks.

figure(1)

plot(T, X)
xlabel('Time (sec)', 'interpreter', 'latex')
ylabel('Reference States', 'interpreter', 'latex')
legend({'$\theta_{ref}$ (rad)','${\dot{\theta}}_{ref}$ (rad/s)'},'interpreter', 'latex')