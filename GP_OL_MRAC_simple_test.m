close all
clear all
clc

%Code for testing a Model Reference Adaptive Control (MRAC) scheme to
%control the pitch of a plane under wing-rock dynamics, accounting for
%modeling errors by using Gaussian Process regression (hence, the algorithm
%is called GP-MRAC). In this case, the oldest data point in the training 
% set will be discarded.

%Wing-rock dynamics modelling error variance
omega_n = 0.01;

%State measurement standard deviation
std_meas = 0.1;

dt = 0.05;                              %Sampling period of 0.05 seconds
T = 0:dt:60;                            %Running time of a minute

Y = zeros(2, length(T));                %State vector
Y_meas = zeros(2, length(T));           %Measured state vector
Xrm = zeros(2, length(T));              %Reference model state vector
%x0 = [4; 5];
x0 = [3.5; 5];                          %Initial state
Y(:,1) = x0;
Y_meas(:,1) = x0;
Xrm(:,1) = x0;

delta_true_vec = zeros(1, length(T));   %Vector of true modeling error values
delta_GP_vec = zeros(1, length(T));     %Vector of estimated modeling error values

delta_GP_cov = zeros(1, length(T));     %Vector of estimated modeling error covariance (estimated by Gaussian Process Regression)

%Since state acceleration will be estimated by a Kalman Filter, the
%estimated and the actual accelerations will be stored in case a comparison
%is required and to test the effect of a perfect acceleration feedback
accel_true_vec = zeros(1, length(T));   %True acceleration vector
accel_est_vec = zeros(1, length(T));    %Estimated acceleration vector via Kalman Filter

%MRAC settings

%Gains for nu_pd feedback (use MRAC_Ideal_Tuner.m for tuning)
%K1 = 20;
%K2 = 20;

%K1 = 2;
%K2 = 2;

K1 = 10;
K2 = 10;

%Kalman filter parameters for acceleration (xdd) estimation

A_k = expm([0 1 0; 0 0 1; 0 0 0]*dt);
H_k = [1 0 0; 0 1 0];
%Q_k = [0.01 0 0; 0 0.01 0; 0 0 10];
Q_k = [0.01 0 0; 0 0.01 0; 0 0 0.5];
R_k = 0.01*eye(2);
x_KF_est = [x0; 0];
P_KF_est = 1*eye(3);

%Gaussian Process Regression settings

%Kernel settings for square exponential kernel
tau = 1;                %Kernel magnitude
%log_width = -0.1;      %Kernel spread
log_width = -0.001;
%Kernel function
kernel_sq = @(x,xp) (tau^2).*exp(-0.5.*((x-xp).^2)./((10^log_width)^2));

%Assumed logarithmic measurement noise (affects spread of estimation covariance)
%log_noise_learn = log10(1e-4);
log_noise_learn = log10(1e0);

pmax = 100;                 %Maximum number of GP training data points
XDATA = zeros(pmax, 2);     %Training state data points
XDATA(1,:) = x0';
YDATA = zeros(pmax, 1);     %Training measurement (modeling uncertainty error) data points
pnum = 0;                   %Current number of GP training data points
last_pnum = 1;              %Index of last data point added to training set

%Tolerance for gamma value. If gamma > eps_tol, then a new data point must
%be added to the Gaussian Process Regression training set.
%eps_tol = 1e-4;
%eps_tol = 1e-7;
eps_tol = 1e-9;

xdd_t = 0;                  %True angle acceleration

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
%deltaXdd = deltaXdd + 6*tripuls(T-15,4);
%deltaXdd = deltaXdd - 6*tripuls(T-30,4);
%deltaXdd = deltaXdd + 6*tripuls(T-45,4);

deltaXdd = deltaXdd + 10*tripuls(T-15,10);
deltaXdd = deltaXdd - 10*tripuls(T-30,10);
deltaXdd = deltaXdd + 10*tripuls(T-45,10);

%nu pseudocontrol input functions 
nu_rm =@(x_rm) [0 1]*A_rm*x_rm;     %nu feedforward term (relies on reference model)
nu_pd = @(e) [K1 K2]*e;             %nu feedback term (relies on error feedback)
nu_pre = 0;                         %previous nu value, used for modeling error estimation

for ii = 1:length(T)-1
    %ii
    %Reference model update
    [~, Yrm] = ode45(@(t, x) A_rm*x + [0; deltaXdd(ii)], [0 dt], Xrm(:,ii));
    Xrm(:,ii+1) = Yrm(end,:)';
    %Error between actual output and reference states
    error = Xrm(:,ii+1) - Y_meas(:,ii);
    
    %Obtaining pseudo control input nu feedworward and feedback terms
    nu_rm_i = nu_rm(Xrm(:,ii+1));
    nu_pd_i = nu_pd(error);
    
    %Acceleration KF and modeling error (Delta x) approximation
    [x_KF_est, P_KF_est] = KF_xdd(Y_meas(:,ii), x_KF_est, P_KF_est, A_k, H_k, R_k, Q_k);
    xdd_est = x_KF_est(3);              %Acceleration estimation
    accel_est_vec(1,ii) = xdd_est;
    
    delta = accel_true_vec(1,ii) - nu_pre;          %Actual modeling error
    delta_prop = accel_est_vec(1,ii) - nu_pre;      %Modeling estimation error using acceleration estimate from Kalman Filter
    
    %Updating Gaussian Process data set
    if ii == 1
        %Initial loop will always add data
        gamma = 1;
    else
        %Otherwise, the value of gamma must be computed to determine
        %whether the new data point state and its modelling error
        %measurement must be added
        gamma = KLI_test(Y_meas(:,ii)', XDATA(1:pnum,:), KZtau, kernel_sq);
    end
    
    if gamma > eps_tol
        %If the training set is full, delete the last added datapoint and
        %its corresponding measurement to add the new elements
        if pnum == pmax
            XDATA(last_pnum, :) = Y_meas(:,ii)';
            %YDATA(last_pnum,:) = delta;             %Use this for perfect (actual) modeling error feedback
            YDATA(last_pnum,:) = delta_prop;        %Use this for estimated modeling error feedback
            last_pnum = last_pnum + 1;
            if last_pnum > pmax
               last_pnum = 1;
            end
        else
            pnum = pnum + 1;
            XDATA(pnum,:) = Y_meas(:,ii)';
            %YDATA(pnum,:) = delta;                  %Use this for perfect (actual) modeling error feedback
            YDATA(pnum,:) = delta_prop;            %Use this for estimated modeling error feedback
        end
        %After adding new data, compute new KZtau and Ctau, and obtain
        %values for the mean and covariance of the modelling error estimate
        %at the current state.
        [delta_mean, delta_cov, KZtau, Ctau] = gpr_full(Y_meas(:,ii)', XDATA(1:pnum,:), YDATA(1:pnum,:), 10^log_noise_learn, kernel_sq);
    else
        %If no new data is added, compute the mean and the covariance of
        %the modeling error estimate using previously calculated Ctau
        [delta_mean, delta_cov] = gpr_update(Y_meas(:,ii)', XDATA(1:pnum,:), YDATA(1:pnum,:), Ctau, kernel_sq);
    end
    
    delta_GP_vec(1, ii+1) = delta_mean;
    delta_GP_cov(1, ii+1) = delta_cov;
    
    %Obtaining pseudo control input nu
    nu_ad_GP_i = delta_mean;                %nu adaptive term to counter modeling error
    nu = nu_rm_i + nu_pd_i - nu_ad_GP_i;    %Adding all nu terms
    
    %Simulating wing-rock dynamics, given u = nu
    %by obtaining Stochastic Differential Equation result
    [res, delta_true] = wrd_sim(Y(:,ii), nu, omega_n, dt);
    Y(:,ii+1) = res;                                %Actual state measurement
    Y_meas(:,ii+1) = res + std_meas*randn(2, 1);    %Noisy state measurement
    delta_true_vec(1, ii+1) = delta_true;           %Actual modeling error
    
    xdd_t = nu + delta_true;                        %Actual state acceleration
    nu_pre = nu;                                    %Saving previous nu value for acceleration estimation
    
    accel_true_vec(1,ii+1) = xdd_t;
end

MSE1 = (1/length(T))*((Y(1,:)-Xrm(1,:))*(Y(1,:)-Xrm(1,:))');
MSE2 = (1/length(T))*((Y(2,:)-Xrm(2,:))*(Y(2,:)-Xrm(2,:))');
MSE3 = (1/length(T))*((delta_GP_vec-delta_true_vec)*(delta_GP_vec-delta_true_vec)');

fprintf("MSE for X1: %d \n", MSE1)          %MSE for angle
fprintf("MSE for X2: %d \n", MSE2)          %MSE for angular velocity
fprintf("MSE for delta X: %d \n", MSE3)     %MSE for uncertainty estimation

%These instructions are used to plot standard deviation areas around the
%estimated Gaussian Process mean.
fillT = [T fliplr(T)];
sigma1 = sqrt(delta_GP_cov);
sigma1t = delta_GP_vec + 2.*sigma1;
sigma1b = delta_GP_vec - 2.*sigma1;
sigma1Fill = [sigma1t fliplr(sigma1b)]; 

%Ploting state tracking and uncertainty estimation results.
%The log_noise_learn has the greatest impact on uncertainty estimation.

figure(1)
hold on
plot(T, Y(1,:), 'b', 'Linewidth', 2)
plot(T, Xrm(1, :), 'r--', 'Linewidth', 2)
xlabel('Time (sec)', 'interpreter', 'latex')
ylabel('$\theta$ (rad)', 'interpreter', 'latex')
legend({'GP-MRAC Output (OL)','Reference Model'},'interpreter', 'latex', 'location', 'best')
hold off

figure(2)
hold on
plot(T, Y(2,:), 'b', 'Linewidth', 2)
plot(T, Xrm(2, :), 'r--', 'Linewidth', 2)
xlabel('Time (sec)', 'interpreter', 'latex')
ylabel('$\dot{\theta}$ (rad/s)', 'interpreter', 'latex')
legend({'GP-MRAC Output (OL)','Reference Model'},'interpreter', 'latex', 'location', 'best')
hold off

figure(3)
hold on
ss2 = fill(fillT, sigma1Fill, 'y');
alpha(ss2, 0.5)
plot(T, delta_GP_vec, 'b', 'Linewidth', 2)
plot(T, delta_true_vec, 'r--', 'Linewidth', 2)
xlabel('Time (sec)', 'interpreter', 'latex')
ylabel('$\Delta x$', 'interpreter', 'latex')
plots=get(gca, 'Children');
legend(plots([2, 1, 3]), {'Estimated $\Delta (x)$ (OL)','True $\Delta x$', '$\Delta_{est} (x) \pm 2\sigma$'}, 'interpreter', 'latex', 'location', 'best')
hold off

%Wing-rock dynamics one-step simulation, using the Euler Maruyama method
%for Stochastic Differential Equations (SDEs)
function [res, dx] = wrd_sim(x, nu, omega_n, dt)

    mean_x = 0.8 + 0.2314*x(1) + 0.6918*x(2) - 0.6245*abs(x(1))*x(2) + ...
        0.0095*abs(x(2))*x(2) + 0.0214*x(1)^3;
    mu = [x(2); nu + mean_x];
    sigma = [0; sqrt(omega_n)];
    dW = randn(1,1)*sqrt(dt);
    res = x + mu*dt + sigma*dW;
    dx = mean_x + sqrt(omega_n)*dW/sqrt(dt);

end

%Kalman filter for estimating state acceleration (angular acceleration)
function [x_est, P_est] = KF_xdd(z_meas, x_est_p, P_est_p, A, H, R, Q)

    x_pred = A*x_est_p;
    P_pred = A*P_est_p*A' + Q;
    
    Kk = (P_pred*H')/(H*P_pred*H' + R);
    x_est = x_pred + Kk*(z_meas - H*x_pred);
    P_est = P_pred - (P_pred*H')*Kk';

end

%Gaussian Process Regression in which the covariance matrices are updated
%and the mean and covariance of the current estimated state (Gaussian Process) are obtained.
function [mean_predict, cov_predict, KZtau, Ctau] = gpr_full(xpred, xtrain, ytrain, noisecov, kernel_Handle)
    KZtau = build_covariance(xtrain, xtrain, kernel_Handle);
    
    %Ctau = inv(KZtau + sqrt(noisecov)*eye(size(xtrain,1)));
    Ctau = inv(KZtau + (noisecov^2)*eye(size(xtrain,1)));
    
    mean_coeff = Ctau*ytrain;
    vec_pred = build_covariance(xpred, xtrain, kernel_Handle);
    mean_predict = vec_pred*mean_coeff;
    
    cov_predict_pre = build_covariance(xpred,xpred,kernel_Handle);
    cov_predict_up = vec_pred*Ctau*vec_pred';
    cov_predict = cov_predict_pre-cov_predict_up;
end

%Gaussian Process Regression in which only the mean and covariance of the requested states are updated.
%Previously obtained covariance matrices are used since no other data
%points have been added to the training set.
function [mean_predict, cov_predict] = gpr_update(xpred, xtrain, ytrain, Ctau, kernel_Handle)
    mean_coeff = Ctau*ytrain;
    vec_pred = build_covariance(xpred, xtrain, kernel_Handle);
    mean_predict = vec_pred*mean_coeff;
    
    cov_predict_pre = build_covariance(xpred,xpred,kernel_Handle);
    cov_predict_up = vec_pred*Ctau*vec_pred';
    cov_predict = cov_predict_pre-cov_predict_up;
end

%Kernel Linear Independence (KLI) test algorithm used to determine whether the
%new data at state xpred should be added to the current training data set.
%If gamma is greater than a threshold value, the data point is added.
function gamma = KLI_test(xpred, xtrain, KZtau, kernel_Handle)
    
    vec_pred = build_covariance(xpred, xtrain, kernel_Handle);
    cov_predict_pre = build_covariance(xpred,xpred,kernel_Handle);
    gamma = cov_predict_pre - vec_pred*(KZtau\vec_pred');

end

%Function for building covariance matrices between two sets of datapoints
%and a kernel function represented by kern_Handle
function out = build_covariance(x,xp,kern_Handle)
    out = zeros(size(x,1),size(xp,1));
    for ii = 1:size(x,1)
        for jj = 1:size(xp,1)
            out(ii,jj) = kern_Handle(x(ii,1),xp(jj,1))*kern_Handle(x(ii,2),xp(jj,2));
        end
    end
end