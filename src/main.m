clc;
clear;
close all; 

%% Importing data from csv file
%filename = fullfile('..','data', 'group_08.csv');
filename = fullfile('..','data', 'group_07.csv');
opts = detectImportOptions(filename, ...
    'Delimiter', ';', ...
    'DecimalSeparator', ',');
T = readtable(filename, opts);

t       = T{:,1};
q1      = T{:,2};
h1      = T{:,3};
h2      = T{:,4};
q2_meas = T{:,5};

%% Identification of model parameters
% Noise parameters
noise_only = q2_meas(1:50);         
E_hat = mean(noise_only);           % Estimate of noise mean value
sig_hat =  var(noise_only);         % Estimate of noise variance
noise_std  = std(noise_only);
fprintf('=============== Noise parameters ==============\n'); 
fprintf('Noise mean value = %.6f\n', E_hat); 
fprintf('Noise variance = %.6f\n', sig_hat);
fprintf('Noise standard deviation = %.6f\n', noise_std);


% 'K' parameter identification (LS)
idx = 51:length(t);          
sqrt_h2_q2 = sqrt(max(h2(idx), 0));
numerator   = sum(sqrt_h2_q2 .* q2_meas(idx));
denominator = sum(sqrt_h2_q2.^2);
k_hat = numerator / denominator;
fprintf('=============== Model parameters ==============\n'); 
fprintf('K_hat = %.6e\n', k_hat);

% 'α', 'S_2' parameters identification (LS)
idx2 = 52:(length(t)-1);
dh2 = (h2(idx2+1) - h2(idx2-1)) / 0.5;
Y = dh2(:);
sqrt_h2_dh2 = sqrt(max(h2(idx2), 0));
Phi = [h1(idx2)-h2(idx2), -sqrt_h2_dh2];
Theta_hat = (Phi' * Phi) \ (Phi' * Y);
a1_hat = Theta_hat(1);
a2_hat = Theta_hat(2);

S2_hat = k_hat / a2_hat;
alpha_hat = a1_hat * S2_hat;

fprintf('S2_hat = %.6e\n', S2_hat);
fprintf('alpha_hat = %.6e\n', alpha_hat);

% 'S_1' parameter identification (LS)
dh1 = (h1(idx2+1) - h1(idx2-1)) / 0.5;
Y = dh1(:);
Psi = q1(idx2) - alpha_hat * (h1(idx2) - h2(idx2));
Psi = Psi(:);
a3_hat = (Psi' * Y) / (Psi' * Psi);
S1_hat = 1 / a3_hat;
fprintf('S1_hat = %.6e\n', S1_hat);

% Model built on estimated parameters 
q2_model  = k_hat * sqrt_h2_q2;
dh2_model = Phi * Theta_hat;
dh1_model = a3_hat * Psi;

%plot_identification_results(t, idx, idx2, h1, h2, q2_meas, sqrt_h2_q2, q2_model, dh2, dh2_model, dh1, dh1_model, Psi)

%% Extended Kalman filter implementation - Task I
% Initialization
Ts = t(2)-t(1);
N = length(t);

Q = diag([1e-6, 1e-6]);     % Process noise covariance
R = sig_hat;                % Measurment noise covariance
x_hat = [h1(1), h2(1)];     % Initial state estimate
P = diag([1e-4, 1e-4]);     % Initial covariance estimate

% Memory allocation
x_hat_hist = zeros(2,N);
x_pred_hist = zeros(2,N);
y_pred_hist = zeros(2,N);
innovation_hist = zeros(1,N);

x_hat_hist(:,1) = x_hat;

for k = 1:N-1
    
    % Current corrected state estimate
    x1 = x_hat(1);
    x2 = max(x_hat(2),0);
    
    % Known nput 
    u = q1(k);

    % State prediction 
    x_pred = [
        x1 + Ts * ((1 / S1_hat) * u - (alpha_hat / S1_hat) * (x1 - x2));
        x2 + Ts * ((alpha_hat / S2_hat) * (x1 - x2) - (k_hat / S2_hat) * sqrt(x2))
    ];

    % Calculation Jacobian of the state model
    A = [
        1 - Ts * (alpha_hat / S1_hat), Ts * (alpha_hat / S1_hat);
        Ts * (alpha_hat / S2_hat), 1 - Ts * (alpha_hat / S2_hat) - Ts * (k_hat / (2 * S2_hat * sqrt(max(x2, 1e-9))))
    ];

    % Covariance prediction
    P_pred = A * P * A' + Q;

    % Measurement prediciton 
    x2_pred = max(x_pred(2), 0);
    y_pred = k_hat * sqrt(x2_pred);
    
    % Calculation Jacobian of the measurement model
    C = [0, k_hat / (2 * sqrt(max(x2_pred, 1e-9)))];
    
    % Innovation 
    innovation = q2_meas(k+1) - y_pred;

    % Innovation covariance
    S = C * P_pred * C' + R;

    % Kalman gain calculation
    K = P_pred * C' / S;  
    
    % State correction
    x_hat = x_pred + K * innovation;

    % Covariance correction 
    P = P_pred - K * S * K';

    % Store results 
    x_hat_hist(:,k+1) = x_hat;
    x_pred_hist(:,k+1) = x_pred;
    y_pred_hist(k+1) = y_pred;
    innovation_hist(k+1) = innovation;

end

plot_ekf_results(t, h1, h2, q2_meas, x_hat_hist, y_pred_hist, innovation_hist, k_hat);