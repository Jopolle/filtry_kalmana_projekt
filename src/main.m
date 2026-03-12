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

plot_identification_results(t, idx, idx2, h1, h2, q2_meas, sqrt_h2_q2, q2_model, dh2, dh2_model, dh1, dh1_model, Psi)