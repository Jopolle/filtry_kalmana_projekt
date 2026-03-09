clc;
clear;
close all; 

%% Zczytywanie danych z pliku csv do tabeli 
filename = fullfile('..','data', 'group_08.csv');
opts = detectImportOptions(filename, ...
    'Delimiter', ';', ...
    'DecimalSeparator', ',');
T = readtable(filename, opts);

t       = T{:,1};
q1      = T{:,2};
h1      = T{:,3};
h2      = T{:,4};
q2_meas = T{:,5};

%% Identyfikacja parametrów modelu 
noise_only = q2_meas(1:50);         % Część czysto szumowa
E_hat = mean(noise_only);           % Estymata średniej szumu pomiarowego
sig_hat =  var(noise_only);         % Estymata wariancji szumu pomiarowego
noise_std  = std(noise_only);
fprintf('Średnia szumu = %.6f\n', E_hat); 
fprintf('Wariancja szumu = %.6f\n', sig_hat);
fprintf('Odchylenie standardowe szumu = %.6f\n', noise_std);
