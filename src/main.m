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


