%% Select working directory interactively and add to the MATLAB path
%selpath = uigetdir();
%addpath(genpath(selpath));

% Load temperature file
load('WeatherTemp.mat');

% Load Einput files 
Eidata = importdata('Einput7.txt');

% Specify alpha values to adjust enzyme activity levels here 
global Vrubusco_adj;
Vrubusco_adj = 1.36; 
global VmaxAdj;
VmaxAdj = 1.12;
global pcfactor;
pcfactor = 1;

% Adjust enzymes to calibrated enzyme activity levels (Vmax) for rice
Eidata.data(1,1) = Eidata.data(1,1) * Vrubusco_adj;
Eidata.data(2:26,1) = Eidata.data(2:26,1) * VmaxAdj;
Ei = Eidata.data;

% Ensure double-counted enzymes have the same activity i.e. V8 = V5 and V10 = V7
Ei(7) = Ei(4); 
Ei(9) = Ei(6);

% Define the folder containing the BestMatrix text files
results_folder = fullfile('Results', 'Enzymes');

% List all text files inside the folder
all_file_list = dir(fullfile(results_folder, '*.txt'));

% Load average optimized Vmax (1-10) modelled under Cc = 130, 250 and 360 to create stacked enzyme scenarios

% Get the list of 10 text files containing enzyme Vmax of optimizations at
% 130/250/360 umol mol-1 using a regular expression to match names containing _Cc_ and numbered suffixes = 1-9 or 10
optimized_130_list = all_file_list(~cellfun('isempty', regexp({all_file_list.name}, '^outputenz_130_(10|[1-9])\.txt$')));
optimized_250_list = all_file_list(~cellfun('isempty', regexp({all_file_list.name}, '^outputenz_250_(10|[1-9])\.txt$')));
optimized_360_list = all_file_list(~cellfun('isempty', regexp({all_file_list.name}, '^outputenz_360_(10|[1-9])\.txt$')));

% Initialize empty cell arrays to store the matrices
optimized_130_matrices = cell(numel(optimized_130_list), 1);
optimized_250_matrices = cell(numel(optimized_250_list), 1);
optimized_360_matrices = cell(numel(optimized_360_list), 1);

% Loop over the files to retrieve data
for i = 1:numel(optimized_130_list)
    % Read the current file using readmatrix
    optimized_130_matrix = readmatrix(fullfile(results_folder, optimized_130_list(i).name));
    optimized_250_matrix = readmatrix(fullfile(results_folder, optimized_250_list(i).name));
    optimized_360_matrix = readmatrix(fullfile(results_folder, optimized_360_list(i).name));
    % Store the matrix in the cell array
    optimized_130_matrices{i} = optimized_130_matrix;
    optimized_250_matrices{i} = optimized_250_matrix;
    optimized_360_matrices{i} = optimized_360_matrix;
end

% Combine all cells in the array 
combined_130_matrix = horzcat(optimized_130_matrices{:,1});
combined_250_matrix = horzcat(optimized_250_matrices{:,1});
combined_360_matrix = horzcat(optimized_360_matrices{:,1});
blank_rows = zeros(1,10);

% Pad blank rows to align with enzymes
combined_130_matrix_new = vertcat(combined_130_matrix(3:8,:), ...
    blank_rows,combined_130_matrix(9,:),blank_rows,combined_130_matrix(10:11,:),blank_rows,combined_130_matrix(12:25,:));

combined_250_matrix_new = vertcat(combined_250_matrix(3:8,:), ...
    blank_rows,combined_250_matrix(9,:),blank_rows,combined_250_matrix(10:11,:),blank_rows,combined_250_matrix(12:25,:));

combined_360_matrix_new = vertcat(combined_360_matrix(3:8,:), ...
    blank_rows,combined_360_matrix(9,:),blank_rows,combined_360_matrix(10:11,:),blank_rows,combined_360_matrix(12:25,:));

% Calculate averages across the columns for each optimisation
optimised_Vmax_130 = mean(combined_130_matrix_new,2);
optimised_Vmax_250 = mean(combined_250_matrix_new,2);
optimised_Vmax_360 = mean(combined_360_matrix_new,2);

% Generate 25% increase for comparison of enzyme stacking for Rubisco
Vmax_plus_25 = Ei*1.25; % Yoon et al 2020

% Generate twofold increase for comparison of enzyme stacking for all other enzymes
twofold_Vmax = Ei*2; 

%%% Create Vmax profiles to model strategies

% To model low Cc under drought, optimise only Rubisco and SBPase at 130 ppm, keep all other enzyme levels non-optimized
Eidata_low = vertcat(optimised_Vmax_130(1),Ei(2:6),Ei(4),optimised_Vmax_130(8),Ei(6),Ei(10:end));
Eidata_low_2 = vertcat(Vmax_plus_25(1),Ei(2:6),Ei(4),twofold_Vmax(8),Ei(6),Ei(10:end));

% To model ambient Cc, optimise Rubisco, SBPase, Aldolase and PRK at Cc = 250 ppm
Eidata_ambient = vertcat(optimised_Vmax_250(1),Ei(2:3),optimised_Vmax_250(4),Ei(5:6),optimised_Vmax_250(4),optimised_Vmax_250(8),Ei(6),optimised_Vmax_250(10),Ei(11:end));
Eidata_ambient_2 = vertcat(Vmax_plus_25(1),Ei(2:3),twofold_Vmax(4),Ei(5:6),twofold_Vmax(4),twofold_Vmax(8),Ei(6),twofold_Vmax(10),Ei(11:end));

% To model elevated Cc, optimise Rubisco, SBPase, Aldolase, PRK, FBPase and TK at Cc = 360 ppm
Eidata_elevated = vertcat(optimised_Vmax_360(1),Ei(2:3),optimised_Vmax_360(4:6),optimised_Vmax_360(4),optimised_Vmax_360(8),optimised_Vmax_360(6),optimised_Vmax_360(10),Ei(11:end));
Eidata_elevated_2 = vertcat(Vmax_plus_25(1),Ei(2:3),twofold_Vmax(4:6),twofold_Vmax(4),twofold_Vmax(8),twofold_Vmax(6),twofold_Vmax(10),Ei(11:end)); 

Einput=ones(37,1); % 
PPFDi = 2000; % Set light intensity
GRNC=0;

% Create output vectors for A to loop over multiple Cc's
GrossAssimilationRate = zeros(26,1);
GrossAssimilationRate_low = zeros(26,1);
GrossAssimilationRate_ambient = zeros(26,1);
GrossAssimilationRate_elevated = zeros(26,1);
GrossAssimilationRate_low_2 = zeros(26,1);
GrossAssimilationRate_ambient_2 = zeros(26,1);
GrossAssimilationRate_elevated_2 = zeros(26,1);

% Loop function to calculate gross assimilation rate A for a range of Cc values using set of enzyme activity levels 
for i = 1:26 % No. of A values
    CO2i = (130:10:380)'; % Set Cc values ranging from 130-380 with stepsize of 10
    % Calculate assimilation if using Vmax obtained from optimising
    % e-Photosynthesis for a given no. of enzymes (low = 2, ambient = 4, elevated = 6)
    GrossAssimilationRate(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Ei);%0
    GrossAssimilationRate_low(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Eidata_low);%2
    GrossAssimilationRate_ambient(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Eidata_ambient);%4
    GrossAssimilationRate_elevated(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Eidata_elevated);%6

    % Calculate assimilation if using twofold increases in Vmax
    GrossAssimilationRate_low_2(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Eidata_low_2);%2
    GrossAssimilationRate_ambient_2(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Eidata_ambient_2);%4
    GrossAssimilationRate_elevated_2(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Eidata_elevated_2);%6
end

% Save assimilation rates to output text files
writematrix(GrossAssimilationRate,'Results\non_optimized_A.txt','Delimiter','space');
writematrix(GrossAssimilationRate_low,'Results\low_A.txt','Delimiter','space');
writematrix(GrossAssimilationRate_ambient,'Results\ambient_A.txt','Delimiter','space');
writematrix(GrossAssimilationRate_elevated,'Results\elevated_A.txt','Delimiter','space');
writematrix(GrossAssimilationRate_low_2,'Results\low2_A.txt','Delimiter','space');
writematrix(GrossAssimilationRate_ambient_2,'Results\ambient2_A.txt','Delimiter','space');
writematrix(GrossAssimilationRate_elevated_2,'Results\elevated2_A.txt','Delimiter','space');