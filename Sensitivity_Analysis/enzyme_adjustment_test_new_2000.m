%% Sensitivity Analysis for determining assimilation rates associated with increases 
%% from non-optimal enzyme increases (sampled between no change from non-opt to random FC between 1-1.25 for Rubisco or 1-4 for other enzymes) 

% Get absolute path of this script
ScriptPath = mfilename('fullpath');

% Locate directory containing this script
ScriptDir = fileparts(ScriptPath);

% Change directory to main e-Photosynthesis code repository
cd(fullfile(ScriptDir,'..'));
ePhotosynthesis_repository = pwd;

% Add repository to path
addpath(genpath(ePhotosynthesis_repository));

% Load original enzyme activity levels for e-Photosynthesis (pre-adjustment)
Eidata = importdata('Einput7.txt');
Eio = Eidata.data(:,1);
load('WeatherTemp.mat');

% Specify scaling factors obtained from running Jmax_adj_simple and Vcmax_adj_simple
global Vrubusco_adj;
Vrubusco_adj = 1.36; 
global VmaxAdj;%
VmaxAdj = 1.12;

% Adjust Einput for rice according to scaling factors (pre-optimisation)
Eio(1) = Eidata.data(1,1)*Vrubusco_adj;
Eio(2:26) = Eidata.data(2:26,1)*VmaxAdj;

% Ensure double-counted enzymes have the same activity i.e. V8=V5 and V10=V7
Eio(7) = Eio(4); 
Eio(9) = Eio(6);

%% Set lower bounds for each enzyme based on FC = 1 from non-opt (i.e. no change in protein from WT)
rubisco_lb = 1;     % FC = 1
aldolase_lb = 1;    % FC = 1
sbpase_lb = 1;      % FC = 1
prk_lb = 1;         % FC = 1
fbpase_lb = 1;      % FC = 1
tk_lb = 1;          % FC = 1

%% Set upper bounds for each enzyme based on highest FC values found in literature from non-opt 
rubisco_ub = 1.25;  % FC = 1.25, Yoon et al (2020) - sense line is 125% of WT
aldolase_ub = 2.8;    % FC = 2.8, Simkin 2015
sbpase_ub = 2.8;      % FC = 2.8,
prk_ub = 2.8;         % FC = 2.8,
fbpase_ub = 2.8;      % FC = 2.8,
tk_ub = 2.8;          % FC = 2.8,

%% Use a random number generator to randomise FC values by which to increase protein of Rubisco Aldolase SBPase and PRK
rng("shuffle");
rng_name = rng;

%% Using rand, start to draw random scalars from a uniform distribution for enzyme fold changes (test with 100 and then increase to 2000)
%% Rubisco SBPase
% RS Rubisco fold change (drawn between 1 and rubisco_250_ub)
RS_rubisco_FC = (rubisco_ub-rubisco_lb).*abs(rand(2000,1)) + rubisco_lb;
RS_rubisco_FC(1) = 1; % substitute first value with fold change 1

% Check random vector of RS Rubisco FC is between limits
RS_rubisco_range = [min(RS_rubisco_FC) max(RS_rubisco_FC)];

% RS SBPase fold change (drawn between 1 and RS_sbpase_ub)
RS_sbpase_FC = (sbpase_ub-sbpase_lb).*abs(rand(2000,1)) + sbpase_lb;
RS_sbpase_FC(1) = 1; % substitute first value with fold change 1

% Check random vector of RS SBPase FC is between limits
RS_sbpase_range = [min(RS_sbpase_FC) max(RS_sbpase_FC)];

%% Rubisco-SBPase-Aldolase-PRK
% RSAP Rubisco fold change (drawn between 1 and rubisco_250_ub)
RSAP_rubisco_FC = (rubisco_ub-rubisco_lb).*abs(rand(2000,1)) + rubisco_lb;
RSAP_rubisco_FC(1) = 1; % substitute first value with fold change 1

% Check random vector of RSAP Rubisco FC is between limits
RSAP_rubisco_range = [min(RSAP_rubisco_FC) max(RSAP_rubisco_FC)];

% RSAP SBPase fold change (drawn between 1 and sbpase_250_ub)
RSAP_sbpase_FC = (sbpase_ub-sbpase_lb).*abs(rand(2000,1)) + sbpase_lb;
RSAP_sbpase_FC(1) = 1; % substitute first value with fold change 1

% Check random vector of RSAP SBPase FC is between limits
RSAP_sbpase_range = [min(RSAP_sbpase_FC) max(RSAP_sbpase_FC)];

% RSAP Aldolase fold changes (drawn between 1 and aldolase_250_ub)
RSAP_aldolase_FC = (aldolase_ub-aldolase_lb).*abs(rand(2000,1)) + aldolase_lb;
RSAP_aldolase_FC(1) = 1; % substitute first value with fold change 1 

% Check random vector of RSAP aldolase FC is between limits
RSAP_aldolase_range = [min(RSAP_aldolase_FC) max(RSAP_aldolase_FC)];

% RSAP PRK fold change (drawn between 1 and prk_250_ub)
RSAP_prk_FC = (prk_ub-prk_lb).*abs(rand(2000,1)) + prk_lb;
RSAP_prk_FC(1) = 1; % substitute first value with fold change 1

% Check random vector of RSAP PRK FC is between limits
RSAP_prk_range = [min(RSAP_prk_FC) max(RSAP_prk_FC)];

%% Rubisco-SBPase-Aldolase-PRK-FBPase-TK
% RSAPFT Rubisco fold change (drawn between 1 and rubisco_360_ub)
RSAPFT_rubisco_FC = (rubisco_ub-rubisco_lb).*abs(rand(2000,1)) + rubisco_lb;
RSAPFT_rubisco_FC(1) = 1; % substitute first value with fold change 1

% Check random vector of RSAPFT Rubisco FC is between limits
RSAPFT_rubisco_range = [min(RSAPFT_rubisco_FC) max(RSAPFT_rubisco_FC)];

% RSAPFT SBPase fold change (drawn between 1 and sbpase_360_ub)
RSAPFT_sbpase_FC = (sbpase_ub-sbpase_lb).*abs(rand(2000,1)) + sbpase_lb;
RSAPFT_sbpase_FC(1) = 1; % substitute first value with fold change 1

% Check random vector of RSAPFT SBPase FC is between limits
RSAPFT_sbpase_range = [min(RSAPFT_sbpase_FC) max(RSAPFT_sbpase_FC)];

% RSAPFT Aldolase fold changes (drawn between 1 and aldolase_360_ub)
RSAPFT_aldolase_FC = (aldolase_ub-aldolase_lb).*abs(rand(2000,1)) + aldolase_lb;
RSAPFT_aldolase_FC(1) = 1; % substitute first value with fold change 1 

% Check random vector of RSAPFT aldolase FC is between limits
RSAPFT_aldolase_range = [min(RSAPFT_aldolase_FC) max(RSAPFT_aldolase_FC)];

% RSAPFT PRK fold change (drawn between 1 and prk_360_ub)
RSAPFT_prk_FC = (prk_ub-prk_lb).*abs(rand(2000,1)) + prk_lb;
RSAPFT_prk_FC(1) = 1; % substitute first value with fold change 1

% Check random vector of RSAPFT PRK FC is between limits
RSAPFT_prk_range = [min(RSAPFT_prk_FC) max(RSAPFT_prk_FC)];

%RSAPFT FBPase fold changes (drawn between 1 and fbpase_360_ub)
RSAPFT_fbpase_FC = (fbpase_ub-fbpase_lb).*abs(rand(2000,1)) + fbpase_lb;
RSAPFT_fbpase_FC(1) = 1; % substitute first value with fold change 1

% Check random vector of RSAPFT FBPase FC is between limits
RSAPFT_fbpase_range = [min(RSAPFT_fbpase_FC) max(RSAPFT_fbpase_FC)];

% RSAPFT TK fold changes (drawn between 1 and tk_360_ub)
RSAPFT_tk_FC = (tk_ub-tk_lb).*abs(rand(2000,1)) + tk_lb;
RSAPFT_tk_FC(1) = 1; % substitute first value with fold change 1

% Check random vector of RSAPFT TK FC is between limits
RSAPFT_tk_range = [min(RSAPFT_tk_FC) max(RSAPFT_tk_FC)];
%% Pre-allocate outputs for percentage change
[percent_change_RS_rubisco,...
 percent_change_RS_sbpase, ...
 percent_change_RSAP_rubisco,...
 percent_change_RSAP_aldolase,...
 percent_change_RSAP_sbpase, ...
 percent_change_RSAP_prk,...
 percent_change_RSAPFT_rubisco,...
 percent_change_RSAPFT_aldolase,...
 percent_change_RSAPFT_sbpase, ...
 percent_change_RSAPFT_prk,...
 percent_change_RSAPFT_fbpase,...
 percent_change_RSAPFT_tk]= deal(zeros(2000,1));
%% Compute percentage changes from each FC distribution
for p = 1:length(percent_change_RS_rubisco)
percent_change_RS_rubisco(p) = (RS_rubisco_FC(p)-RS_rubisco_FC(1))/RS_rubisco_FC(1)*100;
percent_change_RS_sbpase(p) = (RS_sbpase_FC(p)-RS_sbpase_FC(1))/RS_sbpase_FC(1)*100;
%
percent_change_RSAP_rubisco(p) = (RSAP_rubisco_FC(p)-RSAP_rubisco_FC(1))/RSAP_rubisco_FC(1)*100;
percent_change_RSAP_aldolase(p) = (RSAP_aldolase_FC(p)-RSAP_aldolase_FC(1))/RSAP_aldolase_FC(1)*100;
percent_change_RSAP_sbpase(p) = (RSAP_sbpase_FC(p)-RSAP_sbpase_FC(1))/RSAP_sbpase_FC(1)*100;
percent_change_RSAP_prk(p) = (RSAP_prk_FC(p)-RSAP_prk_FC(1))/RSAP_prk_FC(1)*100;
%
percent_change_RSAPFT_rubisco(p) = (RSAPFT_rubisco_FC(p)-RSAPFT_rubisco_FC(1))/RSAPFT_rubisco_FC(1)*100;
percent_change_RSAPFT_aldolase(p) = (RSAPFT_aldolase_FC(p)-RSAPFT_aldolase_FC(1))/RSAPFT_aldolase_FC(1)*100;
percent_change_RSAPFT_sbpase(p) = (RSAPFT_sbpase_FC(p)-RSAPFT_sbpase_FC(1))/RSAPFT_sbpase_FC(1)*100;
percent_change_RSAPFT_prk(p) = (RSAPFT_prk_FC(p)-RSAPFT_prk_FC(1))/RSAPFT_prk_FC(1)*100;
percent_change_RSAPFT_fbpase(p) = (RSAPFT_fbpase_FC(p)-RSAPFT_fbpase_FC(1))/RSAPFT_fbpase_FC(1)*100;
percent_change_RSAPFT_tk(p) = (RSAPFT_tk_FC(p)-RSAPFT_tk_FC(1))/RSAPFT_tk_FC(1)*100;
end
%% Pre-allocate outputs for increased Vmax
[increased_vmax_RS_rubisco,...
 increased_vmax_RS_sbpase,...
 increased_vmax_RSAP_rubisco,...
 increased_vmax_RSAP_sbpase,...
 increased_vmax_RSAP_aldolase,...
 increased_vmax_RSAP_prk,...
 increased_vmax_RSAPFT_rubisco,...
 increased_vmax_RSAPFT_sbpase,...
 increased_vmax_RSAPFT_aldolase,...
 increased_vmax_RSAPFT_prk,...
 increased_vmax_RSAPFT_fbpase,...
 increased_vmax_RSAPFT_tk] = deal(zeros(2000,1));
%% Starting from non-optimized activity, increase all Vmax values by new FC of non-optimized Vmax 
for j=1:length(increased_vmax_RS_rubisco)
increased_vmax_RS_rubisco(j) = (RS_rubisco_FC(j)*Eio(1))';
increased_vmax_RS_sbpase(j) = (RS_sbpase_FC(j)*Eio(8))';
increased_vmax_RSAP_rubisco(j) = (RSAP_rubisco_FC(j)*Eio(1))';
increased_vmax_RSAP_aldolase(j) = (RSAP_aldolase_FC(j)*Eio(4))';
increased_vmax_RSAP_sbpase(j) = (RSAP_sbpase_FC(j)*Eio(8))';
increased_vmax_RSAP_prk(j) = (RSAP_prk_FC(j)*Eio(10))';
increased_vmax_RSAPFT_rubisco(j) = (RSAPFT_rubisco_FC(j)*Eio(1))';
increased_vmax_RSAPFT_aldolase(j) = (RSAPFT_aldolase_FC(j)*Eio(4))';
increased_vmax_RSAPFT_sbpase(j) = (RSAPFT_sbpase_FC(j)*Eio(8))';
increased_vmax_RSAPFT_prk(j) = (RSAPFT_prk_FC(j)*Eio(10))';
increased_vmax_RSAPFT_fbpase(j) = (RSAPFT_fbpase_FC(j)*Eio(5))';
increased_vmax_RSAPFT_tk(j) = (RSAPFT_tk_FC(j)*Eio(6))';
end

% Transpose variables
increased_vmax_RS_rubisco = increased_vmax_RS_rubisco';
increased_vmax_RS_sbpase = increased_vmax_RS_sbpase';
increased_vmax_RSAP_rubisco = increased_vmax_RSAP_rubisco';
increased_vmax_RSAP_aldolase = increased_vmax_RSAP_aldolase';
increased_vmax_RSAP_sbpase = increased_vmax_RSAP_sbpase';
increased_vmax_RSAP_prk = increased_vmax_RSAP_prk';
increased_vmax_RSAPFT_rubisco = increased_vmax_RSAPFT_rubisco';
increased_vmax_RSAPFT_aldolase = increased_vmax_RSAPFT_aldolase';
increased_vmax_RSAPFT_sbpase = increased_vmax_RSAPFT_sbpase';
increased_vmax_RSAPFT_prk = increased_vmax_RSAPFT_prk';
increased_vmax_RSAPFT_fbpase = increased_vmax_RSAPFT_fbpase';
increased_vmax_RSAPFT_tk = increased_vmax_RSAPFT_tk';
%% Create 100 (increased to 2000) rows of non-optimized protein 
non_opt_vmax = repmat(Eio,1,2000);
%% Re-specify new list of values for Einput by concatenating combinations of 2 or 4 increased Vmax with non-optimized Vmax for enzymes in Einput vector
%% Non-optimal control 
Einput_nonopt = non_opt_vmax(1:end,1);
%% Optimal Rubisco and SBPase
Einput_RS = vertcat(increased_vmax_RS_rubisco,non_opt_vmax(2:7,1:end),increased_vmax_RS_sbpase,non_opt_vmax(9:end,1:end));
%% Optimal Rubisco, Aldolase, SBPase, PRK % ADD IN Aldolase and PRK
Einput_RSAP = vertcat(increased_vmax_RSAP_rubisco,non_opt_vmax(2:3,1:end),increased_vmax_RSAP_aldolase,non_opt_vmax(5:6,1:end),increased_vmax_RSAP_aldolase,increased_vmax_RSAP_sbpase,non_opt_vmax(9,1:end),increased_vmax_RSAP_prk,non_opt_vmax(11:end,1:end));
%% Optimal Rubisco, Aldolase, SBPase, PRK, FBPase, TK % ADD IN FBPASE AND TK
Einput_RSAPFT = vertcat(increased_vmax_RSAPFT_rubisco,non_opt_vmax(2:3,1:end),increased_vmax_RSAPFT_aldolase,increased_vmax_RSAPFT_fbpase,increased_vmax_RSAPFT_tk,increased_vmax_RSAPFT_aldolase,increased_vmax_RSAPFT_sbpase,increased_vmax_RSAPFT_tk,increased_vmax_RSAPFT_prk,non_opt_vmax(11:end,1:end));
%% Set conditions
Einput=ones(37,1); 
% Set light intensity 
PPFDi = 2000; 
% Set other variables
global pcfactor;
pcfactor = 1;
%% Pre-allocate outputs for assimilation rate for different strategies at low, ambient and elevated CO2
[GrossAssimilationRate_RS_130,...     % two
 GrossAssimilationRate_RSAP_130,...   
 GrossAssimilationRate_RSAPFT_130,... 
 GrossAssimilationRate_RS_250,...     % four
 GrossAssimilationRate_RSAP_250,... 
 GrossAssimilationRate_RSAPFT_250,...
 GrossAssimilationRate_RS_360,...     % six
 GrossAssimilationRate_RSAP_360,... 
 GrossAssimilationRate_RSAPFT_360] = deal(zeros(2000,1));
%% Set Cc value % Query multiple inputs - 130, 250 or 360 - here OR as inputs to a function
CO2i_low = 130;
CO2i_amb = 250; 
CO2i_high = 360;
%% Calculate net assimilation for control - only changing CO2 input
GrossAssimilationRate_nonopt_130 = EPS_Drive_GRNs(Einput,CO2i_low,PPFDi,WeatherTemp,0,0,Einput_nonopt);
GrossAssimilationRate_nonopt_250 = EPS_Drive_GRNs(Einput,CO2i_amb,PPFDi,WeatherTemp,0,0,Einput_nonopt);
GrossAssimilationRate_nonopt_360 = EPS_Drive_GRNs(Einput,CO2i_high,PPFDi,WeatherTemp,0,0,Einput_nonopt);
%GrossAssimilationRate_opt = EPS_Drive_GRNs(Einput,CO2i,PPFDi,WeatherTemp,0,0,Einput_opt);
%% Calculate net assimilation rate A for increased enzyme Vmax (100 results per 1xquad/6xduo enzymes) - changing CO2 and enzymes in Einput
% Loop and calculate assimilation rates for combinations of increased Vmax enzymes
for n=1:length(GrossAssimilationRate_RS_130)
GrossAssimilationRate_RS_130(n) = EPS_Drive_GRNs(Einput,CO2i_low,PPFDi,WeatherTemp,0,0,Einput_RS(:,n));
GrossAssimilationRate_RSAP_130(n) = EPS_Drive_GRNs(Einput,CO2i_low,PPFDi,WeatherTemp,0,0,Einput_RSAP(:,n));
GrossAssimilationRate_RSAPFT_130(n) = EPS_Drive_GRNs(Einput,CO2i_low,PPFDi,WeatherTemp,0,0,Einput_RSAPFT(:,n));
GrossAssimilationRate_RS_250(n) = EPS_Drive_GRNs(Einput,CO2i_amb,PPFDi,WeatherTemp,0,0,Einput_RS(:,n));
GrossAssimilationRate_RSAP_250(n) = EPS_Drive_GRNs(Einput,CO2i_amb,PPFDi,WeatherTemp,0,0,Einput_RSAP(:,n));
GrossAssimilationRate_RSAPFT_250(n) = EPS_Drive_GRNs(Einput,CO2i_amb,PPFDi,WeatherTemp,0,0,Einput_RSAPFT(:,n));
GrossAssimilationRate_RS_360(n) = EPS_Drive_GRNs(Einput,CO2i_high,PPFDi,WeatherTemp,0,0,Einput_RS(:,n));
GrossAssimilationRate_RSAP_360(n) = EPS_Drive_GRNs(Einput,CO2i_high,PPFDi,WeatherTemp,0,0,Einput_RSAP(:,n));
GrossAssimilationRate_RSAPFT_360(n) = EPS_Drive_GRNs(Einput,CO2i_high,PPFDi,WeatherTemp,0,0,Einput_RSAPFT(:,n));
end
%% On HEC
% Save the work space 
% save(workspacefileName);
save(fullfile("Sensitivity_Analysis","output_enzyme_adjustment_test_2000_new.mat"));