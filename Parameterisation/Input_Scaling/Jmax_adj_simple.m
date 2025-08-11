%% Simplified Jmax_adj for running only the metabolic model (MetaOnly==1)
%% Select working directory interactively and add to the MATLAB path
selpath = uigetdir();
addpath(genpath(selpath));
%% Load parameters
parameters = readtable("CO2_response_fitting/Outputs/complete_Fits_rice_params_all.csv",'PreserveVariableNames', true);
Gstar_vals = readtable("CO2_response_fitting/Outputs/Gstars_rice",'PreserveVariableNames', true);
Kcair_vals = readtable("CO2_response_fitting/Outputs/Kcair_umol_rice",'PreserveVariableNames', true);
temp_vals = readtable("CO2_response_fitting/Outputs/temps_rice",'PreserveVariableNames', true);
Gm_vals = readtable("CO2_response_fitting/Outputs/gm_umol_rice",'PreserveVariableNames', true);

Kcair_vals(:,1) = [];  % Remove the first column
Gm_vals(:,1) = [];  
temp_vals(:,1) = []; 

Kcair_vals.Properties.VariableNames = arrayfun(@num2str, 1:8, 'UniformOutput', false);  
Gm_vals.Properties.VariableNames = arrayfun(@num2str, 1:8, 'UniformOutput', false);  

Vcmax_m = mean(parameters.VcMax);
J = mean(parameters.J);
Gr = mean(Gstar_vals.x);
Rd = mean(parameters.rL);
Gm = mean(Gm_vals{:,:}); 
Kc_air = mean(Kcair_vals{:,:}); 

Lii = 2000;%Light intensity from IRRI
O = 210;%mbar

%%%%%%%%%%%%%%%%%%%%%
% Establish range of Cc as inputs
CA = [400,600,800,1000];

% Before beginning, we need to replot Ci vs A and 
% substitute the mean fitted params Vcmax, J and TPU in the Farq model equations 
% to ascertain where CO2 ranges of Ac, Aj and Ap lie 
% Transition point between RuBP and TPU limitation is around Ci = 798 umol
% equivalent to CA = 1140 umol
% Thus the initial range of CA=[800,1200,1500,1800] looks already to be in the TPU limited range
% A better range for scaling RuBP limitation would be CA = [400,600,800,1000]

global Vrubusco_adj;
Vrubusco_adj = 1.0; % keep same rubisco activity as original Einput
global VmaxAdj;%adjust enzyme activity
global pcfactor;  
ProteinTotalRatio=0;
%pcfactor=1/ProteinTotalRatio;
% Change pcfactor to 1 here to avoid Inf
pcfactor = 1;
%%%%%%%%%%%%%%%%%%%%%

Einput = ones(37,1);%No gene expression data input
Edata = importdata('Einput7.txt');
Eio = Edata.data(:,1);
%MetaOnly=1;% if MetaOnly=1 run only Metabolic model
WeatherTemp = mean(temp_vals.x); %Avg Tleaf, Original = 25C
GRNC = 0; % In EPS_Drive_GRNs, if GRNC==0 cATPsyn,CPSi,cNADPHsyn and cpsii=1 are all set to 1

%% Assimilation rates for Farquhar model 
% Other than 4 different Ci values, everything else is the same - 
% Maybe for indexing in later steps its easier to keep the repeated values of assimilation?
Farq_Matrix_J = zeros(200,5);
k = 1; % Initialize row index
for j = 1:50
    VmaxAdj = 0.5 + j * 0.02; % adjust enzyme activity
    Eio(1) = Edata.data(1,1) * Vrubusco_adj;
    Eio(2:27) = Edata.data(2:27,1) * VmaxAdj;
    % Farquhar model A calculation 
    for i = 1:4
        Air_CO2 = CA(i);
        Ci = Air_CO2 * 0.7; % intercellular CO2 
        b = (J/4) - Rd + (Ci + 2*Gr) * Gm;
        c = ((Ci - Gr) * J/4 - (Ci + 2*Gr) * Rd) * Gm;
        % Previous Net Aj calculation
        % ACI_m=J*(Ci-Gr)/(4*Ci+8*Gr)-Rd; 
        % Net Aj expressed as a function of Cc does not include Rd
        Net_A = ((b - sqrt(b^2 - 4*c)) / 2); 
        % Gross_A is being calculated in EPS_Drive_GRNs instead of Net_A
        % Therefore Gross Aj expressed as a function of Cc should include Rd
        Gross_A = Net_A + Rd; % Gross Aj expressed as a function of Cc includes Rd
        % Use Net_A to calculate Cc
        Cc = (Ci - Net_A/Gm); 
        % Fill in results matrix
        Farq_Matrix_J(k, 1) = VmaxAdj;
        Farq_Matrix_J(k, 2) = CA(i);
        Farq_Matrix_J(k, 3) = Ci;
        Farq_Matrix_J(k, 4) = Gross_A;
        Farq_Matrix_J(k, 5) = Cc; %% Cc
        k = k + 1; % Move to the next row
    end
end
%% Assimilation rates for e-Photosynthesis model
ePhoto_Matrix_J = zeros(200,5);
k = 1; % Initialize row index
for j = 1:50
    VmaxAdj = 0.5 + j * 0.02; % Adjust enzyme activity - start at 0.5, then 1.0 and check minima of SSR vals
    Eio(1) = Edata.data(1,1) * Vrubusco_adj;
    Eio(2:27) = Edata.data(2:27,1) * VmaxAdj;
    % e-Photosynthesis model A calculation 
    for i = 1:4
        Air_CO2 = CA(i);
        Ci = Air_CO2 * 0.7;
        Cc = Farq_Matrix_J(k,5); % check index
        PPFDi = Lii; 
        GrossAssimilation = EPS_Drive_GRNs(Einput,Cc,PPFDi,WeatherTemp,GRNC,0,Eio);
        %GrossAssimilation = EPS_Drive_GRNs(Einput,Ci,PPFDi,WeatherTemp,GRNC,0,Eio);
        %Check assimilation rates from using Cc
        ePhoto_Matrix_J(k, 1) = VmaxAdj;
        ePhoto_Matrix_J(k, 2) = CA(i);
        ePhoto_Matrix_J(k, 3) = Ci;
        ePhoto_Matrix_J(k, 4) = GrossAssimilation;
        ePhoto_Matrix_J(k, 5) = Cc; %% Cc
        k = k + 1; % Move to the next row
    end
end

%Create vector to store Vmax_adj, differences in assimilation rates and corresponding SSRs
Diff_Matrix_J = zeros(200,2);

% Compute differences between two photosynthetic models
for k = 1:length(Farq_Matrix_J)
    for j = 1:50
    Diff_Matrix_J(k,1) = Farq_Matrix_J(k,1);
    Diff_Matrix_J(k,2) = (Farq_Matrix_J(k,4)-ePhoto_Matrix_J(k,4))^2;%the squares of the residuals
        while k<200
        k = k + 1; % Move to the next row
        end
    end
end

% Get sums of squared residuals by summing every four rows 
SSR_Matrix_J = zeros(50,2);
VmaxAdj_vals = Diff_Matrix_J(:,1);
SSR_Matrix_J(:,1) = VmaxAdj_vals(4:4:end);
SSR_Matrix_J(:,2) = sum(reshape(Diff_Matrix_J(:,2), 4, []))';

% Get full list of gross assimilation rates associated with each Vmax_Adj
% sampled and for each Cc tested
numbers_J = reshape(ePhoto_Matrix_J(:,4),[4,50]);
%toc

% Find scaling factor with lowest SSR
[~, scaling_index_J] = min(SSR_Matrix_J(:,2));
a_Enzymes = SSR_Matrix_J(scaling_index_J, 1);
min_SSR_J = SSR_Matrix_J(scaling_index_J, 2);

% Save result 
save Jmax_simple_new_result.mat;

% % Plot scaling factors (optional)
% fig = figure;
% scatter(SSR_Matrix_J(:,1),SSR_Matrix_J(:,2),'MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0]);
% hold on
% scatter(a_Enzymes,min_SSR_J,'MarkerEdgeColor','r','MarkerFaceColor','r')
% xticks(0.5:0.2:1.5);
% xlabel('α_{Enzymes}');
% ylabel('SSR');
% 
% % Set figure size
% set(fig, 'PaperUnits', 'inches');
% set(fig, 'PaperPosition', [0 0 6 4]);        % [left bottom width height]
% set(fig, 'PaperSize', [6 4]);                % Exact size of output file
% 
% % Export figure
% print(fig, 'Scaling_Factor_Optimisation_for_RuBP_regeneration.pdf', '-dpdf', '-r300');

