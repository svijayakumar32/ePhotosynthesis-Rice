%% Simplified Vcmax_adj for running only the metabolic model (MetaOnly==1)
%% Select working directory interactively and add to the MATLAB path
% selpath = uigetdir();
% addpath(genpath(selpath));
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
CA = [100,150,200,250,300];
%Transition point between Rubisco and RuBP limitations is around Ci = 210
%umol/mol or 21 Pa which in terms of CA = 300 umol
%So this range is OK
global Vrubusco_adj;
global VmaxAdj;%adjust enzyme activity
VmaxAdj = 1.12;%adjust enzyme activity i.e. sub in optimal Vmaxadj from Jmax_adj % default 1.3
global pcfactor;  
ProteinTotalRatio = 0;
%pcfactor=1/ProteinTotalRatio;
%21/05/24 Change pcfactor to 1 here to avoid Inf
pcfactor = 1;
%%%%%%%%%%%%%%%%%%%%%

Einput = ones(37,1);%No gene expression data input
Edata = importdata('Einput7.txt');
Eio = Edata.data(:,1);
%MetaOnly=1;% if MetaOnly=1 run only Metabolic model
WeatherTemp = mean(temp_vals.x); %Avg Tleaf, Original = 25C
GRNC = 0; % In EPS_Drive_GRNs, if GRNC==0 cATPsyn,CPSi,cNADPHsyn and cpsii=1 are all set to 1

%% Assimilation rates for Farquhar model 
% Other than 5 different Ci values, everything else is the same - 
% Maybe for indexing in later steps its easier to keep the repeated values of assimilation?
Farq_Matrix_V = zeros(250,5);
k = 1; % Initialize row index
for j = 1:50
    Vrubusco_adj = 0.5 + j * 0.02; % adjust enzyme activity
    Eio(1) = Edata.data(1,1) * Vrubusco_adj;
    Eio(2:27) = Edata.data(2:27,1) * VmaxAdj;
    % Farquhar model A calculation 
    for i = 1:5
        Air_CO2 = CA(i);
        Ci = Air_CO2 * 0.7; % intercellular CO2 
        b = Vcmax_m-Rd+(Ci+Kc_air)*Gm;
        c = ((Ci-Gr)*Vcmax_m-(Ci+Kc_air)*Rd)*Gm;
        %Previous Ac calculation, subbing in Kc_air for Kc*(1+O/Ko)
        %ACI_m=Vcmax_m*(Ci-Gr)/(Ci+Kc_air)-Rd; 
        % Net Ac expressed as a function of Cc does not include Rd
        Net_A = ((b - sqrt(b^2 - 4*c)) / 2); 
        % Gross_A is being calculated in EPS_Drive_GRNs instead of Net_A
        % Therefore Gross Ac expressed as a function of Cc should include Rd
        Gross_A = Net_A + Rd; 
        % Use Net_A to calculate Cc
        Cc = (Ci - Net_A/Gm); 
        % Fill in results matrix
        Farq_Matrix_V(k, 1) = Vrubusco_adj;
        Farq_Matrix_V(k, 2) = CA(i);
        Farq_Matrix_V(k, 3) = Ci;
        Farq_Matrix_V(k, 4) = Gross_A;
        Farq_Matrix_V(k, 5) = Cc; %% Cc
        k = k + 1; % Move to the next row
    end
end

%% Assimilation rates for e-Photosynthesis model
ePhoto_Matrix_V = zeros(250,5);
k = 1; % Initialize row index
for j = 1:50
    Vrubusco_adj = 0.5 + j * 0.02; % Adjust enzyme activity - start at 0.5, then 1.0 and check minima of SSR vals
    Eio(1) = Edata.data(1,1) * Vrubusco_adj;
    Eio(2:27) = Edata.data(2:27,1) * VmaxAdj;
    % e-Photosynthesis model A calculation 
    for i = 1:5
        Air_CO2 = CA(i);
        Ci = Air_CO2 * 0.7; % intercellular CO2 
        Cc = Farq_Matrix_V(k,5); % check index
        PPFDi = Lii; 
        GrossAssimilation = EPS_Drive_GRNs(Einput,Cc,PPFDi,WeatherTemp,GRNC,0,Eio);
        % Fill in results matrix
        ePhoto_Matrix_V(k, 1) = Vrubusco_adj;
        ePhoto_Matrix_V(k, 2) = CA(i);
        ePhoto_Matrix_V(k, 3) = Ci;
        ePhoto_Matrix_V(k, 4) = GrossAssimilation;
        ePhoto_Matrix_V(k, 5) = Cc; 
        k = k + 1; % Move to the next row
    end
end

%Create vector to store Vmax_adj, differences between gross assimilation rates and corresponding SSRs
Diff_Matrix_V = zeros(250,2);

% Compute differences between two photosynthetic models
for k = 1:length(Farq_Matrix_V)
    for j = 1:50
    Diff_Matrix_V(k,1) = Farq_Matrix_V(k,1);
    Diff_Matrix_V(k,2) = (Farq_Matrix_V(k,4)-ePhoto_Matrix_V(k,4))^2;%the squares of the residuals
        while k<250
        k = k + 1; % Move to the next row
        end
    end
end

% Get sums of squared residuals by summing every five rows 
SSR_Matrix_V = zeros(50,2);
Vrubusco_adj_vals = Diff_Matrix_V(:,1);
SSR_Matrix_V(:,1) = Vrubusco_adj_vals(5:5:end);
SSR_Matrix_V(:,2) = sum(reshape(Diff_Matrix_V(:,2), 5, []))';
% Get full list of gross assimilation rates associated with each Vrubisco
% sampled and for each Cc tested
numbers_V = reshape(ePhoto_Matrix_V(:,4),[5,50]);
%toc

% Find scaling factor with lowest SSR
[~, scaling_index_V] = min(SSR_Matrix_V(:,2));
a_Rubisco = SSR_Matrix_V(scaling_index_V, 1);
min_SSR_V = SSR_Matrix_V(scaling_index_V, 2);

% Save results
save Vcmax_simple_new_result.mat;
save("WeatherTemp.mat","WeatherTemp");

% % Plot scaling factors (optional)
% fig = figure;
% scatter(SSR_Matrix_V(:,1),SSR_Matrix_V(:,2),'MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0]);
% hold on
% scatter(a_Rubisco,min_SSR_V,'MarkerEdgeColor','r','MarkerFaceColor','r')
% xticks(0.5:0.2:1.5);
% xlabel('α_{Rubisco}');
% ylabel('SSR');
% 
% % Set figure size
% set(fig, 'PaperUnits', 'inches');
% set(fig, 'PaperPosition', [0 0 6 4]);        % [left bottom width height]
% set(fig, 'PaperSize', [6 4]);                % Exact size of output file
% 
% % Export figure
% print(fig, 'Scaling_Factor_Optimisation_for_Rubisco_carboxylation.pdf', '-dpdf', '-r300');


