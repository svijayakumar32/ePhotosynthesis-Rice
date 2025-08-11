%% Plot A/Cc curve for measured average parameters from IRRI data and using Ac and Aj equations from Vcmax_adj_simple and Jmax_adj_simple
% Vcmax, J and TPU parameters are the averaged values of parameters derived from fitting eight curves 
% using the R package msuRACiFit described in Gregory et al (2021) and is available at: https://github.com/poales/msuRACiFit
%% Farquhar model parameters
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
TPU = mean(parameters.TPU);
Gr = mean(Gstar_vals.x);
Rd = mean(parameters.rL);
Gm = mean(Gm_vals{:,:}); 
Kc_air = mean(Kcair_vals{:,:}); 

Lii=2000;% Light intensity from IRRI
O=210;%mbar

%%%%%%%%%%%%%%%%%%%%%

%% e-Photosynthesis model parameters
%% For running the untuned e-Photosynthesis model
%global Vrubusco_nonopt;
Vrubusco_nonopt=1.0;%adjust enzyme activity to 1
%global VmaxAdj_nonopt;
VmaxAdj_nonopt=1.0;%adjust enzyme activity to 1 

%% For running the tuned e-Photosynthesis model
%global Vrubusco_opt;
Vrubusco_opt = 1.36;%adjust enzyme activity sub in optimal Vrubisco 
%global Vmax_opt;
VmaxAdj_opt = 1.12;%adjust enzyme activity i.e. sub in optimal Vmaxadj

global pcfactor;  
ProteinTotalRatio=0;
%pcfactor=1/ProteinTotalRatio;
%21/05/24 Change pcfactor to 1 here to avoid Inf
pcfactor=1;

Einput=ones(37,1);%No gene expression data input
Edata=importdata('Einput7.txt');
%Eio=Edata.data(:,1);
Eio_nonopt=Edata.data(:,1);
Eio_opt=Edata.data(:,1);

%MetaOnly=1;% if MetaOnly=1 run only Metabolic model
WeatherTemp = mean(temp_vals.x); %Avg Tleaf, Original = 25C
GRNC=0; % In EPS_Drive_GRNs, if GRNC==0 cATPsyn,CPSi,cNADPHsyn and cpsii=1 are all set to 1

%%%%%%%%%%%%%%%%%%%%%

%% Specify range of CAs 
%Ca=(20:20:1200)'; % 60 CA values
Ca_Gr = 67.4;
Ca = (70:20:1200)'; % 60 CA values
%Ca = (70:10:1200)'; % 120 CA values
Ca = vertcat(Ca_Gr,Ca);

%Ca=[0:1:2000]';

%% Preallocate vectors to store outputs including Cc calculated from CA, and
%% constants for Rubisco limited assimilation, RuBP regeneration limited
%% assimilation and TPU limited assimilation 
%% and also e-Photosynthesis model assimilation 
%% using initial (a_Enzymes/Rubisco=1) and adjusted (a_Enzymes=1.12, a_Rubisco=1.32) scaling factors

[Ci, ...
    Ac_b, Ac_c, Aj_b, ...
    Aj_c, Ap_b, Ap_c, ...
    Net_Ac, Net_Aj, Net_Ap, ...
    Cc_Ac, Cc_Aj, Cc_Ap, Cc_ePhoto, ...
    Gross_Ac, Gross_Aj, Gross_Ap, Gross_A_nonopt, Gross_A_opt] = deal(zeros(numel(Ca),1));

%ePhoto_Matrix = zeros(numel(Ca),6);

for i = 1:numel(Ca)
% Calculate Ci from CA and various assimilation rates for Farquhar model
    Ci(i)= 0.7*Ca(i);
    Ac_b(i) = Vcmax_m-Rd+(Ci(i)+Kc_air) * Gm;
    Ac_c(i) = ((Ci(i)-Gr) * Vcmax_m-(Ci(i)+Kc_air) * Rd) * Gm;
    Aj_b(i) = (J/4) - Rd + (Ci(i)+ 2 * Gr) * Gm;
    Aj_c(i) = ((Ci(i) - Gr) * J/4 - (Ci(i) + 2 * Gr) * Rd) * Gm;
% Calculate net assimilation rates under Rubisco limitation
    Net_Ac(i) = ((Ac_b(i) - sqrt(Ac_b(i)^2 - 4 * Ac_c(i))) / 2);
    Cc_Ac(i) = (Ci(i) - Net_Ac(i)/Gm); % Use Net_Ac to calculate Cc
    Gross_Ac(i) = Net_Ac(i) + Rd; 
% Calculate net assimilation rates under RuBP limitation
    Net_Aj(i) = ((Aj_b(i) - sqrt(Aj_b(i)^2 - 4 * Aj_c(i))) / 2);
    Cc_Aj(i) = (Ci(i) - Net_Aj(i)/Gm); % Use Net_Aj to calculate Cc
    Gross_Aj(i) = Net_Aj(i) + Rd; 
% Calculate net assimilation rates under TPU limitation
    Net_Ap(i) = 3 * TPU-Rd;
    Cc_Ap(i) = (Ci(i) - Net_Ap(i)/Gm); % Use Net_Ap to calculate Cc
    Gross_Ap(i) = Net_Ap(i) + Rd;
% Adjust enzyme activities for e_Photosynthesis - nonopt is both alpha
% values = 1, opt is aRubisco = 1.32, aEnzymes = 1.12
    Eio_nonopt(1) = Edata.data(1,1) * Vrubusco_nonopt;      %aRubisco = 1
    Eio_nonopt(2:27) = Edata.data(2:27,1) * VmaxAdj_nonopt; %aEnzymes = 1
    Eio_opt(1) = Edata.data(1,1) * Vrubusco_opt;            %aRubisco = 1.32
    Eio_opt(2:27) = Edata.data(2:27,1) * VmaxAdj_opt;       %aEnzymes = 1.12  
    % e-Photosynthesis model A calculation 
    Air_CO2 = Ca(i);
    Ci = Air_CO2 * 0.7; % intercellular CO2
    % Choose which A is most limited to calculate Cc
    if     Net_Ac(i) < Net_Aj(i)
           Cc_ePhoto(i) = Cc_Ac(i);
    elseif Net_Aj(i) < Net_Ap(i)
           Cc_ePhoto(i) = Cc_Aj(i);
    else   
           Cc_ePhoto(i) = Cc_Ap(i);
    end
%  Calculate Gross Assimilation for e-Photosynthesis model
    PPFDi = Lii; 
    Gross_A_nonopt(i) = EPS_Drive_GRNs(Einput,Cc_ePhoto(i),PPFDi,WeatherTemp,GRNC,0,Eio_nonopt); %aRubisco = 1, aEnzymes = 1
    Gross_A_opt(i) = EPS_Drive_GRNs(Einput,Cc_ePhoto(i),PPFDi,WeatherTemp,GRNC,0,Eio_opt); %aRubisco = 1.32, aEnzymes = 1.12
end

save('Farq_ePhoto_results_new.mat');
%%%%%%%%%%%%%%%%%%%%%
