%   Copyright   Xin-Guang Zhu and Stephen P. Long, University of Illinois 
%   Copyright ©  2007

%   This file is part of CarbonMetabolism.

%    CarbonMetabolism is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.

%    CarbonMetabolism is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License (GPL)
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function gpmain_rice_all(Ci)
% Get absolute path of this script
ScriptPath = mfilename('fullpath');

% Locate directory containing this script
ScriptDir = fileparts(ScriptPath);

% Change directory to main e-Photosynthesis code repository
cd(fullfile(ScriptDir,'..','..'));
ePhotosynthesis_repository = pwd;

% Add repository to path
addpath(genpath(ePhotosynthesis_repository));

% Load temperature file
load("WeatherTemp.mat"); %Avg Tleaf, Original = 25C

% Adjust alpha values here and read in Einput7.txt
global Vrubusco_adj;
Vrubusco_adj = 1.36; 
global VmaxAdj;
VmaxAdj = 1.12;

%Ci=__;% Adjust Ci= 140-420 here in script or
CO2i=160;% Ci given as input to function
PPFDi=2000;
GRNC=0;
Einput=ones(37,1);%No gene expression data input
Edata=importdata('Einput7.txt');
Eio=Edata.data(:,1);
Eio(1)=Edata.data(1,1)*Vrubusco_adj;
Eio(2:26)=Edata.data(2:26,1)*VmaxAdj;
Enzyme=importdata('MW&Kcat.txt');
MWKcat=Enzyme.data;
MWKcat([7,9,12],:) = []; %remove rows corresponding to V8, V10 and V16
%%%%%%%%%%%%%%%%%%% Initialize Variables %%%%%%%%%%%%%%%%%%%
global pcfactor; 
pcfactor=1;
global VmaxNum;
%global NConc; 
% IntContinue = 0;        % This is a flag to indicate whether we need to use the external data or not for the initialization. Default is 0.
% WriteExcel = 0;

% The following variables are used to output data from the simulation. 

global d_plot;
global Tt_plot;

% global NetCarbon_plot;
% global PS_VEL_gpmain;

global GP;
GP = 1;         % This is a parameter used for transfer infomration to the PSRate file. 
global TimeBegin;
TimeBegin = clock;

global Optiona;
Optiona =0;     % The option used in the mutate.m file for creating mutation.

%WARNING population size (PopSize) must be divisible by 4
popSize = 16;   %default=16    
numofGen = 1500;       % The total number of generations default = 1500 %test i=100  
mutatePercentage = 0.02;     % Maximal percentage of changes in Vmax in each generation

% Different options of generating new population. Default is 1. 
generationTransfer = 1; 
VmaxNum = 23;           % The number of enzymes used in the optimization. %originally 26 

factor = 1;
ScaleR = factor;        % This is a factor used to modify the enzyme concentration.

%%%% Initialize Population Array %%%%

% Coeff = 1/30/3600*1000;

pop = zeros(VmaxNum+2,popSize);

% In the pop matrix, the first element is rank; 2rd: the CO2 uptake rate; 3:VmaxNum + 2, different Vmax;

for i = 1:popSize
pop(3:25,i)=Eio([1:6,8,10,11,13:26]); %taking out V8, V10 and V16 in optimization
end


%%%% Initialize BK and MW Array %%%%

global BK;
BK=MWKcat(:,2);
%BK([7,9,12],:) = [];
global MW;
MW=MWKcat(:,3);
%MW([7,9,12],:) = [];


% Calculate the default nitrogen concentration
sumd = 0;
for k = 1:VmaxNum
    if k == 1
    sumd= sumd + (1.1*pop(k+2,1))/BK(k)*MW(k); % For Rubisco, adjust activity to 80% - 10% inhibition
    else
    sumd= sumd + pop(k+2,1)/BK(k)*MW(k);
    end
end
%sumd = sumd- pop(9,1)/BK(7)*MW(7)-pop(11,1)/BK(9)*MW(9)-pop(13,1)/BK(11)*MW(11)-pop(14,1)/BK(12)*MW(12);      % The transketolase is double counted. So, it is corrected at here. 


% NConcDefault = sumd/1000/33 * 0.16;        % This is the default nitrogen concentration for those enzymes in the carbon metabolism pathway. 
% NConc = NConcDefault * factor;
% NConc = NConc * 1000 * 33 /0.16       % mg protein l-1

global NTotal;
%NTotal = sumd * factor
NTotal = sumd
%%%% Mutate Initial Population %%%%
tempn = popSize/8;
temp_pop=pop(:,1:tempn);
pop = mutate(pop, popSize, mutatePercentage);
pop(:,1:tempn) = temp_pop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global gp2condition_RuCon;
CO2PLOT = zeros(3,5);  
% 1: Generation number; 2: Average; 3, Error bar; 

TAVMatrix = zeros(numofGen,VmaxNum+2);
TSTDMatrix = zeros(numofGen,VmaxNum+2);
TBestMatrix = zeros(numofGen,VmaxNum+2);

for i = 1:numofGen
    
    %%%%%%%%%%%%% Test Population %%%%%%%%%%%%%%
    
% Store CO2 Assimilation Rate into pop
    
      disp('generation:');i
%     display('Time passed:'); 
%     Lapse = clock-TimeBegin

% Load constraints for photorespiratory enzymes (= average of values calculated during enzyme adjustment running at 129 ppm)
PR_constraints = importdata('PR_constraints.txt');    
%PR_constraints=PR_constraints_data.data(:,1);

    for j = 1:popSize
%         gp2condition_RuCon = pop(3,j)/2;
%         global gp2V111;
%         gp2V111 = pop(3,j) * 0.24;            % Here the rate of Rubisco oxygenation is adjusted.
            %Temp = CM_Drive(pop, j);
            Eiopop=Eio([1:6,8,10,11,13:66]); %remove three enzymes instead of Eiopop=Eio
            Eiopop(1:23)=pop(3:25,j); %adjust to 23 enzymes
%           Eiopop(7)=Eiopop(4); %These duplicates are already missing from the optimization list so no need to set them as equal
%           Eiopop(9)=Eiopop(6);
% 	        Eiopop(11)=Eio(11);
%           Eiopop(12)=Eio(12);%VATP synthesis not used

            %Fixing minimal enzyme level constraints for enzymes involved in starch and sucrose/photorespiratory metabolism
            %from adjusted enzyme levels in Einput_rice but missing out V8, V10 and V16 from optimization
            if Eiopop(9)<Eio(11) %V23
            disp('Glucose-1-phosphate adenylyltransferase value is below the limit.')
            continue
            elseif Eiopop(17)<Eio(20) %V51
            disp('Fructose-bisphosphate aldolase (C) value is below the limit.')
            continue
            elseif Eiopop(18)<Eio(21) %V52
            disp('Fructose-bisphosphatase (C) value is below the limit.')
            continue
            elseif Eiopop(19)<Eio(22) %V55
            disp('UTP-glucose-1-phosphate uridylyltransferase value is below the limit.')
            continue
            elseif Eiopop(20)<Eio(23) %V56
            disp('Sucrose-phosphate synthase value is below the limit.')
            continue
            elseif Eiopop(21)<Eio(24) %V57
            disp('Sucrose-phosphate phosphatase value is below the limit.')
            continue
            elseif Eiopop(22)<Eio(25) %V58
            disp('Fructose-2,6-bisphosphate 2-phosphatase value is below the limit.')
            continue
            elseif Eiopop(23)<Eio(26) %V59
            disp('6-phosphofructo-2-kinase value is below the limit.')
            continue
			%Fixing minimal enzyme level constraints for photorespiratory enzymes
            elseif Eiopop(10)< PR_constraints(1) %V112  
            disp('Phosphoglycolate phosphatase value is below the limit.')
            continue
            elseif Eiopop(11)< PR_constraints(2) %V113
            disp('Glycerate kinase value is below the limit.')
            continue
            elseif Eiopop(12)< PR_constraints(3) %V121
            disp('(S)-2-hydroxy-acid oxidase & Catalase value is below the limit.')
            continue
            elseif Eiopop(13)< PR_constraints(4) %V122
            disp('Serine-glyoxylate transaminase value is below the limit.')
            continue
            elseif Eiopop(14)< PR_constraints(5) %V123
            disp('Glycerate dehydrogenase value is below the limit.')
            continue
            elseif Eiopop(15)< PR_constraints(6) %V124
            disp('Glycine transaminase value is below the limit.')
            continue
            elseif Eiopop(16)< PR_constraints(7) %V131
            disp('Glycine dehydrogenase (aminomethyl-transferring) value is below the limit.')
            continue
            else
            end % All values within constraints, so analysis can proceed
            %Use Eiopop from optimization loop, setting V7=V4 and V10=V7 
%             Eiopop=vertcat(Eiopop(1:6),Eiopop(4),Eiopop(7),Eiopop(6),Eiopop(8:63)); 
%             Temp=EPS_Drive_GRNs(Einput,CO2i,PPFDi,WeatherTemp,GRNC,0,Eiopop);
            % Set ATP synthase as equal to its starting value in Eio
            Eiopop_phot=vertcat(Eiopop(1:6),Eiopop(4),Eiopop(7),Eiopop(6),Eiopop(8:9),Eio(12),Eiopop(10:63)); 
            Temp=EPS_Drive_GRNs(Einput,CO2i,PPFDi,WeatherTemp,GRNC,0,Eiopop_phot);
            %Check whether the concentrations of metabolites reach steady states#
            sizeT=size(d_plot);
            global tglobal;
            Tcheck = tglobal * 4/5;             
            
            Tindex = find(Tt_plot>Tcheck);    
            IndexMiddle = Tindex(1);
            ddiff=d_plot(sizeT(1),53:87)-d_plot(IndexMiddle,53:87);%metabolites
            dnorm1 = norm(ddiff);
            tdiff = Tt_plot(sizeT(1))-Tt_plot(IndexMiddle);
            sloped = dnorm1/tdiff/0.01;
            
            if sloped < 10^(-3)
                pop(2,j) = Temp;
            else
                pop(2,j) = 0;
            end
            
            
            test = pop(2,j);
   
    end

    %%%%%%%%%%% Rank Population %%%%%%%%%%%
    pop = rankPop(pop,popSize);
    
    %%%%%% Save optimization Statistics %%%%%%%%
    AVMatrix = average(pop,popSize);
    STDMatrix = stdev(pop,popSize);
    BestMatrix = pop(:,1)';
    
    TAVMatrix(i,:)=AVMatrix;
    TSTDMatrix(i,:)=STDMatrix;
    TBestMatrix(i,:)=BestMatrix;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Check if this is the last generation
    CO2PLOT(1,i) = i;
    CO2PLOT(2,i) = AVMatrix(2);
    CO2PLOT(3,i) = STDMatrix(2);

    %mcc nojvm does not use the Java Virtual Machine (JVM) so comment out
    %all figure rendering commands to avoid errors
    %errorbar(CO2PLOT(1,:),CO2PLOT(2,:),CO2PLOT(3,:)); pause(5);
    %xlim([1, i+1]); ylim([0,80]);  
        %Resize Pop Array
        switch generationTransfer
            case 1
                pop = resizePop(pop, popSize);
            case 2
                pop = twoPoint(pop,popSize);
        end
        % Mutate Population
        pop = mutate(pop, popSize, mutatePercentage);
	    if rem(i/100,1)==0
        i
        Ci_str = num2str(CO2i);
        task_id = getenv('SLURM_ARRAY_TASK_ID');
            if isempty(task_id) % if running on local machine and there is no task_id
                task_id = 'local';  % set as 'local'
            end
        % Specify unique filenames for each Ci
        workspacefileName = strcat ("CO2_rice_",Ci_str,"_",task_id,".mat");
        BestMatrixfileName = strcat ("outputenz_",Ci_str,"_",task_id,".txt");
        d_plotfileName = strcat ("d_plot_",Ci_str,"_",task_id,".csv");
   
        % Create directories to save results
        workspaceDir = fullfile(ePhotosynthesis_repository,"Optimisation", "Results", "Workspaces");
        if ~exist(workspaceDir, 'dir')
            mkdir(workspaceDir);
        end

        enzymeDir = fullfile(ePhotosynthesis_repository,"Optimisation", "Results", "Enzymes");
        if ~exist(enzymeDir, 'dir')
            mkdir(enzymeDir);
        end

        metaboliteDir = fullfile(ePhotosynthesis_repository,"Optimisation", "Results", "Metabolites");
        if ~exist(metaboliteDir, 'dir')
            mkdir(metaboliteDir);
        end

        % Debug output
        % disp("Saving workspace to:");
        % disp(fullfile(workspaceDir, workspacefileName));
        % 
        % disp("Saving BestMatrix to:");
        % disp(fullfile(enzymeDir, BestMatrixfileName));
        % 
        % disp("Saving metabolites to:");
        % disp(fullfile(metaboliteDir, d_plotfileName));

	    % Save the workspace 
        save(fullfile(workspaceDir, workspacefileName));
        % Save matrix of optimal enzyme rates to output file, transposing BestMatrix
        writematrix(BestMatrix', fullfile("Optimisation","Results","Enzymes", BestMatrixfileName));
        % Save .csv of metabolite concentrations
        writematrix(d_plot, fullfile("Optimisation","Results","Metabolites", d_plotfileName));
      end
end