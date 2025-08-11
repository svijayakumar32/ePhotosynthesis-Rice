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


% This program ranked the populations of enzymes activity based on the CO2
% uptake value. 

function rankPopulation = rankPop(pop, popSize)

    global VmaxNum;
    
    %%%%%%%% Evaluate Fitness %%%%%%%%%%
%     NewM = zeros(VmaxNum+2,popSize);
    newM = zeros(VmaxNum+2,popSize);

    CO2Array = zeros(popSize,1);

% % Faster rewrite - check it works    
% % Sort the CO2Array in descending order
%     [out,idx] = sort(pop(2,:),'descend');
    
    for j = 1:popSize
        CO2Array(j) = pop(2,j);
    end
    
    for m = 1:popSize
        % Find the index of the highest CO2 concentration. 
        [x,maxindex] = max(CO2Array);
        CO2Array(maxindex) = -1;       
        NewM(1,m) = m;
        NewM(2,m) = x;
        NewM(3:VmaxNum+2,m)=pop(3:VmaxNum+2,maxindex);
    end
    
    rankPopulation = NewM;

    % % Faster rewrite - check it works
    % rankPopulation = pop(:,idx);
    % rankPopulation(1,:) = 1:popSize;
    % rankPopulation(2,:) = out;
