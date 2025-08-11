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
%    along with this program.  If not, see <https://urldefense.com/v3/__http://www.gnu.org/licenses/__;!!DZ3fjg!8kTVQ33-Ja6a2Eqd8SR_oi8FCPu68sOy6wRp53rq8uhfq71Mwo9-pQsTCZRJn-K9u0eJ2EYbTF4pO1agj3AfwuXD_lPEPqeDvPA$ >.



function returnPop = mutate(pop, popSize, mutatePercentage)

%%%% Global Variables %%%%

global BK;
global MW;
global VmaxNum;

global NTotal;

tempPop = zeros((VmaxNum+2),1);   % Create a vector of the size of Vmax+2
for j = 2:popSize
                sum=0;
                for m = 3:(VmaxNum+2)
			        rng shuffle % shuffle random values so there is no repetition of sequence from the random number generator
                    %This loop sets a non-negativity constraint when generating randval using absolute values for V23,V51,V52,V55,V56,V57,V58,V59 (sucrose/starch metabolism enzymes)
                    if (m == 11 ||m >= 19) 
                        randval = abs(1 - 2 * rand(1))*0.01; % Change mutatePercentage for sucrose/starch enzymes from 0.02 to 0.0002 through multiplying by 0.01
                    else
                        randval = 1 - 2 * rand(1); 
                    end

					tempPop(m) = (1+randval * mutatePercentage) * pop(m,j);

					if m == 3 
						sum = sum + (1.1*tempPop(m)/BK(m-2))*MW(m-2);  %For Rubisco, adjust activity by 80% - 10% inhibition
					else
                        sum = sum + (tempPop(m)/BK(m-2))*MW(m-2);       % mg protein l-1
					end
                end

                Ratio = NTotal/sum;
                pop(:,j) = tempPop * Ratio; 
end

%%%
% % Faster rewrite w/o loops - check it works, especially with Rubisco readjustment - saved as mutate_mod
% % Generate random values for the population size
% randvals = rand((VmaxNum+2),popSize); 
% randvals=1- randvals*2; 
% % Change mutatePercentage for sucrose/starch enzymes from 0.02 to 0.0002 through multiplying by 0.01
% randvals(11) = abs(randvals(11)) * 0.01; 
% randvals(19:VmaxNum+2) = abs(randvals(19:VmaxNum+2))*0.01; 
% % Multiply the random values by the mutatePercentage
% randvals(:,:) = (randvals(:,:)+1) * mutatePercentage; 
% 
% temp_pop = randvals(3:(VmaxNum+2),:) .* pop(3:(VmaxNum+2),:);
% 
% items = temp_pop(:,:);
% items(1,:) = items(1,:)*1.1; % previously 1.2
% Ratio = repmat(NTotal,[1,popSize])./sum(items./repmat(BK(1:(VmaxNum)),[1,popSize]).* repmat(MW(1:(VmaxNum)),[1,popSize]),1);
% 
% pop(3:(VmaxNum+2),2:popSize) = temp_pop(:,2:popSize) .* Ratio(2:popSize); 
%%%

returnPop = pop;
