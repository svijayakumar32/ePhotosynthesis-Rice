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

function resizePopulation = resizePop(pop, popSize)

for m = 1:3
    temp0 = popSize/4;
    temp1 = popSize/4*m + 1;
    temp2 = popSize/4*(m+1);
    pop(:,temp1:temp2)= pop(:,1:temp0);
end

resizePopulation = pop;

% % Faster rewrite of entire function - check it works 
% function pop = resizePop(pop, popSize)
% temp0 = popSize/4;
% for m = 1:3
%     temp1 = popSize/4*m + 1;
%     temp2 = popSize/4*(m+1);
%     pop(:,temp1:temp2)= pop(:,1:temp0);
% end
