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


function twopt = twoPoint(pop,popSize);
% This routine is not used in the current model. 

% Global variables
global VmaxNum;
global NConc;
global BK;
global MW;
t = 1;

VmaxNumF = VmaxNum+2;

dominantPop = zeros(VmaxNumF,popSize/4);
recessivePop = zeros(VmaxNumF,popSize/4);

offspring1 = zeros(VmaxNumF,popSize/4);
offspring2 = zeros(VmaxNumF,popSize/4);
offspring3 = zeros(VmaxNumF,popSize/4);
offspring4 = zeros(VmaxNumF,popSize/4);

pop = orderPop(pop,popSize);
tempPop = zeros(VmaxNumF,popSize/2);
finalPop = zeros(VmaxNumF,popSize);



%store top 50% in an array
for j = 1:popSize/2
    for c = 1:VmaxNumF
        tempPop(c,j)=pop(c,j);
    end
end

pop = tempPop;


%extract dominant mates
t = 1;
x = 1;
while t <=(popSize/2)
    for c = 1:VmaxNumF
        dominantPop(c,x)= pop(c,t);
    end         
    x = x + 1;
    t = t + 2;
end



%extract recessive mates
t = 2;
x = 1;
while t <=(popSize/2)
    for c = 1:VmaxNumF
        recessivePop(c,x)= pop(c,t);
    end         
    x = x + 1;
    t = t + 2;
end


%First Offspring of all the Mates

points = zeros(2,2);

%Generates random permutation of 12 numbers and stores into an array

%Dominant 

for y = 1:popSize/4
    num = randperm(12);
    points(1,1) = num(1) + 2;
    points(2,1) = num(2) + 2;
    points(1,2) = recessivePop(num(1) + 2,y);
    points(2,2) = recessivePop(num(2) + 2,y);
    offspring1(num(1) + 2,y) = points(1,2);
    offspring1(num(2) + 2,y) = points(2,2);
    for c = 1:VmaxNum
        if (c ~= (num(1) + 2)) && (c ~= (num(2) + 2))
            offspring1(c,y)=dominantPop(c,y);
        end
    end

    num = randperm(12);
    points(1,1) = num(1) + 2;
    points(2,1) = num(2) + 2;
    points(1,2) = recessivePop(points(1,1),y);
    points(2,2) = recessivePop(points(2,1),y);
    offspring2(points(1,1),y) = points(1,2);
    offspring2(points(2,1),y) = points(2,2);
    for c = 1:VmaxNum
        if (c ~= (num(1) + 2)) && (c ~= (num(2) + 2))
            offspring2(c,y)=dominantPop(c,y);
        end
    end
   
end




for y = 1:popSize/4
    num = randperm(12);
    points(1,1) = num(1) + 2;
    points(2,1) = num(2) + 2;
    points(1,2) = dominantPop(num(1) + 2,y);
    points(2,2) = dominantPop(num(2) + 2,y);
    offspring3(num(1) + 2,y) = points(1,2);
    offspring3(num(2) + 2,y) = points(2,2);
    
    for c = 1:VmaxNum
        if (c ~= (num(1) + 2)) && (c ~= (num(2) + 2))
            offspring3(c,y)=recessivePop(c,y);
        end
    end
    

    num = randperm(12);
    points(1,1) = num(1) + 2;
    points(2,1) = num(2) + 2;
    points(1,2) = dominantPop(points(1,1),y);
    points(2,2) = dominantPop(points(2,1),y);
    offspring4(points(1,1),y) = points(1,2);
    offspring4(points(2,1),y) = points(2,2);
    
    for c = 1:VmaxNum
        if (c ~= (num(1) + 2)) && (c ~= (num(2) + 2))
            offspring4(c,y)=recessivePop(c,y);
        end
    end
end



temp = 1;
check = zeros(2);
for y = 1: popSize/2
    if y == 1
        temp = 1;
    else
        if check(1) == 1 && check(2) == 1
            temp = temp + 1;
            check = zeros(2);
        end
    end
    for j = 1:VmaxNum 
        if j == 1 || j == 2
            finalPop(j,y) = 0;
        else
            if mod(y,2) ~=0
                finalPop(j,y) = offspring1(j,temp);
                check(1) = 1;
            else
                finalPop(j,y) = offspring2(j,temp);
                check(2) = 1;
            end
        end
    end

end


temp = 1;
check = zeros(2);
for y = 1: popSize/2
    if y == 1
        temp = 1;
    else
        if check(1) == 1 && check(2) == 1
            temp = temp + 1;
            check = zeros(2);
        end
    end
    
    for j = 1:VmaxNum 
        if j == 1 || j == 2
            finalPop(j,y+popSize/2) = 0;
        else
            if mod(y,2) ~=0
                finalPop(j,y+popSize/2) = offspring3(j,temp);
                check(1) = 1;
            else
                finalPop(j,y+popSize/2) = offspring4(j,temp);
                check(2) = 1;
            end
        end
    end

end

twopt = finalPop;
