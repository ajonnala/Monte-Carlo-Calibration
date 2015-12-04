function [] = graphExpStrikes(tensor, col, t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

a = tensor(1,col);
b = tensor(2,col);
exp = a{1,1};
strikes = b{1,1};

row = size(exp, 1);
%plotStyle = {'b-','k:','r.'};
colors = jet(row); 
hold on;
for i = 1:row
    plot(strikes(i,:), exp(i,:), 'color', colors(i,:)); %plotStyle{mod(i,3) +1}
    legendInfo{i} = ['TimeSteps = ' num2str(i)]; 
end
legend(legendInfo)
title(t);
end
