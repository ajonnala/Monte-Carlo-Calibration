function [ r ] = computeR( k,beta,theta,r_0 ,M,N)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

delta = 1/N;

% computing the Es
epsilon = normrnd(0,1,M,N);

% computing the r_t
r = ones(M,N);
r(:,1) = r_0; %initializes all r_0's 
for n = 2:N %makes all negative values 0
    val = r(:,n-1) + k*(theta - r(:,n-1))*delta + beta*(sqrt(delta)*epsilon(:,n));
    valIndic = (val < 0);
    val(valIndic) = 0;
    r(:,n) = val;
end


end

