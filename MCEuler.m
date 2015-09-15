function [ r ] = MCEuler( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

k = .4;
theta = .05;
beta = .05;
r_0 = .03;

N = 20;
M = 20;
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

% computing the zcb interest approximations
interest_approx = ones(1,10);
for i = 1:10
    sum_int = exp(-delta * sum(r(:,1:i*(N/10)),2));
    interest_approx(i) = sum(sum_int)/M;
end









end

