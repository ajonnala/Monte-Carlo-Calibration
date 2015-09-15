function [ price ] = ManyZCBCall(beta )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

k = .4;
theta = .05;
r_0 = .03;

N = 100;
M = 100;
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


interest_approx = exp(-delta * sum(r(:,1:N),2));

price = ones(1,N);
for i = 1:N
    price(i) = interest_approx(i)*ZCBCallOption(r(i,N),beta);
end
 
price = mean(price);
end

