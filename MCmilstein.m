function [ interest_approx ] = MCmilstein( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

k = .4;
theta = .05;
beta = .05;
r_0 = .03;

N = 1000;
M = 1500;
delta = 1/N;

% computing the Es
epsilon = normrnd(0,sqrt(delta),M,N); % should this be delta or 1?

% computing the r_t
r = ones(M,N);
r(:,1) = r_0; %initializes all r_0's 
for n = 2:N
    r(:,n) = r(:,n-1) + k*(theta - r(:,n-1))*delta + beta*(epsilon(n));
end

% computing the zcb interest approximations
interest_approx = ones(1,10);
for i = 1:10
    sum_int = exp(-delta * sum(r(:,1:i*(N/10)),2));
    interest_approx(i) = sum(sum_int)/M;
end








end


