function [ interest_approx ] = MCEuler( k,beta,theta,r_0 ,M,N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

delta = 1/N;

r = computeR(k,beta,theta,r_0 ,M,N*10);
% computing the zcb interest approximations
interest_approx = ones(1,10);
for i = 1:10
    sum_int = exp(-delta * sum(r(:,1:i*(N)),2));
    interest_approx(i) = sum(sum_int)/M;
end









end

