function [ interest_approx ] = EulerFeller( k,beta,theta,r_0,M,N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


delta = 1/N;

% computing the Es
epsilon = normrnd(0,1,M,N);

% computing the r_t
r = ones(M,N);
r(:,1) = r_0; %initializes all r_0's 
for n = 2:N %makes all negative values 0
    val = r(:,n-1) + k*(theta - r(:,n-1))*delta + beta*(sqrt(r(:,n-1)).*sqrt(delta).*epsilon(:,n));
    valIndic = (val < 0);
    val(valIndic) = 0;
    r(:,n) = val;
end



% computing the zcb interest approximations
interest_approx = ones(1,1);
for i = 1:1
    sum_int = exp(-delta * sum(r(:,1:(N)),2));
    interest_approx(i) = sum(sum_int)/M;
end




end

