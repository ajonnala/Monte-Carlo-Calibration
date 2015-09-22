function [ k,beta,theta ] = Minimization( k, beta, M ,N )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

err = .05;
[theta,e] = thetaMinimization(k,beta,M,N);
dif = Inf;
count = 0;
while dif < err ;
    [k,beta,e] = kBetaMinimization(theta,20,20);
    [thetaNew,e] = thetaMinimization(k,beta,M,N);
    dif = abs(theta - thetaNew);
    count = count + 1
end   

end

