function [ min_theta,min_error ] = thetaMinimization( k,beta,M,N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

r_0 = .05;

min_error = Inf;
min_theta = -1;

x = ones(1,100);
y = ones(1,100);
i = 1;
for theta = 0:.01:1 %minimization
    approx = MCEuler(k,beta,theta,r_0,M,N);
    derived = ZCBderived(k,beta,theta,r_0);
    
    error = sum((approx - derived).^2);
    x(1,i) = theta;  
    y(1,i) = error;
    i=i+1;
    if (error <= min_error)
        min_theta = theta;
        min_error = error;
    end
    
end

plot(x,y)

end

