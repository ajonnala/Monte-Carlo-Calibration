function [ min_k,min_beta,min_error ] = kBetaMinimization( theta,M,N )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


r_0 = .05;
min_error = Inf;
min_k = -1;
min_beta = -1;

x = ones(1,100);
y = ones(1,100);
z = ones(1,100);
i = 1;
for k = 0:.01:1
    for beta = 0:.01:1
        approx = MCEuler(k,beta,theta,r_0,M,N);
        derived = ZCBderived(k,beta,theta,r_0);
    
    
        error = sum((approx - derived).^2);
        x(1,i) = k;
        y(1,i) = beta;
        z(1,i) = error;
        
        if (error <= min_error)
            min_k = k;
            min_beta = beta;
            min_error = error;
        end
    end
end

figure
plot3(x,y,z)

end

