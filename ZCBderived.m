function [ ZCB ] = ZCBderived( k,beta,theta,r_0 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



b = inline('(1/k)*(1 - exp(-s*k))','k','s');
a = inline('(t - (beta^2/(2*k^2))) * (s - bs) + (beta^2/(4*k))*(bs^2)','t','beta','k','bs','s');

ZCB = ones(1,10);
for i = 1:10
    ZCB(i) = exp(-a(theta,beta,k,b(k,i),i) - b(k,i)*r_0);
end




end

