function [ ZCB] = ZCBCallOption( r_0,beta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

k = .4;
theta = .05;
K = .9;

b = inline('(1/k)*(1 - exp(-s*k))','k','s');
a = inline('(t - (beta^2/(2*k^2))) * (s - bs) + (beta^2/(4*k))*(bs^2)','t','beta','k','bs','s');

val = exp(-a(theta,beta,k,b(k,1),1) - b(k,1)*r_0) - K ;
if (val < 0) 
    ZCB = 0;
else
    ZCB = val;
       
    




end

