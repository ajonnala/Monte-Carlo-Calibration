function [ ans ] = ZCBAnalyticOU( kappa, theta, beta, r_t, t, T )

ans = exp(-aAnalytic(kappa, theta, beta, r_t, T-t) - ... 
          bAnalytic(kappa, theta, beta, r_t, T-t) * r_t); 


end

function [ ans ] = aAnalytic( kappa, theta, beta, r_t, s )

first = theta - (beta^2)/(2*kappa^2);
second = s - bAnalytic(kappa, theta, beta, r_t, s);
third = (beta^2)/(4*kappa) * bAnalytic(kappa, theta, beta, r_t, s)^2;
ans = first * second + third;

end


function [ ans ] = bAnalytic( kappa, theta, beta, r_t, s )

numerator = 1 - exp(-s * kappa);
denominator = kappa;
ans = numerator/denominator;

end