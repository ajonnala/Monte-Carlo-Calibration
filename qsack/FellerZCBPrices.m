
function [prices] = FellerZCBPrices()
prices=ones(22,11);
for i=1:10
  prices(1, i+1) = 0.1*i;
end

for i=1:21
  prices(i+1, 1) = 0.89 + 0.01*i;
end


for T=2:11
  T
  for K=2:22
    prices(K, T) = ZCBCallPriceFeller(.4, .05, .05, 0.03, prices(K, 1), prices(1, T), prices(1, T) + 1);
  end
end

end


function [ bond_price ] = ZCBCallPriceFeller ( kappa, theta, beta, r_0, K, t, T )

RES = 100;
N = round(T*RES); % number of time steps
M = 100000; % number of MC paths
delta = T/N; % time step 

% epsilon = samples of normal distribution
% M x NT array.
epsilon = normrnd(0, 1, M, N);

% compute the interest rate paths.
% each row describes the evolution of one rate path over time.
r = ones(M, N);
r(:, 1) = r_0;
for n = 2:N
    % interest rates evolve according to SDE dynamics
    nextStep = r(:, n-1) + kappa * (theta - r(:, n-1)) * delta + ...
               beta * sqrt(delta * r(:, n-1)) .* epsilon(:, n) + ...
               (1/4) * beta^2 * delta * (epsilon(:, n).^2 - 1);
    % if interest rate cannot be negative according to this SDE, then: 
    nextStep(nextStep < 0) = 0;
    r(:,n) = nextStep;
end

% integral of rate from t to T to find bond price at time t
r_bond = r(:, t*RES+1:end);
bond_price_t = exp(-delta * sum(r_bond, 2));
option_price_t = (bond_price_t - K);
option_price_t(option_price_t < 0) = 0;
%discount the bond price
r_discount = r(:, 1:t*RES);
discount = exp(-delta * sum(r_discount, 2));
bond_price = mean(option_price_t .* discount);






end


