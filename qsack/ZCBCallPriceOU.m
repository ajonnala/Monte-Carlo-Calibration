function [ bond_price ] = ZCBCallPriceOU ( kappa, theta, beta, r_0, K, t, T )

RES = 100;
N = round(T*RES); % number of time steps
M = 100000; % number of MC paths
delta = T/N; % time step 

% epsilon = samples of normal distribution
% M x N array.
epsilon = normrnd(0, 1, M, N);

% compute the interest rate paths.
% each row describes the evolution of one rate path over time.
r = ones(M, N);
r(:, 1) = r_0;
for n = 2:N
    % interest rates evolve according to SDE dynamics
    r(:, n) = r(:, n-1) + kappa*(theta - r(:, n-1))*delta + beta*sqrt(delta)*epsilon(:, n);
    % if interest rate cannot be negative according to this SDE, then: 
    %nextStep(nextStep < 0) = 0;
    %r(:,n) = nextStep;
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

