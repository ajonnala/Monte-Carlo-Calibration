function [ ans ] = ZCBMilsteinFeller( kappa, theta, beta, r_0, T )


% parameters
%kappa = .4;
%theta = .05;
%beta = .05;
%r_0 = .03;

N = 100; % number of time steps
M = 1000000; % number of MC paths
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
    nextStep = r(:, n-1) + kappa * (theta - r(:, n-1)) * delta + ...
               beta * sqrt(delta * r(:, n-1)) .* epsilon(:, n) + ...
               (1/4) * beta^2 * delta * (epsilon(:, n).^2 - 1);
    % if interest rate cannot be negative according to this SDE, then: 
    nextStep(nextStep < 0) = 0;
    r(:,n) = nextStep;
end

% approximate the integral
integrals = exp(-delta * sum(r, 2));
ans = mean(integrals);

