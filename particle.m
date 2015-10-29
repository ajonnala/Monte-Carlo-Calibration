function [expectation] = particle(initialPrice, timeSteps, timeHorizon, ...
            numParticles, numStrikes, numEtaPaths )

% parameters
lambda = 1; % for the eta paths


deltaT = 1/(timeSteps) * timeHorizon; 


% matrix of eta paths
% dimensions: t x n'
% n' = num of paths, t' = timesteps
etaPaths = generateEtaPaths(timeSteps, deltaT, numEtaPaths, lambda);
etaSquared = etaPaths .^ 2;

% brownian paths
% dimensions: t x n
% n = num of paths, t = timesteps
w = generateBrownian(timeSteps, numParticles);

% particle step

expectation = zeros(timeSteps, numStrikes);
beta = zeros(timeSteps, numStrikes);
pricePaths = zeros(timeSteps, numParticles);

strikes = generateStrikes(timeSteps, numStrikes);
averageStrikes = genAverageStrikes(strikes);

% initialize the first row of pricePaths --> current spot
pricePaths(1, :) = initialPrice;
% initialize first row of expectation ---> 1
expectation(1, :) = 1;


for i=2:timeSteps
    % compute next price
    currentBeta = getBeta(pricePaths(i-1,:), strikeAverages(i-1), beta(i-1, :));
    pricePaths(i, :) = pricePaths(i-1, :) .* lambdaD(i, numParticles);
    pricePaths(i, :) = pricePaths(i, :) .* exp( sqrt(deltaT) * w(i-1, :) .* ...
            currentBeta .* etaPaths(i-1, :) - (1 / 2) * (currentBeta).^2 .* ... 
            etaSquared(i-1, :) * deltaT );
    
    % compute expectation
    
    deltaRow = delta(pricePaths(i), strikes(i), averageStrikes(i));
    numerator = etaSquared(i, strikeIdx) * deltaRow;
    denominator = sum(deltaRow);
    expectation(i, :) = numerator ./ denominator;
    
    % update the betas
    % be careful about K dependence for sigma^2
    % TODO: is this guaranteed to be >= 0?
    beta(i) = sqrt(getDupireVol(i, strike) ./ expectation(i, :));

    
    
end

end


% kernel wrapper
% returns vector
function delta(prices, strikes, avgStrikes)

output = zeros(1, length(avgStrikes));

for i = 1:length(avgStrikes)
    h = strikes(i+1) - strikes(i);
    %TODO FIX THIS
    output(i) = 1/h * kernel((prices(i) - avgStrikes(i))/h);
end

return output;

end


function kernel(x)

if x>(1/2) && x<=(-1/2)
    return 0;
end

return 1;

function genAverageStrikes(strikes)
% value at (T_i, S_j) = the strike for which beta(T_i, S_j) is computed
timeSteps = size(strikes, 1);
numAvgStrikes = size(strikes, 2) - 1; %-1 because we're only using avg points
avgStrikes = zeros(timeSteps, numAvgStrikes);
for row = 1:timeSteps
    for col = 1:numAvgStrikes
        avgStrikes(row, col) = (1/2)*(strikes(row, col) + strikes(row, col+1));
    end
end

return avgStrikes;

end

function getDupireVol()
    
end

function getBeta(prices, strikeAverages, betaRow)

% TODO: fix extrapolation. Right now, it extrapolates with the same
%    method as interpolation.
betaSq = interp1(strikeAverages, betaRow, prices, 'linear', 'extrap');
% TODO:remove this
%betaSq = zeros(size(betaRow));
%for i=1:length(prices);   
%     leftIdx = search(prices(i));
%     betaSq(i) = interpolate(betaRow, leftIdx, price(i));
end

return betaSq;

end

% TODO: remove this if not necessary
function interpolate(betaRow, leftInx)



end

function search(x)

end

% drift factor
%TODO: what is this supposed to be?
function lambdaD(time, n)
% return vector
return ones(1, n);

end

function generateStrikes(t, n)

% for now, uniformly generate strikes between 0 and M
M = 10;

strikes = zeros(1, n);
for i=1:n
    strikes(i) = M/(n-1)*(i-1);
end
% uniform across time for now
return repmat(strikes, t, 1);

end


function generateBrownian(t, N)

return normrnd(0, 1, t, N);

end


function generateEtaPaths(timeSteps, deltaT, numEtaPaths, lambda)

% dimension: time x paths
W = generateBrownian(timeSteps, numEtaPaths);

% forward euler
Z = zeros(timeSteps, numEtaPaths);
etaPaths = zeros(timeSteps, numEtaPaths);

for n = 2:timeSteps
    Z(:, n) = Z(:, n-1) -lambda*r(:, n-1)*deltaT + ... 
        gammaM((n-1)*deltaT)*sqrt(deltaT)*W(:, n);
end

% compute variances. variance(i) = variance at time i
variance = zeros(timeSteps, 1);
V = zeros(timeSteps, 1);
for i = 2:timeSteps
    % left endpoint approximation for the integral
    V(i) = V(i-1) + (gammaM((i-1)*deltaT))^2 * ... 
            exp(2*lambda) * deltaT;
    variance(i) = exp(-2*lambda*i*deltaT) * V(i);
end

% compute the eta paths finally
etaPaths = exp(Z - repmat(variance, 1, numEtaPaths));

return etaPaths;

end



function gammaM(t)

return 1;

end