function [expectation, averageStrikes] = particle(initialPrice, timeSteps, timeHorizon, ...
            numParticles, numStrikes, numEtaPaths )




% parameters
lambda = 1; % for the eta paths

% constant for the equidistant bins function. We can eventually get rid of
%   this. 
M = 10;

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

% numStrikes-1 because we're evaluating at midpoints
expectation = zeros(timeSteps, numStrikes);
beta = zeros(timeSteps, numStrikes);
pricePaths = zeros(timeSteps, numParticles);

% initialize strikes
bins = zeros(timeSteps, numStrikes+1);
averageStrikes = zeros(timeSteps, numStrikes);
bins(1, :) = equidistantBins(numStrikes+1, M);
averageStrikes(1, :) = genAverageBin(bins(1, :));


% initialize the first row of pricePaths --> current spot
pricePaths(1, :) = initialPrice;
% initialize first row of expectation ---> 1
expectation(1, :) = 1;


% initialize initial beta
beta(1, :) = getDupireVol(0, averageStrikes(1, :));

for i=2:timeSteps
    % compute next price
    currentBeta = getBeta(pricePaths(i-1,:), averageStrikes(i-1, :), beta(i-1, :));
    pricePaths(i, :) = pricePaths(i-1, :) .* lambdaD(i, numParticles);
    pricePaths(i, :) = pricePaths(i, :) .* exp( sqrt(deltaT) * w(i-1, :) .* ...
            currentBeta .* etaPaths(i-1, :) - (1 / 2) * (currentBeta).^2 .* ... 
            etaSquared(i-1, :) * deltaT );
    
    % update strikes
    bins(i, :) = mostProbableBins(numStrikes+1, pricePaths(i, :));
    averageStrikes(i, :) = genAverageBin(bins(i, :));
    
    
    % compute expectation
    
    % procede one strike at a time.
    %TODO: parallelize?
    
    numerator = zeros(1, numStrikes);
    denominator = zeros(1, numStrikes);
    
    for k=1:numStrikes
        dist = bins(i, k+1) - bins(i, k);
        deltaRow = delta(pricePaths(i, :), averageStrikes(i, k), dist);
        numerator(k) = etaSquared(i, :) * deltaRow';
        denominator(k) = sum(deltaRow);
    end
 
    expectation(i, :) = numerator ./ denominator;
    
    % update the betas
    % be careful about K dependence for sigma^2
    % TODO: is this guaranteed to be >= 0?
    beta(i, :) = sqrt(getDupireVol(i, averageStrikes(i, :)) ./ expectation(i, :));

    
    
end

end





% generates the grid using equidistant points between 0 and maxVal
function [bins] = equidistantBins(numBins, maxVal)

bins = zeros(1, numBins);
for i=1:numBins
    bins(i) = maxVal/(numBins-1)*(i-1);
end

end


function [avgBins] = genAverageBin(bins)

numAvgBins = length(bins) - 1; %-1 because we're only using avg points
avgBins = zeros(1, numAvgBins);
for i = 1:numAvgBins
    avgBins(i) = (1/2)*(bins(i) + bins(i+1));
end

end


% generates the grid using the most probable points.
% TODO: only run a few particles for this, just to get an approximation? 
%    this will speed up actual runtime (but not asymptotic runtime)
%    This might not be possible for multiple timesteps.
function [bins] = mostProbableBins(numBins, prices)

bins = zeros(1, numBins);

% f: CDF of prices evaluated at x.
% note that x and f are sorted.
[f, x] = ecdf(prices);

% if P = numParticles, B = numBins, then this takes O(B logP). Can this
%    be improved?
%    If we were to use a naive linear scan, we can do this in O(P).
% get bin cutoffs for CDFs from 0 to 1 in increments of 1/numBins.
bins(1) = x(1);
for i=1:numBins-1
    cdfProb = i/(numBins); 
    
    % binary search O(log n)
    lower = 1;
    upper = length(x);
    % +1 to allow termination when lower and upper are just 1 off
    while (lower+1) < upper
        % fix: round towards 0
        mid = fix(lower + (upper-lower)/2);
        if f(mid) < cdfProb
            lower = mid;
        elseif f(mid) > cdfProb
            upper = mid;
        else
            lower = mid;
            upper = mid;
        end
    end
    
    bins(i+1) = x(mid);
    
    
end



end



% kernel wrapper
% returns vector
function [output] = delta(prices, strike, dist)
% note that kernel takes in a vector as input.
% TODO: change dist to a vector.

output = 1/dist * kernel((prices - strike)/dist);

end


% "naive" kernel.
% output: 1 if x is in (-0.5, 0.5], 0 otherwise. 
% vectorized. 
function [n] = kernel(x)

n = x>-0.5 & x<= 0.5;

end


%vectorized Gaussian kernel. 
function [n] = gaussianKernel(x)
c = 1 / sqrt(2*pi);
n = c * exp(-1/2 * x.^2);
end


% OLD.
%TODO: delete?
function [avgStrikes] =genAverageStrikes(strikes)
% value at (T_i, S_j) = the strike for which beta(T_i, S_j) is computed
timeSteps = size(strikes, 1);
numAvgStrikes = size(strikes, 2) - 1; %-1 because we're only using avg points
avgStrikes = zeros(timeSteps, numAvgStrikes);
for row = 1:timeSteps
    for col = 1:numAvgStrikes
        avgStrikes(row, col) = (1/2)*(strikes(row, col) + strikes(row, col+1));
    end
end


end

function [vol] = getDupireVol(t, strikes)
    vol = ones(size(strikes)); % all 1's. Change this to the real dupire vol
end

function [betaInterp] = getBeta(prices, strikeAverages, betaRow)

% TODO: fix extrapolation. Right now, it extrapolates with the same
%    method as interpolation.
betaInterp = interp1(strikeAverages, betaRow, prices, 'linear', 'extrap');
% TODO:remove this
%betaSq = zeros(size(betaRow));
%for i=1:length(prices);   
%     leftIdx = search(prices(i));
%     betaSq(i) = interpolate(betaRow, leftIdx, price(i));
%end


end

% TODO: remove this if not necessary
function interpolate(betaRow, leftInx)



end

%TODO: remove this if not necessary
function search(x)

end

% drift factor
%TODO: what is this supposed to be?
function [lambda] = lambdaD(time, n)
% return vector
lambda = ones(1, n);

end


% OLD.
%TODO: delete?
function [strikes] = generateStrikes(t, n)

% for now, uniformly generate strikes between 0 and M
M = 10;

strikes = zeros(1, n);
for i=1:n
    strikes(i) = M/(n-1)*(i-1);
end
% uniform across time for now
strikes = repmat(strikes, t, 1);

end


function [w] = generateBrownian(t, N)

w = normrnd(0, 1, t, N);

end


function [etaPaths] = generateEtaPaths(timeSteps, deltaT, numEtaPaths, lambda)

% dimension: time x paths
W = generateBrownian(timeSteps, numEtaPaths);

% forward euler
Z = zeros(timeSteps, numEtaPaths);
etaPaths = zeros(timeSteps, numEtaPaths);

for n = 2:timeSteps
    Z(n, :) = Z(n-1, :) -lambda*Z(n-1, :)*deltaT + ... 
        gammaM((n-1)*deltaT)*sqrt(deltaT).*W(n, :);
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



end



function [n] = gammaM(t)

n = 1;

end
