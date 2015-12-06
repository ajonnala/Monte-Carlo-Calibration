%>>>>>>>>>>>>>>>>>>>>>>>
%
% particle.m
% Nov 2016
%
% Runs the particle method.
% This version minimizes RAM usage
%
%
%<<<<<<<<<<<<<<<<<<<<<<<
function [expectation, averageStrikes, strikeBins] = particle(initialPrice, timeSteps, timeHorizon, ...
            numParticles, numStrikes )


% first, run the particle algorithm with no path resampling to get an
%   estimate for the bin cutoffs.
[~, ~, strikeBins] = runParticle(initialPrice, timeSteps, timeHorizon, ...
        numParticles, {1, numStrikes});
% now, run the particle algorithm again, but this time with the bin cutoffs
%   from before.
[expectation, averageStrikes, ~] = runParticle(initialPrice, timeSteps, timeHorizon, ...
        numParticles, {0, numStrikes, strikeBins});

end


% strikeCell: {generate (bool), numStrikes (int), bins (array)}.
% If generate == 0, then we're passing in bins that we want to use.
%   Otherwise, we are generating bins at each step according to numStrikes.
function [expectation, averageStrikes, bins] = runParticle(initialPrice, timeSteps, ...
        timeHorizon, numParticles, strikeCell )

% constant for the equidistant bins function. We can eventually get rid of
%   this.
M = 10;

generateOrNo = strikeCell{1};
numStrikes = strikeCell{2};
if (generateOrNo == 0)
    bins = strikeCell{3};
else
    % initialize bins
    bins = zeros(timeSteps, numStrikes+1);
    bins(1, :) = equidistantBins(numStrikes+1, M);
end


% parameters
lambda = 1; % for the eta paths



deltaT = 1/(timeSteps) * timeHorizon;



% particle step


expectation = zeros(timeSteps, numStrikes);
beta = zeros(timeSteps, numStrikes);
etaPath = ones(1, numParticles);
etaSquared = ones(1, numParticles);
% Z, V is used in the eta path
Z = zeros(1, numParticles);
V = zeros(timeSteps, 1);

% initialize strikes
averageStrikes = zeros(timeSteps, numStrikes);
averageStrikes(1, :) = genAverageBin(bins(1, :));

% initialize the first row of pricePaths --> current spot
pricePath = ones(1, numParticles) * initialPrice;
% initialize first row of expectation ---> 1
expectation(1, :) = 1;

% initialize initial beta
beta(1, :) = getDupireVol(0, averageStrikes(1, :));




for i=2:timeSteps
    w = generateBrownian(1, numParticles);
    %%% compute next price
    currentBeta = getBeta(pricePath, averageStrikes(i-1, :), beta(i-1, :));
    pricePath = pricePath .* lambdaD(i, numParticles);
    pricePath = pricePath .* exp( sqrt(deltaT) * w .* ...
            currentBeta .* etaPath - (1 / 2) * (currentBeta).^2 .* ...
            etaSquared * deltaT );

	clear w; % free up some ram

    %%% update strikes (and bins, if necessary)
    if (generateOrNo ~= 0)
        bins(i, :) = mostProbableBins(numStrikes+1, pricePath);
    end
    averageStrikes(i, :) = genAverageBin(bins(i, :));

    %%% update eta
    w = generateBrownian(1, numParticles);
    %TODO is this right?
    Z = Z -lambda*Z*deltaT + gammaM((i-1)*deltaT)*sqrt(deltaT).*w;
    clear w; % free up some ram
    V(i) = V(i-1) + (gammaM((i-1)*deltaT))^2 * exp(2*lambda) * deltaT;
    variance = exp(-2*lambda*i*deltaT) * V(i);
    etaPath = exp(Z - repmat(variance, 1, numParticles));
    etaSquared = etaPath .^2;

    %%% compute expectation

    % procede one strike at a time.
    %TODO: parallelize?

    numerator = zeros(1, numStrikes);
    denominator = zeros(1, numStrikes);

    for k=1:numStrikes
        dist = bins(i, k+1) - bins(i, k);
        deltaRow = delta(pricePath, averageStrikes(i, k), dist);
        numerator(k) = etaSquared * deltaRow';
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

% Note that as x is sorted, then the jth percentile value is at index
%    (j/100)*numParticles.
% We ignore everything after the 99th percentile.
bins(1) = x(1);
maxIndex = ceil(0.99*size(prices, 2));
bins(numBins) = x(maxIndex);

for i=2:numBins-1
    cdfProb = (i-1)/(numBins - 1);
    index = ceil(cdfProb * maxIndex);
    bins(i) = x(index);
end


end



% kernel wrapper
% returns vector
function [output] = delta(prices, strike, dist)
% note that kernel takes in a vector as input.
% TODO: change dist to a vector.

output = 1/dist *hatKernel((prices - strike)/dist);

end


% "naive" kernel.
% output: 1 if x is in (-0.5, 0.5], 0 otherwise.
% vectorized.
function [n] = kernel(x)

n = x>-0.5 & x<= 0.5;

end

function [n] = hatKernel(x)
ind = (abs(x) >= .5);
ind2 = (abs(x) < .5);
x(ind) = 0;
x(ind2) = -abs(2*x(ind2)) + 1;
n = x;
end

function [n] = epanKernel(x)
ind = (abs(x) <= 1);
ind2 = (abs(x) > 1);
x(ind2) = 0;
x(ind) = (3/4) * ( 1 - x(ind) .^2);
n = x;
end

function [n] = quarticKernel(x)
ind = (abs(x) <= 1);
ind2 = (abs(x) > 1);
x(ind2) = 0;
x(ind) = (15/16)* ( 1 - x(ind) .^2).^2;
n = x;
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


% used in eta path
function [n] = gammaM(t)

n = 1;

end
