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

% numStrikes-1 because we're evaluating at midpoints
expectation = zeros(timeSteps, numStrikes-1);
beta = zeros(timeSteps, numStrikes-1);
pricePaths = zeros(timeSteps, numParticles);

strikes = generateStrikes(timeSteps, numStrikes);
averageStrikes = genAverageStrikes(strikes);

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
    
    % compute expectation
    
    % procede one strike at a time.
    %TODO: parallelize?
    
    numerator = zeros(1, numStrikes-1);
    denominator = zeros(1, numStrikes-1);
    
    for k=1:numStrikes-1
        dist = strikes(i, k+1) - strikes(i, k);
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




% generates the grid using the most probable points.
% TODO: only run a few particles for this, just to get an approximation? 
%    this will speed up actual runtime (but not asymptotic runtime)
function [bins] mostProbableBins()


end





% kernel wrapper
% returns vector
function [output] = delta(prices, strike, dist)
% note that kernel takes in a vector as input.
% TODO: change dist to a vector.

output = 1/dist * kernel((prices - strike)/dist);

end


% "naive" kernel.
% vectorized. 
function [n] = kernel(x)

n = x>-0.5 & x<= 0.5;

end


%vectorized Gaussian kernel. 
function [n] = gaussianKernel(x)
c = 1 / sqrt(2*pi);
n = c * exp(-1/2 * x.^2);
end


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
