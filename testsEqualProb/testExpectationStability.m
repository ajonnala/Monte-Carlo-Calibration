function testExpectationStability(initialPrice, timeSteps, timeHorizon, ...
    numStrikes)

%  initialPrice = 0.5;
%  timeSteps = 2; % keep this as 2 for a one-period model
%  timeHorizon = 1;
  %numParticles = 1000000;
  %numStrikes = 50;

  filePath = 'expStabilityLog';
  filePath = strcat(filePath, num2str(numStrikes));
  filePath = strcat(filePath, '.txt');
  file = fopen(filePath, 'a');
  
  numIters = 10;
  
  for numParticles = 100000:100000:10000000
    if((numParticles > 1000000) && (mod(numParticles, 1000000) ~= 0))
        continue;
    end
        E = zeros(numIters, numStrikes);
        S = zeros(numIters, numStrikes);
       tStart = cputime;
       for i=1:numIters
        [expectation, strikes] = particle(initialPrice, timeSteps, timeHorizon, ...
          numParticles, numStrikes, numParticles);
         E(i, :) = expectation(2, :);
         S(i, :) = strikes(2, :);
       end
       tEnd = cputime;
       
        

        totalTime = (tEnd - tStart);
        totalZero = sum(sum(isnan(S))); 
      
      
      fprintf(file, 'particles:%u\nstrikes:%u\ntime:%f\nzeros:%d\nnumIters:%d\n', ...
              numParticles, numStrikes, totalTime, totalZero, numIters);
      % now print out num and denom

      fprintf(file, 'expectation:\n');
  
      fprintf(file, [repmat('%.4f ', 1, size(E, 2)) '\n'], E);
      
      fprintf(file, '\nstrikes:\n');
      
      fprintf(file, [repmat('%.4f ', 1, size(E, 2)) '\n'], S);
      
      fprintf(file, '\n\n');

    
    
  end
  fclose(file);
  exit;
  
end
