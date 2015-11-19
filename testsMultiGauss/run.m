function runParticleAndLog(initialPrice, timeSteps, timeHorizon, ...
    numStrikes)

%  initialPrice = 0.5;
%  timeSteps = 2; % keep this as 2 for a one-period model
%  timeHorizon = 1;
  %numParticles = 1000000;
  %numStrikes = 50;

  file = fopen('particle_100m_5_steps.txt', 'a');
  
  
  
numParticles = 100000000;
for initialPrice = 0.2:0.2:4
       tStart = cputime;
        [expectation, strikes] = particle(initialPrice, timeSteps, timeHorizon, ...
          numParticles, numStrikes);
        tEnd = cputime;
       
        

        totalTime = (tEnd - tStart);
        totalZero = sum(sum(isnan(strikes))); 
      
      
      fprintf(file, 'particles:%u\nstrikes:%u\ntime:%f\nzeros:%d\ninitPrice:%f\n', ...
              numParticles, numStrikes, totalTime, totalZero, initialPrice);
      % now print out num and denom

      fprintf(file, 'expectation:\n');
  
    for t=2:timeSteps 
      for i=1:numStrikes
        fprintf(file, '%.4f ', expectation(t, i));
      end
      fprintf(file, '\n');
    end

    fprintf(file, 'strikes:\n');
    
    for t=2:timeSteps
      for i=1:numStrikes
        fprintf(file, '%.4f ', strikes(t, i));
      end
      fprintf(file, '\n');
    end    
      fprintf(file, '\n');

end 
    
  fclose(file);
  exit;
  
end
