function runParticleAndLog()

  initialPrice = 0.5;
  timeSteps = 2; % keep this as 2 for a one-period model
  timeHorizon = 1;
  %numParticles = 1000000;
  %numStrikes = 50;

  file = fopen('particle_log.txt', 'a');
  
  
  
  for numParticles = 100:100:100000000
    for numStrikes = 10:10:10000
      % run 10 trials each
      totalTime = 0;
      totalNan = 0;
      for i=1:10
        tStart = cputime;
        expectation = particle(initialPrice, timeSteps, timeHorizon, ...
          numParticles, numStrikes, numParticles);
        tEnd = cputime;
        
        totalTime = totalTime + (tEnd - tStart);
        totalNan = totalNan + sum(sum(isnan(expectation)));
      end
      
      
      fprintf(file, 'particles:%u\nstrikes:%u\ntime:%f\nNaN:%d\n\n', ...
              numParticles, numStrikes, totalTime, totalNan);
    
    end
  end
  fclose(file);
  exit;
  
end
