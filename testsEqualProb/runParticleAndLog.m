function runParticleAndLog(initialPrice, timeSteps, timeHorizon, ...
    numStrikes)

%  initialPrice = 0.5;
%  timeSteps = 2; % keep this as 2 for a one-period model
%  timeHorizon = 1;
  %numParticles = 1000000;
  %numStrikes = 50;

  filePath = 'particle_log';
  filePath = strcat(filePath, num2str(numStrikes));
  filePath = strcat(filePath, '.txt');
  file = fopen(filePath, 'a');
  
  
  
  for numParticles = 100000:100000:10000000
    if((numParticles > 1000000) && (mod(numParticles, 1000000) ~= 0))
        continue;
    end

       tStart = cputime;
        [expectation, strikes] = particle(initialPrice, timeSteps, timeHorizon, ...
          numParticles, numStrikes, numParticles);
        tEnd = cputime;
       
        

        totalTime = (tEnd - tStart);
        totalZero = sum(sum(isnan(strikes))); 
      
      
      fprintf(file, 'particles:%u\nstrikes:%u\ntime:%f\nzeros:%d\n', ...
              numParticles, numStrikes, totalTime, totalZero);
      % now print out num and denom

      fprintf(file, 'expectation:\n');
  
      for i=1:numStrikes
        fprintf(file, '%.4f ', expectation(2, i));
      end
      fprintf(file, '\nstrikes:\n');
      for i=1:numStrikes
        fprintf(file, '%.4f ', strikes(2, i));
      end    
      fprintf(file, '\n\n');

    
    
  end
  fclose(file);
  exit;
  
end
