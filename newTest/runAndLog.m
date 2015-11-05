function runAndLog()

  initialPrice = .5;
  %timeSteps = 2; % keep this as 2 for a one-period model
  timeHorizon = 1;
  %numParticles = 1000000;
  %numStrikes = 50;

  file = fopen('particle_10.txt', 'a');
  
 
  numStrikes = 10; 
  
    for numParticles = 100000:100000:10000000
%for numStrikes = 10:10:1000
      % for large strike prices, increase the jump size
      %if((numStrikes > 100) && (mod(numStrikes, 100) ~= 0))
        if((numParticles > 1000000) && (mod(numParticles, 1000000) ~= 0))
          continue;
        %end
      end

      % run 10 trials each
      %totalTime = 0;
      %totalNan = 0;
      %for i=1:10
        tStart = cputime;
        [numerator, denominator] = particle(initialPrice, timeHorizon, ...
          numParticles, numStrikes, numParticles);
        tEnd = cputime;
       
        

        totalTime = (tEnd - tStart);
        totalZero = sum(denominator(:) == 0); 
        %totalNan = totalNan + sum(sum(isnan(expectation)));
      %end
      
      
      fprintf(file, 'particles:%u\nstrikes:%u\ntime:%f\nzeros:%d\n', ...
              numParticles, numStrikes, totalTime, totalZero);
      % now print out num and denom

      fprintf(file, 'numerator:\n');
  
      for i=1:numStrikes-1
        fprintf(file, '%.4f  ', numerator(i));
      end
      fprintf(file, '\ndenominator:\n');
      for i=1:numStrikes-1
        fprintf(file, '%.4f  ', denominator(i));
      end    
      fprintf(file, '\n\n');

    %end
  end
  fclose(file);
  exit;
  
end
