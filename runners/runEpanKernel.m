function runEpanKernel()


  timeHorizon = 50;
  numStrikes = 50;
  timeSteps = 261;
  initialPrice = .5;
  numParticles = 1000000;
  s = strcat('~/private/Monte-Carlo-Calibration/Data/Epan_TH:',int2str(timeHorizon),'strikes:',int2str(numStrikes),'timesteps:',int2str(timeSteps),'numParticles:',int2str(numParticles),'.txt');
  file = fopen(s, 'a');

%for initialPrice = 0.2:0.2:4
       tStart = cputime;
        [expectation, strikes] = particleEpan(initialPrice, timeSteps, timeHorizon, ...
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

%end

  fclose(file);
  exit;

end
