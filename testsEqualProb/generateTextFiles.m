function generateTextFiles(initialPrice, timeSteps, timeHorizon)

strikes = [3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 1000];
len = length(strikes);

for i=1:len
    
    parameters = '';
    parameters = strcat(parameters, num2str(initialPrice));
    parameters = strcat(parameters, ', ');
    parameters = strcat(parameters, num2str(timeSteps));
    parameters = strcat(parameters, ', ');
    parameters = strcat(parameters, num2str(timeHorizon));
    parameters = strcat(parameters, ', ');
    parameters = strcat(parameters, num2str(strikes(i)));
    command = 'runParticleAndLog(';
    command = strcat(command, parameters);
    command = strcat(command, ');');
  
    filePath = 'command';
    filePath = strcat(filePath, num2str(strikes(i)));
    filePath = strcat(filePath, '.txt');
    file = fopen(filePath, 'a');
    
    fprintf(file, command);
    fclose(file);
  

end