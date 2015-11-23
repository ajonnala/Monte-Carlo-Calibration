function [dupireData] = sheetParse(fileName, numDates)
  % 41 strikes
  % 260 dates
  numStrikes = 41;
  
  M = dlmread(fileName, ',', 1, 1);
  
  time = M(1:numDates, 1)';
  fwd = M(1:numDates, 2)';
  atmVol = M(1:numDates, 4)';
  
  strikeStdDev = dlmread(fileName, ',', [0, 5, 0, 4+numStrikes]);
  
  data = dlmread(fileName, ',', [1, 5, numDates, 4 + numStrikes]);
  
  strikes = repmat(fwd', 1, numStrikes) .* exp(repmat(strikeStdDev, numDates, 1) ...
      .* repmat(atmVol', 1, numStrikes) .* repmat(sqrt(time'), 1, numStrikes));
  
  dupireData = cell(3, 1);
  dupireData{1} = time;
  dupireData{2} = strikes;
  dupireData{3} = data;
  

end