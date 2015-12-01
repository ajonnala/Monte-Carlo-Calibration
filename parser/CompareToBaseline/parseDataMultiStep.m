function [tensor] = parseDataMultiStep(fileName, numStrikes)
  % NOTE: numItems is 9 for the 10m particle data, and 7 for the incomplete
  %     100m data.l
  numItems = 7; % :%s/particles//gn  in vim to find this number
  headers = {'particles','strikes','time','zeros','initPrice','timeSteps'};
  dataTags ={'expectation','strikes'};
  numHeaders = length(headers);
  numData = length(dataTags);
  fid = fopen(fileName, 'r');
 
  tensor = cell(numData, numItems);
  %numTimeSteps = zeros(1, numItems);

  for i=1:numItems
      
    % drop the next lines completely
    initialHeaderLines = numHeaders-1; 
    data = textscan(fid, '%s', initialHeaderLines, 'Delimiter', '\n');
    
    % grab the last line
    strFormat = strcat(headers{6},':%f');
    data = textscan(fid, strFormat);
    %numTimeSteps(i) = data{1};
    
    
    %get the data
    for j=1:numData
      % clear data tag
      data = textscan(fid, '%s', 1, 'Delimiter', '\n');
      strFormat = repmat('%f', 1, numStrikes);
      data = textscan(fid, strFormat, 'Delimiter', ' ');
      tensor{j, i} = cell2mat(data);
       
      % note: no need to clear out the trailing newlines here.
    end

  end

  fclose(fid);

end
