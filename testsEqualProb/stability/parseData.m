function [tensor, numParticles] = parseData(fileName, numStrikes)
  numItems = 19; % :%s/particles//gn  in vim to find this number
  headers = {'particles','strikes','time','zeros','numIters'};
  dataTags ={'expectation','strikes'};
  numHeaders = length(headers);
  numData = length(dataTags);
  fid = fopen(fileName, 'r');
 
  tensor = cell(numData, numItems);
  numParticles = zeros(1, numItems);

  for i=1:numItems
    % grab the first line
    strFormat = strcat(headers{1},':%f');
    data = textscan(fid, strFormat);
    numParticles(i) = data{1};

    % drop the next lines completely
    remainingHeaderLines = numHeaders-1; 
    data = textscan(fid, '%s', remainingHeaderLines, 'Delimiter', '\n');
    
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
