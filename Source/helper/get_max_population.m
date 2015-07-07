
function [maxN] = get_max_population(),
global SIMOPTS;
basic_map_size = SIMOPTS.basic_map_size;
if ~exist('SIMOPTS.op'), op = 0.25;  else, op = SIMOPTS.op;  end
[~,X,Y] = landscape_measures(basic_map_size);
h = op*sqrt(3)/2;
if X>=Y,  
  max_width = floor(X/op) +1; %209
  max_height = floor(Y/h) +1; %209
  if mod(max_height,2)==1,  
    long_rows = max_height/2;
    short_rows = long_rows;
  else, 
    long_rows = ceil(max_height/2); %105
    short_rows = floor(max_height/2); %104
  end
else, 
  max_width = floor(Y/op) +1;
  max_height = floor(Y/h) +1;
  if mod(max_height,2)==1,  
    long_rows = max_height/2;
    short_rows = long_rows;
  else, 
    long_rows = ceil(max_height/2);
    short_rows = floor(max_height/2);
  end
end
nL = max_width*long_rows;
nS = (max_width-1)*short_rows;
maxN = floor(nL +nS);
% fprintf('%1.0f\n',maxN);
end