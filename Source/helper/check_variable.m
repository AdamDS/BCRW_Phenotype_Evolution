%% check_variable.m
% function [cdim,cmin,cmax] = check_variable(go,variable,standards), 
% Checks the dimensions, minimum, & maximum of a variable according to standards.
function [crow,ccol,cmin,cmax] = check_variable(go,variable,standards), 
crow = 0;
ccol = 0;
cmin = 0;
cmax = 0;
if go==1, 
  srow = standards(1);
  scol = standards(2);
  smin = standards(3);
  smax = standards(4);

  [row,col] = size(variable); 
  if srow,  if srow==row, crow = 1; end, end
  if scol,  if scol==col, ccol = 1; end, end
  if smin,  minv = min(min(variable));  if smin==minv, cmin = 1; end, end
  if smin,  maxv = max(max(variable));  if smax==maxv, cmax = 1; end, end
end
      
end