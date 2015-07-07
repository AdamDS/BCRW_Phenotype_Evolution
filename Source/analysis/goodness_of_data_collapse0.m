%% goodness_of_data_collapse0.m
%
% [crit_exp,dexp] = goodness_of_data_collapse0(exp1,exp2,x,y,I,J,q,eta)
%
% Bhatarcharjee & Seno 2001: goodness of data-collapse
function [Pb] = goodness_of_data_collapse0(x,y,J,q), 
dj = J(2)-J(1);
Pb = 0;
Nover = 1;  %track number of pairs
pb = zeros(length(J));
D = 1:length(J);
d = 0;
for p = J, %need to interpolate
  d = d +1;
  cat_dy = [];
  for j = J, %current function (gives horizontal range)
    dy = [];
    if j==p,  
      cat_dy = cat_row(cat_dy,0);
%       fprintf('\t\t\tp = %1.0f\t\t\tj = %1.0f\n',p,j);
    else, 
      if numel(x{j})>0 && numel(x{p}),  %make sure there is data to work with
        interpmin = max([min(x{j}) min(x{p})]); %x greatest lower bound of interpolation range
        interpmax = min([max(x{j}) max(x{p})]); %x lowest upper bound of interpolation range
        xjk = x{j}(x{j}>interpmin & x{j}<interpmax);
        xpk = x{p}(x{p}>=interpmin & x{p}<=interpmax);
        if length(xjk)>2 && length(xpk)>2,  
          yint = interp1(xjk,y{j}(find(xjk)),xpk,'spline'); %interpolate y{i,k}
          dy = abs(y{p}(isnan(yint)==0)-yint(isnan(yint)==0)).^q;
          Nover = length(dy);
          cat_dy = cat_row(cat_dy,dy./Nover);
        else, 
          cat_dy = cat_row(cat_dy,0);
        end
      else
        fprintf('no data i=%1.0f\t\tj=%1.0f\t\tk=%1.0f\n',j,p);
      end
    end
  end
%   disp([p j; d c]);
%   fprintf('\t\t\t\t\t\trxc = %1.0f x %1.0f\n',size(cat_dy,1),size(cat_dy,2));
  if length(cat_dy(:,1)~=0)<length(J), 
    pb(:,d) = -1;
  else, 
    pb(:,d) = sum(cat_dy,2);  %sum i,over (x values)
  end
end
% goodness value for this pair of exp1 & exp2

Pb = sum(sum(pb,2)).^(1/q); %sum p( sum j~=p )
end