%% goodness_of_data_collapse.m
% Bhattacharjee & Seno 2001 arXiv: goodness of data-collapse
% [exp,Pb,iminPb] = goodness_of_data_collapse(exponents,x,y,J,q,eta)
%
% INPUTS:
% exponents = vector of attempted exponents (1D)
% x = universal scaling function argument values (horizontal axis)
% y = universal scaling function values (vertical axis)
% J = cells to collapse from x & y
% q = difference moment-like value (same as the q in the paper)
% eta = estimated error (same as the eta in the paper)
%
% OUTPUTS:
% exp = best [exponent-error exponent exponent+error]
% Pb = goodness of fit measure for each of the exponents
% iminPb = location of minimum in Pb
%
% -ADS 11*27*13
function [exp,Pb,iminPb] = goodness_of_data_collapse(exponents,x,y,J,q,eta), 
exp = [0 0 0];
C = 'bgrym';
color = [C C C C];
de = mean(diff(exponents));
oode = 1/de; %One Over Difference between Exponents
nexp = length(exponents);
dj = J(2)-J(1);
Pb = zeros(1,nexp);
% counter = 0;
for k = 1:nexp,  %for each exponent attempted
  spb = zeros(length(J),1);
  Nover = 0;  %track number of pairs
  pb = zeros(length(J));
  d = 0;
  for j = J, %need to interpolate (NOTE: j is p in Bhatacharjee & Seno)
    d = d +1;
    not_d = J(j~=J);
    dy = [];
    cat_dy = [];
    c = 0;
    for i = not_d, %(NOTE: i is j~=p in Bhatacharjee & Seno)
      c = c +1;
%       counter = debug_round(counter);
      if numel(x{i,k})>0 && numel(x{j,k}),  %make sure there is data to work with
        interpmin = max([min(x{i,k}) min(x{j,k})]); %x greatest lower bound of interpolation range
        interpmax = min([max(x{i,k}) max(x{j,k})]); %x lowest upper bound of interpolation range
%         counter = debug_round(counter);
        tik = find(x{i,k}>=interpmin & x{i,k}<=interpmax);
        xik = x{i,k}(tik);
        yik = y{i,k}(tik);
        tjk = find(x{j,k}>=interpmin & x{j,k}<=interpmax);
        xjk = x{j,k}(tjk);
        yint = interp1(xik,yik,xjk,'spline'); %interpolate y{i,k} at xjk points
        dy = abs(y{j,k}(isnan(yint)==0)-yint(isnan(yint)==0)).^q;
      else
        fprintf('no data i=%1.0f\t\tj=%1.0f\t\tk=%1.0f\n',i,j,k);
      end
      cat_dy = cat_row(cat_dy,dy./length(dy));  %length(dy) = Nover in paper
    end
    pb(not_d,d) = sum(cat_dy,2);  %sum i,over (x values)
  end
  Pb(k) = sum(sum(pb)).^(1/q); %sum p( sum j~=p )
end

% plot(exponents,Pb,'-x'); xlabel('critical exponent'); ylabel('Goodness measure');
% title('Goodness of data-collapse');

[minPb,iminPb] = min(Pb);

exp(2) = exponents(iminPb);

plus = (1+eta)*exp(2);
minus = (1-eta)*exp(2);
base_range = find(exponents<=(plus+de) & exponents>=(minus-de));
base_Pb = Pb(base_range);

if minus>exponents(1),  int_exp = [exponents(base_range(1)) minus ...
    exponents(base_range(2:(end-1)))]; mbuff = 1;
else, int_exp = exponents(1:end-1); mbuff = 0; end

if plus<exponents(end), int_exp = [int_exp plus ...
    exponents(base_range(end))]; pbuff = -1;
else, int_exp = [int_exp exponents(base_range(end))]; pbuff = 0;  end

int_Pb = interp1(exponents(base_range),base_Pb,int_exp,'spline');
exp(1) = exp(2)-eta*exp(2)*sqrt(2*log(int_Pb(1+mbuff)/minPb)).^-1;
exp(3) = exp(2)+eta*exp(2)*sqrt(2*log(int_Pb(end+pbuff)/minPb)).^-1;
% plus = (1+eta)*exp(2);
% rex = abs(exponents-round(plus*oode)/oode);
% Pbp = Pb(rex==min(rex)); %get nearest Pb(exponents) for plus
% exp(1) = exp(2)-eta*exp(2)*sqrt(2*log(Pbp(1)/minPb)).^-1;
% 
% minus = (1-eta)*exp(2);
% rex = abs(exponents-round(minus*oode)/oode);
% Pbm = Pb(rex==min(rex)); %get nearest Pb(exponents) for minus
% exp(3) = exp(2)+eta*exp(2)*sqrt(2*log(Pbm(1)/minPb)).^-1;

end