%% goodness_of_data_collapse2.m
%
% [crit_exp,dexp] = goodness_of_data_collapse2(exp1,exp2,x,y,I,J,q,eta)
%
% Bhatarcharjee & Seno 2001: goodness of data-collapse
function [crit_exp,dexp,global_min_iPb] = goodness_of_data_collapse2(exp1,exp2,x,y,J,q,eta), 
crit_exp = [0 0];  dexp = [0 0; 0 0];
C = 'bgrym';
color = [C C C C];
oode1 = 1/mean(diff(exp1)); %One Over Difference between Exponents
oode2 = 1/mean(diff(exp2)); %One Over Difference between Exponents
nexp1 = length(exp1);
nexp2 = length(exp2);
dj = J(2)-J(1);
Pb = zeros(nexp1,nexp2);
figure(5150); hold off;
for e1 = 1:nexp1, 
  for e2 = 1:nexp2,  %for each exponent attempted
    spb = zeros(length(J),1);
    Nover = 0;  %track number of pairs
    pb = zeros(length(J)-1);
    d = 0;
    for p = J, %need to interpolate
      d = d +1;
      dy = [];
      cat_dy = [];
      c = 0;
      for j = J(p~=J), %current function (gives horizontal range)
        if numel(x{j,e1,e2})>0 && numel(x{p,e1,e2}),  %make sure there is data to work with
          c = c +1;
          interpmin = max([min(x{j,e1,e2}) min(x{p,e1,e2})]); %x greatest lower bound of interpolation range
          interpmax = min([max(x{j,e1,e2}) max(x{p,e1,e2})]); %x lowest upper bound of interpolation range
          xjk = x{j,e1,e2}(x{j,e1,e2}>interpmin & x{j,e1,e2}<interpmax);
          xpk = x{p,e1,e2}(x{p,e1,e2}>=interpmin & x{p,e1,e2}<=interpmax);
          yint = interp1(xjk,y{j,e1,e2}(find(xjk)),xpk); %interpolate y{i,k}
          dy = abs(y{p,e1,e2}(isnan(yint)==0)-yint(isnan(yint)==0)).^q;
          cat_dy = cat_row(cat_dy,dy./length(dy));
        else
          fprintf('no data i=%1.0f\t\tj=%1.0f\t\tk=%1.0f\n',j,p,e1,e2);
        end
      end
      pb(:,d) = sum(cat_dy,2);  %sum i,over (x values)
    end
    % goodness value for this pair of exp1 & exp2
    Pb(e1,e2) = sum(sum(pb,2)).^(1/q); %sum p( sum j~=p )
  end
  [minPb(e1),iminPb(e1)] = min(Pb(e1,:)); %minimum Pb for this exp1
  best_exp2(e1) = exp2(iminPb(e1));
end

figure(5150);
plot(exp1,best_exp2,'-x');


figure(5151); surf(Pb); 
xlabel('critical exponent (row)'); ylabel('critical exponent (col)');  zlabel('Goodness measure');
title('Goodness of data-collapse');

[global_min_Pb,global_min_iPb(2)] = min(minPb);
global_min_iPb(1) = iminPb(global_min_iPb(2));

% figure(5152); plot3(exp1,exp2,iminPb);

%% exponent 1
crit_exp(1) = exp1(global_min_iPb(1));

plus = (1+eta)*crit_exp(1);
rex = abs(exp1-round(plus*oode1)/oode1);
Pbp = Pb(rex==min(rex),global_min_iPb(2)); %get nearest Pb(exponents) for plus
dexp(1,1) = eta*crit_exp(1)*sqrt(2*log(Pbp(1)/global_min_Pb)).^-1;

minus = (1-eta)*crit_exp(1);
rex = abs(exp1-round(minus*oode1)/oode1);
Pbm = Pb(rex==min(rex),global_min_iPb(2)); %get nearest Pb(exponents) for minus
dexp(2,1) = eta*crit_exp(1)*sqrt(2*log(Pbm(1)/global_min_Pb)).^-1;

%% exponent 2
crit_exp(2) = exp2(global_min_iPb(2));

plus = (1+eta)*crit_exp(2);
rex = abs(exp2-round(plus*oode2)/oode2);
Pbp = Pb(global_min_iPb(1),rex==min(rex)); %get nearest Pb(exponents) for plus
dexp(1,2) = eta*crit_exp(2)*sqrt(2*log(Pbp(1)/global_min_Pb)).^-1;

minus = (1-eta)*crit_exp(1);
rex = abs(exp2-round(minus*oode2)/oode2);
Pbm = Pb(global_min_iPb(1),rex==min(rex)); %get nearest Pb(exponents) for minus
dexp(2,2) = eta*crit_exp(2)*sqrt(2*log(Pbm(1)/global_min_Pb)).^-1;
end

% numerical errors from rounding
% for i=1:length(np),
% cep(i) = crit_exp(i)*(1+eta);
% cem(i) = crit_exp(i)*(1-eta);
% nearp(i) = round(cep(i)*oode)/oode;
% nearm(i) = round(cem(i)*oode)/oode;
% rounddiffp(i) = (nearp(i)-crit_exp(i))/crit_exp(i);
% rounddiffm(i) = (nearm(i)-crit_exp(i))/crit_exp(i);
% end
% plot(np,rounddiffp,'-xr',np,rounddiffm,'-xb');
% hold on;
% plot(np,eta*ones(size(np)),'k',np,-eta*ones(size(np)),'k')
% the red & blue lines are about |eta|, so deviations introduce numerical
% errors depending on the crit_exp (horizontal axis) found 


%% Old code:
% for k = 1:nnp,  %for each nu_para attempted
%   spb = zeros(length(J),1);
% %   figure(7357); hold off;
%   Nover = 0;  %track number of pairs
%   pb = zeros(length(J)-1);
%   d = 0;
%   for j = J, %need to interpolate (NOTE: j is p in Bhatacharjee & Seno)
%     d = d +1;
% %     dy = cell(length(J)-1,1);
%     dy = [];
%     cat_dy = [];
%     c = 0;
%     for i = J(j~=J), %(NOTE: i is j~=p in Bhatacharjee & Seno)
% %       if i~=j,  
% %         Nover = Nover +1;
%         if numel(x{i,k})>0 && numel(x{j,k}),  %make sure there is data to work with
%           c = c +1;
% %           if min(T([i j],1))==T(i,1), tail_min = T(i,1);
% %           else, tail_min = T(j,1);  end
% %           if max(T([i j],2))==T(i,2), tail = tail_min:T(i,2); 
% %           else, tail = tail_min:T(j,2); end
%           interpmin = max([min(x{i,k}) min(x{j,k})]); %x greatest lower bound of interpolation range
%           interpmax = min([max(x{i,k}) max(x{j,k})]); %x lowest upper bound of interpolation range
%           %get range of x{i,k} within bounds
% %           mxik = x{i,k}(x{i,k}>=interpmin);
% %           Mxik = x{i,k}(x{i,k}<=interpmax);
%           %get range of x{j,k} within bounds
% %           mxjk = x{j,k}(x{j,k}>=interpmin);
% %           Mxjk = x{j,k}(x{j,k}<=interpmax);
%           xik = x{i,k}(x{i,k}>interpmin & x{i,k}<interpmax);
%           xjk = x{j,k}(x{j,k}>=interpmin & x{j,k}<=interpmax);
%           yint = interp1(xik,y{i,k}(find(xik)),xjk); %interpolate y{i,k}
%           dy = abs(y{j,k}(isnan(yint)==0)-yint(isnan(yint)==0)).^q;
% %           loglog(x{i,k},y{i,k},'k',x{j,k},yint,color(c)); hold on;
% %           disp([size(x{j,k}); size(yint); size(x{i,k}); size(y{i,k}); size(dy{c})]);
% %           fprintf('yint # !NaN = %1.0f\t\t y{%1.0f,%1.0f} # !NaN = %1.0f\n',...
% %             length(~isnan(yint)),j,k,length(~isnan(y{j,k})));
%         else
%           fprintf('no data i=%1.0f\t\tj=%1.0f\t\tk=%1.0f\n',i,j,k);
%         end
% %       end
%       cat_dy = cat_row(cat_dy,dy./length(dy));
%     end
% %     spb = sum(cat(1,dy{:}),2))./Nover;  %sum i,over (x values)
% %     fprintf('\t\tpb # ~NaN = %1.0f\n',length(~isnan(mpb)));
% %     pb(d) = sum(spb(isnan(spb)==0));
%     pb(:,d) = sum(cat_dy,2);  %sum i,over (x values)
%   end
% %   Pb(k) = sum(sum(pb,2)./Nover).^(1/q); %sum p( sum j~=p )
%   Pb(k) = sum(sum(pb,2)).^(1/q); %sum p( sum j~=p )
% end