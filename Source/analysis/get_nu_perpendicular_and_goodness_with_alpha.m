%% Determine correlation time exponent: rho*t^alpha = f(t*Delta^nu_para)
these = find(Delta~=0);
t = 1:NGEN;  %maximum time span to fit
limited = fSmM(1):fSmM(2);
% limited = 5:1000;

dnr = 0.01;
nu_perps = 0.58:dnr:0.89;  %expect 0.734 ~ 0.73
Nnr = length(nu_perps);

x = cell(M,Nnr);
y = cell(M,Nnr);

thisp = zeros(1,2);

Pb = zeros(Q,Nnr);
best_Pb = zeros(Q,1);
best_nrb = zeros(Q,3);
best_nra = zeros(Q,3);

fignum = fignum +1;
leg = cell(Q,1);
legm = cell(M,1);

for q = 1:Q,  
  alpha = mean(alpha_crit(:));%alpha_crit(end);
%   nu_para = mean([nu_parab(:,2)]);% nu_paraa(end,2)]);
  nu_para = mean(nu_parac);
  for j = 1:Nnr,  %nu_paras along columns
    nu_perp = nu_perps(j);
    % construct attempted time series for data-collapse
    for m = 1:M,  %construct x & y to fit
      L = min([xl(m) yl(m)]);
      if length(RHOX{m})<max(limited), 
        lim = 1:length(RHOX{m});
      else, 
        lim = limited;  
      end
      tau = t(t>lim(1) & t<=lim(end));  %set time span to fit
      x{m,j} = tau(RHOX{m}(tau)~=0)./(L.^(nu_perp./nu_para));
      y{m,j} = RHOX{m}(tau).*(tau(RHOX{m}(tau)~=0).^alpha);
      legm{m} = int2str(L);
    end
%         fprintf('k = %1.0f\t\tj = %1.0f\n',k,j);
  end
  J = 1:M;
  [nr(q,:),Pb(q,:),iminPb(q)] = goodness_of_data_collapse(nu_perps,x,y,J,q,eta);
  best_Pb(q) = Pb(q,iminPb(q));
  fprintf('q = %1.0f\t\teta = %1.3f\n',q,eta);
  fprintf('\t\t-nu_perp = %1.2f\t\tnu_perp = %1.2f\t\t+nu_perp = %1.2f\n', ...
    nr(q,1),nr(q,2),nr(q,3));
  
  if do_plot_Pb,  
    figure(fignum); hold on;
    plot(nu_perps,Pb(q,:),c(q));
    leg{q} = int2str(q);
  end
end

if do_plot_Pb,  
  legend(leg)
  title('Finite-size goodness of fit');
  xlabel('\nu_{\perp}');
  ylabel('Pb');
end

disp([best_Pb nr(:,2)]);

fprintf('Best nu_perps with %1.0f percent estimated errors\n',eta*100);
[~,best_q] = min(best_Pb);
nu_perp = nr(best_q,:);
disp(nu_perp);

%% Plot the best exponents overall Pb of q
if do_plot_collapse,  
  lmtd = 1:NGEN;
  fignum = fignum +1;
  z = nu_perp./nu_para;
  plot_finite_size([z alpha],RHOX,J,min([xl;yl]),t,lmtd,fignum);
  set(gca,'xscale','log','yscale','log');
  xlabel('log_{10}( t L^{z} )');
  ylabel('log_{10}( \rho t^\alpha )');
  title(['\alpha = ' num2str(alpha) ' & \nu_{\perp} = ' num2str(nr(best_q,2))]);
  legend(legm);
end

write_nu_perp;

% for i = Jb,  %for each time series
%   if length(RHO{i})<max(limited), lim = limited(1):length(RHO{i});
%   else, lim = limited;  end
%   tau = t(t>lim(1) & t<=lim(end));  %set time span to fit
%   x{i} = tau(RHO{i}(tau)~=0).*(abs(Delta(i)).^b(2)); %control parameter x critical exponent
%   y{i} = RHO{i}(tau).*(tau(RHO{i}(tau)~=0).^expb(1));
%   fitted = find(tau==tfit(1) | tau==tfit(2));
%   mm = x{i}(fitted); Mm = y{i}(fitted);
%   plot([mm(1) mm(1)],[1 Mm(1)],'k');
%   plot([mm(2) mm(2)],[1 Mm(2)],'k');
% end
% if exist('Ja','var')>0,  
%   for i = Ja,  %for each time series
%     if length(RHO{i})<max(limited), lim = limited(1):length(RHO{i});
%     else, lim = limited;  end
%     tau = t(t>lim(1) & t<=lim(end));  %set time span to fit
%     x{i} = tau(RHO{i}(tau)~=0).*(abs(Delta(i)).^expa(2)); %control parameter x critical exponent
%     y{i} = RHO{i}(tau).*(tau(RHO{i}(tau)~=0).^expa(1));
%     fitted = find(tau==tfit(1) | tau==tfit(2));
%     mm = x{i}(fitted); Mm = y{i}(fitted);
%     plot([mm(1) mm(1)],[1 Mm(1)],'k');
%     plot([mm(2) mm(2)],[1 Mm(2)],'k');
%   end
% end