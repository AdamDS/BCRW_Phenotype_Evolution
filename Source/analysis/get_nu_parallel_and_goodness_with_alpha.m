%% Determine correlation time exponent: rho*t^alpha = f(t*Delta^nu_para)
these = find(abs(Delta)>=0.0001);
t = 1:NGEN;  %maximum time span to fit
if bms==8,  
  limitedb = fTmM(1):fTmM(2);
  limiteda = fTmM(1):(fTmM(2));
elseif bms==10, 
  limitedb = fTmM(1):(fTmM(2));
  limiteda = fTmM(1):(fTmM(2)+100);
elseif bms==12, 
  limitedb = fTmM(1):(fTmM(2));
  limiteda = fTmM(1):(fTmM(2)+200);
end

% limited = limitedb;

% eta = 0.03;  %estimated 10% error bars on exponents

dnp = 0.01;
nu_paras = 1.00:dnp:1.60;  %expect 1.295 ~ 1.30
Nnp = length(nu_paras);

x = cell(N,Nnp);
y = cell(N,Nnp);

alpha = alpha_crit(m);

thisp = zeros(1,2);

Q = 3;
Pbb = zeros(Q,Nnp);
Pba = zeros(Q,Nnp);
best_Pbb = zeros(Q,1);
best_Pba = zeros(Q,1);
best_npb = zeros(Q,3);
best_npa = zeros(Q,3);

fignum = fignum +1;
leg = cell(Q,1);

for q = 1:Q,  
  for j = 1:Nnp,  %nu_paras along columns
    nu_para = nu_paras(j);
    % construct attempted time series for data-collapse
    for i = 1:N,  %for each time series
      if i<iscrit, 
        % determine available time span
        if length(RHO{i})<max(limitedb), lim = limitedb(1):length(RHO{i});
        else, lim = limitedb;  end
      elseif i>iscrit,  
        if length(RHO{i})<max(limiteda), lim = limiteda(1):length(RHO{i});
        else, lim = limiteda;  end
      end
      % no need for collapsing on critical
      if i~=iscrit, 
        tau = t(t>lim(1) & t<=lim(end));  %set time span to fit
        x{i,j} = tau(RHO{i}(tau)~=0).*(abs(Delta(i)).^nu_para); %control parameter x critical exponent
        y{i,j} = RHO{i}(tau).*(tau(RHO{i}(tau)~=0).^alpha);
      end
    end
%     fprintf('k = %1.0f\t\tj = %1.0f\n',k,j);  %debugging print out
  end
  
  % get below critical goodness measure
  Jb = 1:(lessthan(end-offset));
  [best_npb(q,:),Pbb(q,:),iminPbb(q)] = goodness_of_data_collapse(nu_paras,x,y,Jb,q,eta);
  % get above critical goodness measure
  if greaterthan(offset+1)~=greaterthan(end), %if enough data exists
    Ja = (greaterthan(offset+1)):greaterthan(end);
    [best_npa(q,:),Pba(q,:),iminPba(q)] = goodness_of_data_collapse(nu_paras,x,y,Ja,q,eta);
  end
  best_Pbb(q) = Pbb(q,iminPbb(q));
  best_Pba(q) = Pba(q,iminPba(q));
  comb_Pb(q,:) = Pbb(q,:)+Pba(q,:);
  [best_comb_Pb(q),imincombPb(q)] = min(comb_Pb(q,:));
  
  fprintf('q = %1.0f\t\t\teta = %1.2f\t\t\tbms = %1.0f\n',q,eta,bms);
  fprintf('\t\t-nu_parab = %1.2f\t\tnu_parab = %1.2f\t\t+nu_parab = %1.2f\n', ...
    best_npb(q,1),best_npb(q,2),best_npb(q,3));
  fprintf('\t\t-nu_paraa = %1.2f\t\tnu_paraa = %1.2f\t\t+nu_paraa = %1.2f\n', ...
    best_npa(q,1),best_npa(q,2),best_npa(q,3));
  
  if do_plot_Pb,  
    figure(fignum);  
    subplot(3,1,1); plot(nu_paras,Pbb(q,:),c(q));  hold on;  
    subplot(3,1,2); plot(nu_paras,Pba(q,:),c(q));  hold on;  
    subplot(3,1,3); plot(nu_paras,comb_Pb(q,:),c(q));  hold on;  
    leg{q} = int2str(q);
  end
end

if do_plot_Pb,  
  figure(fignum);  subplot(3,1,1);  title('Pb below critical');
  xlabel('\nu_{||}');  ylabel('Pb');  legend(leg);
  figure(fignum);  subplot(3,1,2);  title('Pb above critical');
  xlabel('\nu_{||}');  ylabel('Pb');  legend(leg);
  figure(fignum);  subplot(3,1,3);  title('combined Pb critical');
  xlabel('\nu_{||}');  ylabel('Pb');  legend(leg);
end

[~,best_qb] = min(best_Pbb);
[~,best_qa] = min(best_Pba);
[~,best_comb_q] = min(best_comb_Pb);
disp([best_npb(best_qb,:)]);
disp([best_npa(best_qa,:)]);
disp([nu_paras(imincombPb(best_comb_q))]);
% disp([best_Pbb(best_qb) best_npb(best_qb,:)]);
% disp([best_Pba(best_qa) best_npa(best_qa,:)]);
% disp([best_comb_Pb(best_comb_q) nu_paras(best_comb_q,imincombPb(best_comb_q))]);

% get the qth row of the minimum Pb among Pb(q)
[~,best_qb] = min(best_Pbb); 
[~,best_qa] = min(best_Pba);
[~,best] = min(min(comb_Pb));

% 
nu_parab(m,:) = best_npb(best_qb(1),:);
nu_paraa(m,:) = best_npa(best_qa(1),:);
nu_parac(m) = nu_paras(best);

fprintf('\nBest nu_paras below & above with %1.0f percent estimated errors\n\n',eta*100);
disp([nu_parab(m,:); nu_paraa(m,:)]);

%% Plot the best exponents overall Pb of q
if do_plot_collapse,  
  lmtd = 1:NGEN;
  fignum = fignum +1;
  plot_data_collapse([alpha best_npb(best_qb,2)],[alpha best_npa(best_qa,2)], ...
    RHO,Jb,Ja,Delta,t,lmtd,fignum);
  set(gca,'xscale','log','yscale','log');
  xlabel('log_{10}( t \Delta^{\nu_{||}} )');
  ylabel('log_{10}( \rho t^\alpha )');
  title(['\alpha = ' num2str(alpha) ... %best alpha
    ' & -\nu_{||} = ' num2str(best_npb(best_qb,2)) ... %best nu_para below
    ' & +\nu_{||} = ' num2str(best_npa(best_qa,2))]);  %best nu_para above
end

% fignum = fignum +1;
% figure(fignum);
% plot_data_collapse([alpha nu_para(m)],[alpha nu_para(m)],RHO,Jb,Ja,Delta,t,lmtd,fignum);
% set(gca,'xscale','log','yscale','log');
% xlabel('log_{10}( t \Delta^{\nu_{||}} )');
% ylabel('log_{10}( \rho t^\alpha )');
% title(['\alpha = ' num2str(alpha) ... %best alpha
%   ' & combined Pb \nu_{||} = ' num2str(nu_para(m))]); %best combined Pb figured nu_para
% 
% fignum = fignum +1;
% figure(fignum);
% plot_data_collapse([alpha 1.3],[alpha 1.3],RHO,Jb,Ja,Delta,t,lmtd,fignum);
% set(gca,'xscale','log','yscale','log');
% xlabel('log_{10}( t \Delta^{\nu_{||}} )');
% ylabel('log_{10}( \rho t^\alpha )');
% title(['\alpha = ' num2str(alpha) ... %best alpha
%   ' & DP \nu_{||} = ' num2str(1.3)]);  %DP nu_para

write_nu_para_beta;