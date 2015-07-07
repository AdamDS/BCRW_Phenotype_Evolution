%% Determine correlation time exponent: rho*t^alpha = f(t*Delta^nu_para)
these = find(Delta~=0);
t = 10:10^3;  %maximum time span to fit
limited = fTmM(1):fTmM(2);

dnr = 0.001;
nu_perps = 0.50:dnr:1.00;  %expect 0.734 ~ 0.73
Nnr = length(nu_perps);

da = 0.01;
% alphas = 0.30:da:0.60;  %expect 0.451 ~ 0.45

Na = length(alphas);

thisp = zeros(1,2);
Pbb = zeros(Na,Nnr);
Pba = zeros(Na,Nnr);

Q = 6;
best_Pbb = zeros(Q,1);
best_Pba = zeros(Q,1);
best_expb = zeros(Q,2);
best_expa = zeros(Q,2);
darb = zeros(Q,2);
dara = zeros(Q,2);
dnrb = zeros(Q,2);
dnra = zeros(Q,2);

nu_para = nu_parab;
nu_para = nu_paraa;
% alpha = alphab;
PICKS = [lessthan(end-1) greaterthan(2)];
picks = 0;
for i = PICKS,  %for each time series
  picks = picks +1;
  for q = 1:Q,  
    for k = 1:Na,  %alphas along rows
  %     k = 1;
      alpha = alphas(k);
      alpha = alphab(end);
      for j = 1:Nnr,  %nu_paras along columns
        nu_perp = nu_perps(j);
        % construct attempted time series for data-collapse
        x = cell(M,1);
        y = cell(M,1);
        for m = 1:M,  %construct x & y to fit
          L = min([xl(picks) yl(picks)]);
          if length(RHOX{m,picks})<max(limited), 
            lim = 1:length(RHOX{m,picks});
          else, 
            lim = limited;  
          end
          tau = t(t>lim(1) & t<=lim(end));  %set time span to fit
          x{m} = tau(RHOX{m,picks}(tau)~=0)./(L.^(nu_perp./nu_para(m)));
          y{m} = RHOX{m,picks}(tau).*(tau(RHOX{m,picks}(tau)~=0).^alpha);
        end
%         fprintf('k = %1.0f\t\tj = %1.0f\n',k,j);
        if picks==1,  Jb = 1:m;
          Pbb(k,j) = goodness_of_data_collapse0(x,y,Jb,q);
        end
        if picks==2,  Ja = 1:m;
          Pba(k,j) = goodness_of_data_collapse0(x,y,Ja,q);
        end
      end
    end

    [arbestb(q,picks),nrbestb(q,picks)] = find(Pbb==min(min(Pbb)),1);
    best_Pbb(q,picks) = Pbb(arbestb(q,picks),nrbestb(q,picks));
    best_expb(q,:) = [alphas(arbestb(q,picks)) nu_perps(nrbestb(q,picks))];

    [arbesta(q,picks),nrbesta(q,picks)] = find(Pba==min(min(Pba)),1);
    best_Pba(q,picks) = Pba(arbesta(q,picks),nrbesta(q,picks));
    best_expa(q,:) = [alphas(arbesta(q,picks)) nu_perps(nrbesta(q,picks))];

    % Get exponent error
    eta = 0.05;
    darb(q,1) = get_data_collapse_exponent_error( ...
      best_expb(q,1),alphas,eta,Pbb(:,nrbestb(q,picks)),best_Pbb(q,picks),-1);
    darb(q,2) = get_data_collapse_exponent_error( ...
      best_expb(q,1),alphas,eta,Pbb(:,nrbestb(q,picks)),best_Pbb(q,picks),+1);

    dara(q,1) = get_data_collapse_exponent_error( ...
      best_expa(q,1),alphas,eta,Pba(:,nrbesta(q,picks)),best_Pba(q,picks),-1);
    dara(q,2) = get_data_collapse_exponent_error( ...
      best_expa(q,1),alphas,eta,Pba(:,nrbesta(q,picks)),best_Pba(q,picks),+1);


    dnrb(q,1) = get_data_collapse_exponent_error( ...
      best_expb(q,2),nu_perps,eta,Pbb(arbestb(q,picks),:),best_Pbb(q,picks),-1);
    dnrb(q,2) = get_data_collapse_exponent_error( ...
      best_expb(q,2),nu_perps,eta,Pbb(arbestb(q,picks),:),best_Pbb(q,picks),+1);

    dnra(q,1) = get_data_collapse_exponent_error( ...
      best_expa(q,2),nu_perps,eta,Pba(arbesta(q,picks),:),best_Pba(q,picks),-1);
    dnra(q,2) = get_data_collapse_exponent_error( ...
      best_expa(q,2),nu_perps,eta,Pba(arbesta(q,picks),:),best_Pba(q,picks),+1);

  %   % Plot with best exponents
  %   lmtd = 1:NGEN;
  %   fignum = fignum +1;
  %   plot_data_collapse(best_expb(q,:),best_expa(q,:),RHO,Jb,Ja,Delta,t,trans,lmtd,fignum);
    fprintf('q = %1.0f\t\teta = %1.3f\n',q,eta);
    fprintf('\t\talphab = %1.2f\t\tnu_perpb = %1.2f\n',best_expb(q,1),best_expb(q,2));
    fprintf('\t\t -dab = %1.3f\t\t -dnrb = %1.3f\n',darb(q,1),dnrb(q,1));
    fprintf('\t\t +dab = %1.3f\t\t +dnrb = %1.3f\n',darb(q,2),dnrb(q,2));
    fprintf('\t\talphaa = %1.2f\t\tnu_perpa = %1.2f\n',best_expa(q,1),best_expa(q,2));
    fprintf('\t\t -daa = %1.3f\t\t -dnra = %1.3f\n',dara(q,1),dnra(q,1));
    fprintf('\t\t +daa = %1.3f\t\t +dnra = %1.3f\n',dara(q,2),dnra(q,2));
  end
  disp([best_Pbb(:,picks) best_expb best_expa best_Pba(:,picks)]);
  %% Plot the best exponents overall Pb of q
  lmtd = 1:NGEN;
  fignum = fignum +1;
  eb = best_expb(find(best_Pbb(:,picks)==min(best_Pbb(:,picks)),1),:);
  ea = best_expa(find(best_Pba(:,picks)==min(best_Pba(:,picks)),1),:);
  plot_data_collapse(eb,ea,RHOX,Jb,[],Delta,t,lmtd,fignum);
  set(gca,'xscale','log','yscale','log');
  xlabel('log_{10}( t L^{z} )');
  ylabel('log_{10}( \rho t^\alpha )');
  title(['\alphab = ' num2str(best_expb(find(best_Pbb(:,picks)==min(best_Pbb(:,picks)),1),1)) ...
    ' & \nu_{\perp}b = ' num2str(best_expb(find(best_Pbb(:,picks)==min(best_Pbb(:,picks)),1),2)) ...
    '     \alphaa = ' num2str(best_expa(find(best_Pba(:,picks)==min(best_Pba(:,picks)),1),1)) ...
    ' & \nu_{\perp}a = ' num2str(best_expa(find(best_Pba(:,picks)==min(best_Pba(:,picks)),1),2))]);
  
%   fignum = fignum +1;
%   figure(fignum);
%   plot(best_expb(:,1),best_expb(:,2),'r',best_expa(:,1),best_expa(:,2),'b');

  use = best_Pbb(best_Pbb~=0);
  if length(use)>0, [~,best_qb(picks)] = min(use);  
  else, best_qb(picks) = 1; end
  
  use = best_Pba(best_Pba~=0);
  if length(use)>0, [~,best_qa(picks)] = min(use);
  else, best_qa(picks) = 1; end

  alpharb(picks) = best_expb(best_qb(picks),1);
  nu_perpb(picks) = best_expb(best_qb(picks),2);
  alphara(picks) = best_expa(best_qa(picks),1);
  nu_perpa(picks) = best_expa(best_qa(picks),2);

  arminb(picks) = alpharb(picks)-darb(best_qb(picks),1);
  armina(picks) = alphara(picks)-dara(best_qb(picks),2);
  nrminb(picks) = nu_perpb(picks)-dnrb(best_qb(picks),1);
  nrmina(picks) = nu_perpa(picks)-dnra(best_qa(picks),2);

  armaxb(picks) = alpharb(picks)+darb(best_qb(picks),1);
  armaxa(picks) = alphara(picks)+dara(best_qa(picks),2);
  nrmaxb(picks) = nu_perpb(picks)+dnrb(best_qb(picks),1);
  nrmaxa(picks) = nu_perpa(picks)+dnra(best_qa(picks),2);

  disp([arminb(picks) alpharb(picks) armaxb(picks) nrminb(picks) nu_perpb(picks) nrmaxb(picks); ...
    armina(picks) alphara(picks) armaxa(picks) nrmina(picks) nu_perpa(picks) nrmaxa(picks)]);
end


%% Old code: 
% figure(fignum); fignum = fignum +1;
% surf(log10(Pbb));
% xlabel('\alpha');
% ylabel('\nu_{||}');
% zlabel('Pb-');
% 
% figure(fignum); fignum = fignum +1;
% surf(log10(Pba));
% xlabel('\alpha');
% ylabel('\nu_{||}');
% zlabel('Pb+');

% fignum = fignum +1;
% expb = [0.451 1.295]; expa = expb;
% plot_data_collapse(expb,expa,RHO,Jb,Ja,Delta,t,trans,limited,fignum);
% 
% fignum = fignum +1;
% expb = best_expb; expa = expb;
% plot_data_collapse(expb,expa,RHO,Jb,Ja,Delta,t,trans,limited,fignum);
% 
% fignum = fignum +1;
% expb = best_expa; expa = expb;
% plot_data_collapse(expb,expa,RHO,Jb,Ja,Delta,t,trans,limited,fignum);
% 
% fignum = fignum +1;
% expb = [best_expb(1) best_expa(2)]; expa = expb;
% plot_data_collapse(expb,expa,RHO,Jb,Ja,Delta,t,trans,limited,fignum);
% 
% fignum = fignum +1;
% expb = [best_expa(1) best_expb(2)]; expa = expb;
% plot_data_collapse(expb,expa,RHO,Jb,Ja,Delta,t,trans,limited,fignum);

% q = 1;
% eta = 0.1;
% I = lessthan(end);
% Jb = 1:(I-1);
% [best_exp,pm_best_exp,ibest_exp] = goodness_of_data_collapse2( ...
%   alphas,nu_paras,x,y,Jb,q,eta);

%   % Bhatarcharjee & Seno 2001: goodness of data-collapse
%   q = 1;  
%   eta = 0.2;  %20% level of error
%   % below critical
%   I = lessthan(end);%lessthan(end-1);%lessthan(end); %10, 12
%   Jb = 1:(I-1);%lessthan(1:(end-1)); %1:9, 15:13
%   [npb(m,k),asdf,Pbb(k,:),thisp(1,k)] = goodness_of_data_collapse( ...
%     nu_paras,x,y,Jb,q,eta);
%   dnpbm(m,k) = asdf(1);
%   dnpbp(m,k) = asdf(2);
%   % above critical
%   I = greaterthan(1);%greaterthan(2);%greaterthan(1);  
%   Ja = N:(-1):(I+1);
%   if length(Ja)>1, 
%     [npa(m,k),asdf,Pba(k,:),thisp(2,k)] = goodness_of_data_collapse( ...
%       nu_paras,x,y,Ja,q,eta);
%     dnpbm(m,k) = asdf(1);
%     dnpbp(m,k) = asdf(2);
%   else, 
%     npa(m,k) = 0; dnpam(m,k) = 0; dnpap(m,k) = 0; thisp(2,k) = 0;
%   end
%   % npb = nu_para below;  npa = nu_para above
%   betab(m,k) = alpha*npb(m,k);
%     dbetabm(m,k) = alpha*(npb(m,k)-dnpbm(m,k));
%     dbetabm(m,k) = alpha*(npb(m,k)+dnpbp(m,k));
% 
%   betaa(m,k) = alpha*npa(m,k);
%     dbetaam(m,k) = alpha*(npa(m,k)-dnpam(m,k));
%     dbetaap(m,k) = alpha*(npa(m,k)+dnpap(m,k));

%   if dnpbm(m,k)~=0 || abs(dnpbm(m,k))~=Inf, check(1) = dnpbm(m,k);  end
%   if dnpbp(m,k)~=0 || abs(dnpbp(m,k))~=Inf, check(2) = dnpbp(m,k);  end
%   if dnpam(m,k)~=0 || abs(dnpam(m,k))~=Inf, check(3) = dnpam(m,k);  end
%   if dnpap(m,k)~=0 || abs(dnpap(m,k))~=Inf, check(4) = dnpap(m,k);  end
%   
%   [~,best_error(k)] = min(check);
%   if best_error<3,  plotp = thisp(1); else, plotp = thisp(2); end
% end

% figure(9011); 
% surf(Pbb);
% figure(9012);
% surf(Pba);
% figure(9013);
% plot(alphas,npb(m,:),'b',alphas,npa(m,:),'k');
% [min(min(npb)) min(min(npa))]