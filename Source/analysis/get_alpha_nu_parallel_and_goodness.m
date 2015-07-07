%% Determine correlation time exponent: rho*t^alpha = f(t*Delta^nu_para)
these = find(Delta~=0);
t = 1:10^3;%NGEN;  %maximum time span to fit
limited = fTmM(1):fTmM(2);

dnp = 0.01;
nu_paras = 0.85:dnp:1.75;  %expect 1.295 ~ 1.30
Nnp = length(nu_paras);

da = 0.01;
alphas = 0.40:da:0.50;  %expect 0.451 ~ 0.45
% a = -alpha_decay(m);
% alphas = (a*0.67):da:(a*1.33);
Na = length(alphas);

thisp = zeros(1,2);
Pbb = zeros(Na,Nnp);
Pba = zeros(Na,Nnp);

Q = 4;
best_Pbb = zeros(Q,1);
best_Pba = zeros(Q,1);
best_expb = zeros(Q,2);
best_expa = zeros(Q,2);
dab = zeros(Q,2);
daa = zeros(Q,2);
dnpb = zeros(Q,2);
dnpa = zeros(Q,2);

for q = 1:Q,  
  no_data = 0;
  for k = 1:Na,  %alphas along rows
    alpha = alphas(k);
%     alpha = -alpha_decay(m);
    for j = 1:Nnp,  %nu_paras along columns
      nu_para = nu_paras(j);
      % construct attempted time series for data-collapse
      x = cell(N,1);
      y = cell(N,1);
      for i = 1:N,  %for each time series
        if sum(i==these), %no need for collapsing on critical
          if length(RHO{i})<max(limited), 
            lim = limited(1):length(RHO{i});
          else, 
            lim = limited;  
          end
          tau = t(t>lim(1) & t<=lim(end));  %set time span to fit
          x{i} = tau(RHO{i}(tau)~=0).*(abs(Delta(i)).^nu_para); %control parameter x critical exponent
          y{i} = RHO{i}(tau).*(tau(RHO{i}(tau)~=0).^alpha);
        end
      end
  %     fprintf('k = %1.0f\t\tj = %1.0f\n',k,j);
      eta = 0.1;
      Jb = 1:(lessthan(end-offset));
      Pbb(k,j) = goodness_of_data_collapse0(x,y,Jb,q);
      if greaterthan(offset+1)~=greaterthan(end), 
        Ja = (greaterthan(offset+1)):greaterthan(end);
        Pba(k,j) = goodness_of_data_collapse0(x,y,Ja,q);
      else, Pba(k,j) = Inf; no_data = no_data +1; end
    end
  end

%   figure, mesh(Pbb);  title(['Pbb q = ' int2str(q)]);
%   figure, mesh(Pba);  title(['Pba q = ' int2str(q)]);
  [abestb(q),npbestb(q)] = find(Pbb==min(min(Pbb)),1);
  best_Pbb(q) = Pbb(abestb(q),npbestb(q));
  best_expb(q,:) = [alphas(abestb(q)) nu_paras(npbestb(q))];

%   if no_data==(Na+Nnp), 
    [abesta(q),npbesta(q)] = find(Pba==min(min(Pba)),1);
    best_Pba(q,:) = Pba(abesta(q),npbesta(q));
    best_expa(q,:) = [alphas(abesta(q)) nu_paras(npbesta(q))];
%   end
  
  % Get exponent error
  eta = 0.05;
  dab(q,1) = get_data_collapse_exponent_error( ...
    best_expb(q,1),alphas,eta,Pbb(:,npbestb(q)),best_Pbb(q),-1);
  dab(q,2) = get_data_collapse_exponent_error( ...
    best_expb(q,1),alphas,eta,Pbb(:,npbestb(q)),best_Pbb(q),+1);
  
  daa(q,1) = get_data_collapse_exponent_error( ...
    best_expa(q,1),alphas,eta,Pba(:,npbesta(q)),best_Pba(q),-1);
  daa(q,2) = get_data_collapse_exponent_error( ...
    best_expa(q,1),alphas,eta,Pba(:,npbesta(q)),best_Pba(q),+1);
  
  
  dnpb(q,1) = get_data_collapse_exponent_error( ...
    best_expb(q,2),nu_paras,eta,Pbb(abestb(q),:),best_Pbb(q),-1);
  dnpb(q,2) = get_data_collapse_exponent_error( ...
    best_expb(q,2),nu_paras,eta,Pbb(abestb(q),:),best_Pbb(q),+1);
  
  dnpa(q,1) = get_data_collapse_exponent_error( ...
    best_expa(q,2),nu_paras,eta,Pba(abesta(q),:),best_Pba(q),-1);
  dnpa(q,2) = get_data_collapse_exponent_error( ...
    best_expa(q,2),nu_paras,eta,Pba(abesta(q),:),best_Pba(q),+1);
  
%   % Plot with best exponents
%   lmtd = 1:NGEN;
%   fignum = fignum +1;
%   plot_data_collapse(best_expb(q,:),best_expa(q,:),RHO,Jb,Ja,Delta,t,trans,lmtd,fignum);
  fprintf('q = %1.0f\t\t\teta = %1.2f\t\t\tbms = %1.0f\n',q,eta,bms);
  fprintf('\t\talphab = %1.2f\t\tnu_parab = %1.2f\n',best_expb(q,1),best_expb(q,2));
  fprintf('\t\t -dab = %1.3f\t\t -dnpb = %1.3f\n',dab(q,1),dnpb(q,1));
  fprintf('\t\t +dab = %1.3f\t\t +dnpb = %1.3f\n',dab(q,2),dnpb(q,2));
  fprintf('\t\talphaa = %1.2f\t\tnu_paraa = %1.2f\n',best_expa(q,1),best_expa(q,2));
  fprintf('\t\t -daa = %1.3f\t\t -dnpa = %1.3f\n',daa(q,1),dnpa(q,1));
  fprintf('\t\t +daa = %1.3f\t\t +dnpa = %1.3f\n',daa(q,2),dnpa(q,2));
end
disp([best_Pbb best_expb best_expa best_Pba]);

%% Plot the best exponents overall Pb of q
lmtd = 1:NGEN;
fignum = fignum +1;
eb = best_expb(find(best_Pbb==min(best_Pbb)),:);
ea = best_expa(find(best_Pba==min(best_Pba)),:);
plot_data_collapse(eb,ea,RHO,Jb,Ja,Delta,t,lmtd,fignum);
set(gca,'xscale','log','yscale','log');
xlabel('log_{10}( t \Delta^{\nu_{||}} )');
ylabel('log_{10}( \rho t^\alpha )');
title(['-\alpha = ' num2str(best_expb(find(best_Pbb==min(best_Pbb),1),1)) ... %best alpha below
  ' & -\nu_{||} = ' num2str(best_expb(find(best_Pbb==min(best_Pbb),1),2)) ... %best nu_para below
  '     +\alpha = ' num2str(best_expa(find(best_Pba==min(best_Pba),1),1)) ... %best alpha above
  ' & +\nu_{||} = ' num2str(best_expa(find(best_Pba==min(best_Pba),1),2))]);  %best nu_para above

% fignum = fignum +1;
% figure(fignum);
% plot(best_expb(:,1),best_expb(:,2),'r',best_expa(:,1),best_expa(:,2),'b');

[~,best_qb] = min(best_Pbb(best_Pbb~=0));
[~,best_qa] = min(best_Pba(best_Pba~=0));

alphab(m) = best_expb(best_qb,1);
nu_parab(m) = best_expb(best_qb,2);
alphaa(m) = best_expa(best_qa,1);
nu_paraa(m) = best_expa(best_qa,2);

aminb(m) = alphab(m)-dab(best_qb,1);
amina(m) = alphaa(m)-daa(best_qb,2);
npminb(m) = nu_parab(m)-dnpb(best_qb,1);
npmina(m) = nu_paraa(m)-dnpa(best_qa,2);

amaxb(m) = alphab(m)+dab(best_qb,1);
amaxa(m) = alphaa(m)+daa(best_qa,2);
npmaxb(m) = nu_parab(m)+dnpb(best_qb,1);
npmaxa(m) = nu_paraa(m)+dnpa(best_qa,2);

disp([aminb(m) alphab(m) amaxb(m) npminb(m) nu_parab(m) npmaxb(m); ...
  amina(m) alphaa(m) amaxa(m) npmina(m) nu_paraa(m) npmaxa(m)]);

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