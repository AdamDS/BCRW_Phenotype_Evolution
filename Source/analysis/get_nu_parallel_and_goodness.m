%% Determine correlation time exponent: rho*t^alpha = f(t*Delta^nu_para)
these = find(Delta~=0);
dt = 1;
t = 1:dt:NGEN;  %1xNGEN
limited = fT(1):fT(end);
dnp = 0.1;
np = 0.10:dnp:2.00; %expect 1.295 ~ 1.30
nnp = length(np);
thisp = zeros(1,2);
x = cell(N,nnp);
y = cell(N,nnp);
for j = 1:nnp,  %for each nu_para
  nu_para = np(j);
  for i = 1:N,  %for each time series
    if sum(i==these), %no need for collapsing on critical
      tau = t(t>limited(1) & t<=length(RHO{i}));
%       x{i,j} = tau(RHO{i}(tau)~=0).*(abs(Delta(i)).^nu_para); %control parameter x critical exponent
%       y{i,j} = RHO{i}(tau).*(tau(RHO{i}(tau)~=0).^alpha(m));
      x{i,j} = tau(RHO{i}(tau)~=0).*(abs(Delta(i)).^nu_para); %control parameter x critical exponent
      y{i,j} = RHO{i}(tau).*(tau(RHO{i}(tau)~=0).^alpha_decay(m));
    end
  end
end

% alpha = -alpha_decay(greaterthan(1)-1);
% Bhatarcharjee & Seno 2001: goodness of data-collapse
eta = 0.1;  %20% level of error
for q = 1:4,  
  % below critical
  I = lessthan(end);%lessthan(end-1);%lessthan(end); %10, 12
  Jb = 1:(I-1);%lessthan(1:(end-1)); %1:9, 15:13
  [npb(m,q),dnpb(m,q,:),Pb,thisp(1)] = goodness_of_data_collapse(np,x,y,Jb,q,eta);
  % above critical
  I = greaterthan(offset);
  Ja = N:(-1):(I+1);
  if length(Ja)>1, 
    [npa(m,q),dnpa(m,q,:),Pb,thisp(2)] = goodness_of_data_collapse(np,x,y,Ja,q,eta);
  else, 
    npa(m,q) = 0;
    dnpa(m,q,:) = [0 0];
    thisp(2) = 0;
  end
  best_of(q) = find(Pb==min(Pb),1);
  best_Pb(q) = Pb(best_of(q));
end
best = find(best_Pb==min(best_Pb),1);
betab(m) = alpha*npb(m,best);
dbetab(m,:) = alpha*[npb(m,best)-dnpb(m,best,1) npb(m,best)+dnpb(m,best,2)];
betaa(m) = alpha*npa(m,best);
dbetaa(m,:) = alpha*[npa(m,best)-dnpa(m,best,1) npa(m,best)+dnpa(m,best,2)];

% figure(fignum); 
% fignum = fignum +1;
% figure(fignum); hold on;
% xlabel('\nu_{||}','FontSize',14); ylabel('Goodness','FontSize',14);
% set(gca,'yscale','log');
% hold off;

dnpbm = dnpb(m,best,:);
dnpam = dnpa(m,best,:);
[~,best_error] = min([dnpbm(dnpbm~=0) dnpam(dnpam~=0)]);
if best_error<3,  plotp = thisp(1); else, plotp = thisp(2); end
figure(m+6474); 
for i = Jb,  
%   plot(x{i,thisp(1)}(y{i,thisp(1)}~=0),y{i,thisp(1)}(y{i,thisp(1)}~=0),'r');  hold on;
  plot(x{i,plotp}(y{i,plotp}~=0),y{i,plotp}(y{i,plotp}~=0),'r');  hold on;
end
if thisp(2), 
  for i = Ja,  
%     plot(x{i,thisp(2)}(y{i,thisp(2)}~=0),y{i,thisp(2)}(y{i,thisp(2)}~=0),'b');
    plot(x{i,plotp}(y{i,plotp}~=0),y{i,plotp}(y{i,plotp}~=0),'b');
  end
end
set(gca,'xscale','log','yscale','log');
xlabel('t \Delta^{\nu_{||}}','FontSize',14);
ylabel('\rho t^{\alpha}','FontSize',14);
hold off;
% X{m} = x{greaterthan(1),best_error};
% Y{m} = y{greaterthan(1),best_error};