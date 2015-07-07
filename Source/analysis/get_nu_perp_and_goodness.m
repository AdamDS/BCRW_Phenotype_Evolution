%% Determine correlation length exponent: rho*t^alpha = f(t/L^z)
  dz = 0.1;
%   midz = 0.567./mean([mean(npa) mean(npb)]);  %expect 0.567 ~ 0.57
  z = 0.10:dz:1.50;
  nz = length(z);
  x = cell(m,nz);
  y = cell(m,nz);
  nr = zeros(2,length(npb));
  dnrm = nr;
  dnrp = dnrm;
  nu_per = dnrp;
  nu_perm = nu_per;
  nu_perp = nu_perm;
  nu_pe = nu_perp;
  nu_pem = nu_pe;
  nu_pep = nu_pem;
  for pick = 1:size(RHOX,2),
    for j = 1:nz,  %for each nu_perp
      Z = z(j);
      for i = 1:m,  %for each time series
        L = min([xl(i) yl(i)]);
  %       pick = lessthan(ceil(length(lessthan)*rand));
%         pick = 4;
        tau = t(t>trans & t<=length(RHOX{i,pick}));
        x{i,j} = tau(RHOX{i,pick}(tau)~=0)./(L.^Z);
        y{i,j} = RHOX{i,pick}(tau).*(tau(RHOX{i,pick}(tau)~=0).^alpha(m));
      end
    end

    % Bhatarcharjee & Seno 2001: goodness of data-collapse
    q = 1;  
    eta = 0.2;  %20% level of error
    thisz = zeros(1,2);
    figure(pick+907);
    J = 1:m;
    [zz(pick),dzz(pick,:),best] = goodness_of_data_collapse(z,x,y,TT(PICKS(pick),:),J,q,eta);  
    % above critical
  %   I = greaterthan(1); Ja = N:(-1):(I+1);
  %   if length(Ja)>1, 
  %     [za(m),dnpa(m,:),thisp(2)] = goodness_of_data_collapse(z,X,Y,I,Ja,q,eta);
  %   else, 
  %     za(m) = 0;
  %     dnza(m,:) = [0 0];
  %     dis(2) = 0;
  %   end
    % z = nu_perp/nu_para
    if pick==1, 
      nr(pick,:) = zz(pick)*npb;%*npb';
      dnrm(pick,:) = (zz(pick)-dzz(pick,1)).*(npb-dnpb(pick,1));
      dnrp(pick,:) = (zz(pick)+dzz(pick,2)).*(npb+dnpb(pick,2));
    elseif pick==2, 
      nr(pick,:) = zz(pick)*npa;%*npa'
      dnrm(pick,:) = (zz(pick)-dzz(pick,1)).*(npa-dnpa(pick,1));
      dnrp(pick,:) = (zz(pick)+dzz(pick,2)).*(npa+dnpa(pick,2));
    end
    nu_per(pick,:) = zz(pick)*nu_par;
    nu_perm(pick,:) = (zz(pick)-dzz(pick,1)).*(nu_par-enu_par);
    nu_perp(pick,:) = (zz(pick)+dzz(pick,2)).*(nu_par+enu_par);
    nu_pe(pick,:) = zz(pick)*nu_pa;
    nu_pem(pick,:) = (zz(pick)-dzz(pick,1)).*(nu_pa-enu_pa);
    nu_pep(pick,:) = (zz(pick)+dzz(pick,2)).*(nu_pa+enu_pa);
      
    figure(pick+907);  
    xlabel('z','FontSize',14); ylabel('Goodness','FontSize',14);
    set(gca,'yscale','log');
    hold off;

    figure(pick+6414); 
    for i = J,  
      plot(x{i,best}(y{i,best}~=0),y{i,best}(y{i,best}~=0),'r');  hold on;
    end
    set(gca,'xscale','log','yscale','log');
    xlabel('t/L^{z}','FontSize',14);
    ylabel('\rho t^{\alpha}','FontSize',14);
    hold off;
  end