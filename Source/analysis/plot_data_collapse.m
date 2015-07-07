%% plot_data_collapse.m
% [] = plot_data_collapse(expb,expa,RHO,Jb,Ja,Delta,t,limited,fignum)
%
function [] = plot_data_collapse(expb,expa,RHO,Jb,Ja,Delta,t,limited,fignum), 
figure(fignum); hold off;
c = 'mykbrgc';  c = [c c c c];
for i = Jb,  %for each time series
  if length(RHO{i})<max(limited), lim = limited(1):length(RHO{i});
  else, lim = limited;  end
  tau = t(t>lim(1) & t<=lim(end));  %set time span to fit
  x{i} = tau(RHO{i}(tau)~=0).*(abs(Delta(i)).^expb(2)); %control parameter x critical exponent
  y{i} = RHO{i}(tau).*(tau(RHO{i}(tau)~=0).^expb(1));
  plot(x{i},y{i},c(i)); hold on;
end
if exist('Ja','var')>0,  
  for i = Ja,  %for each time series
    if length(RHO{i})<max(limited), lim = limited(1):length(RHO{i});
    else, lim = limited;  end
    tau = t(t>lim(1) & t<=lim(end));  %set time span to fit
    x{i} = tau(RHO{i}(tau)~=0).*(abs(Delta(i)).^expa(2)); %control parameter x critical exponent
    y{i} = RHO{i}(tau).*(tau(RHO{i}(tau)~=0).^expa(1));
    plot(x{i},y{i},c(i)); hold on;
  end
end
end