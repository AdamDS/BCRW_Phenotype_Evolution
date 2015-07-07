%% plot_finite_size.m
% [] = plot_finite_size(exp,RHO,J,L,t,limited,fignum)
%
function [] = plot_finite_size(exp,RHO,J,L,t,limited,fignum), 
figure(fignum); hold off;
c = 'mykbrgc';  c = [c c c c];
for i = J,  %for each time series
  if length(RHO{i})<max(limited), lim = limited(1):length(RHO{i});
  else, lim = limited;  end
  tau = t(t>lim(1) & t<=lim(end));  %set time span to fit
  x{i} = tau(RHO{i}(tau)~=0)./(L(i).^exp(2)); %control parameter x critical exponent
  y{i} = RHO{i}(tau).*(tau(RHO{i}(tau)~=0).^exp(1));
  plot(x{i},y{i},c(i)); hold on;
end
end