%% plot_order_v_control.m
% plot_order_v_control(control,order,stdo,fn,xname,yname,c,base_name)
% control = vector of control parameter
% order = vector of order parameter
% stdo = vector of std of order parameter
% fn = figure number
% xname = xlabel
% yname = ylabel
% c = color character
% base_name = as it sounds
%
function [] = plot_order_v_control(control,order,stdo,fn,xname,yname,c,base_name), 
nz = find(order);
%% Plot populations by control parameter
figure(fn);
tn = make_title_name(generalize_base_name(base_name),'');
% if length(nz),  
  errorbar(control(nz),order(nz),stdo(nz),c);  title(tn,'FontSize',16);
  xlabel(xname,'FontSize',14);  ylabel(yname,'FontSize',14);
  if length(control)>1,  xlim([min(control) max(control)]);  end
  if stdo, 
    figure(fn+1);
    tn = [make_title_name(generalize_base_name(base_name),'')];
    plot(control(nz),stdo(nz),[c '-x']); title(tn);
    xlabel(xname);  ylabel(['\sigma_{' yname '}']);
    if length(control)>1,  xlim([min(control) max(control)]);  end
  end
% end
end
 %ylim([0 10000]); 
%   set(gca,'XTick',[min(mutability):((max(mutability)-min(mutability))/5):max(mutability)],...
%     'Box','on','FontSize',14,'TickDir','in','TickLength',[0.015,0.025]);
%   p = min(mutability):((max(mutability)-min(mutability))/5):max(mutability);  
%   P = {num2str(p(}
%   set(gca,'XTick',[p],'XTickLabel',{num2str(p)});