%% linear_fit.m
% [m,b,sigm,sigb,chi2] = linear_fit(x,y,sigma,do_plot)
% Plots the data x & y along with the line of best fit. The title of the
% plot reports the slope and intercept values for the line of best fit. The
% fit is derived from least squares fitting.
% Inputs:
% x - vector of horizontal variable
% y - vector of vertical variable
% sigma - error of y (can pass [] if not available)
% do_plot - option to create the figure of plotted data (figure 491)
% Outputs:
% m - slope of best fit line
% b - y-intercept of best fit line
% sigm - error in the computed slope, m
% sigb - error in the computed y-intercept, b
% chi2 - Chi-Squared measure
%
% FOR QUICK USE, you can just do [m,b] = linear_fit(x,y,[],0) to get the
% slope and y-intercepts without their errors or the plot.
%
% -ADS
function [m,b,sigm,sigb,chi2] = linear_fit(x,y,sigma,do_plot,fignum),  
% x = [1:1:10];
% y = rand(1)*[0.8:1:9.8];
% sigma = 0.1*rand(1,length(x));
if ~exists(fignum),  fignum = 491+length(x);  end

sizex = size(x);
sizey = size(y);
sizesigma = size(sigma);
if sizesigma(1)==0, sigma = ones(size(x));  end
    
S = sum(sigma.^(-2));
sumx = sum(x.*(sigma.^-2));
sumy = sum(y.*(sigma.^-2));
sumx2 = sum(x.*x.*(sigma.^-2));
sumxy = sum(x.*y.*(sigma.^-2));

b = (sumy*sumx2 -sumx*sumxy) /(S*sumx2 -sumx*sumx);
m = (S*sumxy -sumy*sumx) /(S*sumx2 -sumx*sumx);

Y = (b+m.*x);

%% Calculate errors from p144 in Numerical Methods for Physics 2nd Ed. by Garcia
if sizesigma(1), %Equations 5.13, error bars are constant for all data input
  sigb = sqrt(sumx2 /(S*sumx2 -sumx*sumx)); 
  sigm = sqrt(S /(S*sumx2 -sumx*sumx));
else, %Equation 5.16, if no error bars
  sigm = sqrt((1/(N-2)).*sum((y -Y).^2));
  sigb = sigm;
end

chi2 = sum((y -Y)./sigma).^2;

if do_plot, 
  flip = 0;
  if ~ishold(fignum),  hold on;  flip = 1;  end
  figure(fignum);
  plot(x,y,'bx');
  plot(x,Y,'k');
  legend('data','best linear fit');
  title(['best fit parameters: slope = ' num2str(m) '; intercept = ' num2str(b)]);
  if flip,  hold off; end
end
end