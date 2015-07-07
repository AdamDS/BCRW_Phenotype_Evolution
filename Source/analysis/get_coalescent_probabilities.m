% trans = 500;  stop = NGEN;
% crR1 = [];
% crK1dN = [];
% for run=SIMS, 
%   crR1 = cat_row(crR1,R1(trans:stop,run));
%   crK1dN = cat_row(crK1dN,K1dN(trans:stop,run));
% end
% [m,b,sigm,sigb] = linear_fit(crR1,crK1dN,[],1);
x = [3:10000];
X = ones(length(x),1)*x;
m32 = 5.4135*10^-5; b32 = 0.5577;
m34 = 2.7338*10^-5; b34 = 0.51342;
m36 = 3.0865*10^-5; b36 = 0.4607;
m38 = 3.3191*10^-5; b38 = 0.41467;
m40 = 3.5776*10^-5; b40 = 0.38043;
m42 = 3.7274*10^-5; b42 = 0.35407;

m = [m32 m34 m36 m38 m40 m42];
b = [b32 b34 b36 b38 b40 b42];
B = ones(length(x),1)*b;
y = x'*m+B;
no = find(y>1);
y(no) = 1;
plot(x,y(:,1),'k'); hold on;
plot(x,y(:,2),'r');
plot(x,y(:,3),'y');
plot(x,y(:,4),'g');
plot(x,y(:,5),'b');
plot(x,y(:,6),'m');
mu = ([32:2:42]/100)';
legend(num2str(mu));
xlabel('raw population before coalescent death, 2*N_{gen-1}');  
ylabel('coalescent death / raw population');
title('approximate proportion of raw population killed by coalescent death Bacterial\_Mu\_x\_Uniform\_2\_Flatscape');

figure,
plot(mu,m);
xlabel('\mu');
ylabel('slopes');
title('slopes of coalescent probability function vs raw population Bacterial\_Mu\_x\_Uniform\_2\_Flatscape');

figure,
plot(mu,b);
xlabel('\mu');
ylabel('intercepts');
title('intercepts of coalescent probability function vs raw population Bacterial\_Mu\_x\_Uniform\_2\_Flatscape');

figure,
[slp,intrcpt,sigm,sigb] = linear_fit(mu',b,[],1);
xlabel('\mu');
ylabel('intercepts');