%% get_nu_parallel.m *******************************************************
% Uses data collapse.
global SIMOPTS;
m = 0; c = 'ybrmgk'; c = [c c c c c]; fignum = 5150;

num_sims = zeros(length(mutability),1);

win_pops = zeros(length(mutability),NGEN-1);

N = numel(mutability)*numel(death_max)*numel(overpop);
M = numel(basic_map_sizes);

mlt = zeros(N,M);
alphas = zeros(N,M);
sigas = zeros(N,M);
fracsur = zeros(N,M);

alpha_crit = zeros(M,1);
RHOX = cell(M,1);

sim_exists = zeros(N,length(SIMS));

muc = 0.33;

for bms = basic_map_sizes,  
  m = m +1;
  SIMOPTS.basic_map_size = [bms bms]; 
  if bms==20, SIMrange = '1_60';  else, SIMrange = [int2str(SIMS(1)) '_' int2str(SIMS(end))]; end
  maxpop = get_max_population();
  [A(m),xl(m),yl(m)] = landscape_measures(SIMOPTS.basic_map_size);

  fprintf('Landscape sizes: xl=%1.0f\t\tyl=%1.0f\t\tA=%1.0f\n',xl(m),yl(m),A(m));

  get_population_decay_data;  %get population data

  %% Determine critical point
%   muc = 0.33;
  Delta = mutability-muc;
  lessthan = find(Delta<-0.0001);
  [~,iscrit] = min(abs(Delta));
  greaterthan = find(Delta>0.0001);
  
  alpha_crit(m) = -alpha_decay(iscrit,m);

%   RHOX{m} = RHO{iscrit};
  RHOX{m} = RHO{lessthan(end)};

  get_nu_parallel_and_goodness_with_alpha;
  
end %bms

get_nu_perpendicular_and_goodness_with_alpha;
  
%% Old code: 
%   %% Determine critical point
% %   muc(m) = input('What is the critical mutability (want alpha ~ 0.451)?'); %determine critical mu
% %   muc(m) = mutability(find(fracsur(:,m)>0,'last'));% -(dmu/2);
% %   Delta = (mutability-muc(m));
%   muc(m) = 0.33;
%   Delta = mutability-muc(m);
%   lessthan = find(Delta<-0.0001);
%   greaterthan = find(Delta>0.0001);
% %   [mD,imD] = min(abs(Delta(lessthan)));
% %   alpha(m) = -alphas(lessthan(end),m);
% %   alpha_decay(m) = -mean(alphas([lessthan(end) greaterthan(1)],m));
%   btwn = [lessthan(end) greaterthan(1)];
% %   if min(Delta)~=0, 
% %     alpha_decay_interp(:,m) = interp1(btwn,alpha_decay(btwn,m),[btwn(1) muc(m) btwn(end)]);
% %     alpha_crit(m) = -alpha_decay_interp(2,m);
% %   else, 
%     alpha_crit(m) = -alpha_decay(i,lessthan(end)+1);
% %   end
%   
% %   PICKS = [lessthan(end):greaterthan(1)];
% %   RHOX{m,1} = RHO{lessthan(end-offset)};
% %   RHOX{m,2} = RHO{find(min(Delta),1)};
% %   RHOX{m,3} = RHO{greaterthan(offset+1)};
%   RHOX{m} = RHO{find(min(Delta),1)};
%   
% %   get_nu_parallel_and_goodness;
% %   get_alpha_nu_parallel_and_goodness;
%   get_nu_parallel_and_goodness_with_alpha;
%   
% %   [nu_par(m),~,enu_par(m),~] = linear_fit( ...
% %     log10(abs(Delta(lessthan))),log10(T(lessthan,2)'),0.1*ones(size(Delta(lessthan))),1);
% %   [nu_pa(m),~,enu_pa(m),~] = linear_fit( ...
% %     log10(abs(Delta(lessthan))),log10(mlt(lessthan,m))',0.1*ones(size(Delta(lessthan))),1);
% %   bet(m) = alpha(m)*abs(nu_par(m));
% %   dbet(m,1) = alpha(m)*abs(nu_par(m)-enu_par(m));
% %   dbet(m,2) = alpha(m)*abs(nu_par(m)+enu_par(m));
% %   be(m) = alpha(m)*abs(nu_pa(m));
% %   dbe(m,1) = alpha(m)*abs(nu_pa(m)-enu_pa(m));
% %   dbe(m,2) = alpha(m)*abs(nu_pa(m)+enu_pa(m));
% %   
% %   
% %   fprintf(['nu_para-(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n' ...
% %     'nu_para+(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n' ...
% %     'nu_paraR(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n' ...
% %     'nu_paraMR(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n' ...
% %     'beta-(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n' ...
% %     'beta+(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n' ...
% %     'betaR(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n' ...
% %     'betaMR(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n\n\n'], ...
% %     npb(m)-dnpb(m,1),npb(m),npb(m)+dnpb(m,2), ...
% %     npa(m)-dnpa(m,1),npa(m),npa(m)+dnpa(m,2), ...
% %     nu_par(m)-enu_par(m),nu_par(m),nu_par(m)+enu_par(m), ...
% %     nu_pa(m)-enu_pa(m),nu_pa(m),nu_pa(m)+enu_pa(m), ...
% %     betab(m)-dbetab(m,1),betab(m),betab(m)+dbetab(m,2), ...
% %     betaa(m)-dbetaa(m,1),betaa(m),betaa(m)+dbetaa(m,2), ...
% %     dbet(m,1),bet(m),dbet(m,2), ...
% %     dbe(m,1),be(m),dbe(m,2));
% end %bms
% 
% get_alpha_nu_perpendicular_and_goodness;
% 
% fignum = fignum +1;
% figure(fignum);
% plot(1:M,alphab,'-xb',1:M,alphaa,'-xr');
% 
% fignum = fignum +1;
% figure(fignum);
% plot(1:M,nu_parab,'-xb',1:M,nu_paraa,'-xr');



% get_nu_perp_and_goodness;
% 
% fprintf(['nu_perp-(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n' ...
%   'nu_perp+(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n' ...
%   'nu_perpR(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n' ...
%   'nu_perpMR(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n' ...
%   'z-(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n' ...
%   'z+(-=+)\t\t\t%1.4f\t\t%1.4f\t\t%1.4f\n\n\n'], ...
%   dnrm(1,1),nr(1),dnrp(1,2), ...
%   dnrm(2,1),nr(2),dnrp(2,2), ...
%   nu_perm(2,end),nu_per(2,end),nu_perp(2,end), ...
%   nu_pem(2,end),nu_pe(2,end),nu_pep(2,end), ...
%   zz(1)-dzz(1,1),zz(1),zz(1)+dzz(1,2), ...
%   zz(2)-dzz(2,1),zz(2),zz(2)+dzz(2,2));

% toc;
% 
% % fprintf(['betab=%1.3f\t\tnu_parab=%1.3f\t\tnu_perpb=%1.3f\n'...
% %   'betaa=%1.3f\t\tnu_paraa=%1.3f\t\tnu_perpa=%1.3f\n'],...
% %   betab(m),npb(m),nrb(m),betaa(m),npa(m),nra(m));
% 
% disp([' landscape         ' num2str(xl)]);
% disp(['  critical         ' num2str(muc)]);
% disp(['alpha (0.451)      ' num2str(alpha)]);
% disp(['beta- (0.584)      ' num2str(betab)]);
% disp(['beta+ (0.584)      ' num2str(betaa)]);
% disp(['betaR (0.584)      ' num2str(bet)]);
% disp(['betaMR (0.584)     ' num2str(be)]);
% disp(['nu_para- (1.295)   ' num2str(npb)]);
% disp(['nu_para+ (1.295)   ' num2str(npa)]);
% disp(['nu_paraR (1.295)   ' num2str(abs(nu_par))]);
% disp(['nu_paraMR (1.295)  ' num2str(abs(nu_pa))]);
% disp(['nu_perp- (0.7343)  ' num2str(nr(1,:))]);
% disp(['nu_perp+ (0.7343)  ' num2str(nr(2,:))]);
% disp(['nu_perpR- (0.7343)  ' num2str(abs(nu_per(1,:)))]);
% disp(['nu_perpR+ (0.7343)  ' num2str(abs(nu_per(2,:)))]);
% disp(['nu_perpMR- (0.7343) ' num2str(abs(nu_pe(1,:)))]);
% disp(['nu_perpMR+ (0.7343) ' num2str(abs(nu_pe(2,:)))]);
% disp(['     z (0.567)     ' num2str(zz)]);

% disp([xl; muc; alpha; betab; betaa; npb; npa; nrb; nra]);
% disp(zz);

% figure, plot(mu,abs(Delta).^-nr,'r-x')
% hold on;
% plot([muc-0.01 muc+0.01],[77 77],'k')
% plot([muc-0.01 muc+0.01],[45 45],'k')
% plot([muc-0.01 muc+0.01],[37 37],'k')
% plot([muc-0.01 muc+0.01],[29 29],'k')
% plot([muc-0.01 muc+0.01],[21 21],'k')

% Write to xls
% col = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'; 
% for i=1:11, 
%   xlswrite(['RHO_' int2str(use(1)) '_' int2str(diff(use(1:2))) '_' int2str(use(end)) ...
%     '_mc_335_delta_462_nupar_1_' int2str(mod(nu_parallel,1)*100) '_mu_33_001_34'],...
%     RHO{i}',...
%     [col(i) int2str(2) ':' col(i) int2str(length(RHO{i}))]); 
% end


% Old-old code: 
% 
% %% Alternative scaling
% % Tmax = 10^4;
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_33_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 1; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_331_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 2; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_332_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 3; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_333_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 4; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_334_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 5; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_335_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 6; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_336_Uniform_2_Flatscape_1000000_NGEN_1_100.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 7; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_337_Uniform_2_Flatscape_1000000_NGEN_1_100_300_359.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 8; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_338_Uniform_2_Flatscape_1000000_NGEN_300_359.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 9; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_339_Uniform_2_Flatscape_1000000_NGEN_300_359.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 10; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_34_Uniform_2_Flatscape_1000000_NGEN_300_365.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 11; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_341_Uniform_2_Flatscape_1000000_NGEN_305_319_325_329_355_359.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 12; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_342_Uniform_2_Flatscape_1000000_NGEN_305_319.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 13; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_343_Uniform_2_Flatscape_1000000_NGEN_310_314.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 14; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % mu = 0.330:0.001:0.343;
% % mucrit = 0.335; eps = abs(mu(mu<mucrit)-mucrit);
% % figure(1), plot(eps,T(1:length(eps)),'-x');  set(gca,'yscale','log','xscale','log');
% % figure(2), plot(eps,varT(1:length(eps)),'-x'); set(gca,'yscale','log','xscale','log');
% 
% 
% %%
% if only_lt, 
%     
%   else, 
%     %build the simulation names for which you have chosen
%     gen_pops = zeros(length(SIMS),NGEN);
%     i = i +1;
%     r = 0;
%     for run = SIMS, 
%       r = r +1;
%       new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
%       run_name = int2str(run);
%       exp_name = [base_name run_name];
%       pop_name = [new_dir_name 'population_' exp_name];
%       [p,go] = try_catch_load(pop_name,1,1);
%       if go, 
%         num_sims(i) = num_sims(i) +1;
%         sim_exists(i,r) = 1; %track which simulation data is used
%         population = p.population;  clear p
%         g = find(population>=limit);
%         actual_ngen = max(g);
%         ngens(i,r) = actual_ngen;
%         avg_pop(i,run) = mean(population(g));
%         gen_pops(run,1:length(g)) = population(g);
%       end
%     end
%   end
%
%
%   use = find(sim_exists(i,:)); these = 1:length(use);
%   gen_pops(~use,:) = [];
%     
%   [rgood(i),cgood] = size(gen_pops);
%   
%   AVG_POPS(i) = mean(avg_pop(i,:));
%   STD_POPS(i) = std(avg_pop(i,:));
% 
%   for gen=1:max(ngens(i,:)),  
%     alive = find(gen_pops(:,gen));
%     AVG_GEN_POPS(i,gen) = mean(gen_pops(alive,gen),1);
%     STD_GEN_POPS(i,gen) = std(gen_pops(alive,gen));
%   end
%   
% %   for gen = 1:(NGEN-window),  
% %     win_pops(run,gen) = mean(AVG_GEN_POPS(i,gen:(gen+window)),2);
% %   end
% %   for gen = 1:NGEN, 
% %     use = find(gen_pops(:,gen));
% %     if isempty(use),  gen = NGEN-1; 
% %     else, sur_pops(i,gen) = mean(gen_pops(use,gen));  end
% %   end
% 
% 
% %% LOAD STUFF
% % %   color = c(i); 
% % close all;
% % clear all;
% % clc
% % % C = 'kkkkkk';
% % C = 'brykg'; 
% % c = [C C C];
% % mucrit = 0.335;
% % delta = 0.462; %0.462 +- 0.033 from mu=0.335 for 100sims
% % nu_parallel = 1.19; %0.335 1.05
% % A = 37544;%45*45;%/37544;
% % t = 1:10^6;
% % use = 200:200:100000;
% % mutability = 0.33:0.001:0.34;
% % RHO = cell(numel(mutability),1);
% % T = t(use);
% % 
% % i = 1;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_33_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 2;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_331_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 3;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_332_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 4;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_333_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 5;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_334_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 6;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_335_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 7;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_336_Uniform_2_Flatscape_1000000_NGEN_1_100.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 8;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_337_Uniform_2_Flatscape_1000000_NGEN_1_100_300_359.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 9;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_338_Uniform_2_Flatscape_1000000_NGEN_300_359.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 10;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_339_Uniform_2_Flatscape_1000000_NGEN_300_359.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 11;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_34_Uniform_2_Flatscape_1000000_NGEN_300_365.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% %%