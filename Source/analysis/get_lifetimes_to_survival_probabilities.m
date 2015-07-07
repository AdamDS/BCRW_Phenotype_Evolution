%% get_lifetimes_to_survival_probabilities.m
%
%
this_script = 'get_lifetimes_to_survival_probabilities';
fprintf([this_script '\n']);
% global SIMOPTS;
% limit = SIMOPTS.limit;
i = 0;  %index control parameter(s)
use = []; %for debugging
% colors = ['rbkrbk'];  %plot curve colors
j = 0;  %index index plot curve colors

N = length(overpop)*length(death_max)*length(mutability);
M = length(SIMS);
lifetimes = zeros(N,M);

for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  [base_name,dir_name] = NameAndCD(make_dir,do_cd);
  r = 0;  %index SIMS
  i = i +1; 
  for run = SIMS, 
    r = r +1;
    run_name = int2str(run);  
    new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
    fprintf('\t%s%s\n',base_name,run_name);
    [p,go,error] = try_catch_load([new_dir_name 'population_' base_name run_name],1,1);
    if go,  
      population = p.population;  clear p
      lifetimes(i,r) = length(population);
    end
  end

%% Start main algorithm
    nlt = numel(lifetimes); %num data points
    slt = sort(lifetimes); %sorted lifetimes
    csslt = ([length(slt):(-1):1]); %cumulative sum of sorted lifetimes
    prob = ones(nlt +1,1); %initialize probability
    prob(1) = 1; %initialize t = 0 probability
    prob(2:end) = (csslt./nlt); %get probabilities
    loc = ones(nlt +1,1); %initializes time
    loc(1) = 0; %initialize time start
    loc(2:end) = slt; %get rest of times
%     if do_survival_probability_plot,  
      figure(825),  
      hold on;
      j = j +1;
      plot(loc,prob,[colors(j) '-x']);
%     end
%% End main algorithm
  end %exist lifetimes or population
end %for mu
end %for dm
end %for op
leg = cellstr(num2str(mutability(use)'));
legend(leg);

%% Old code: 
%   for run = SIMS, 
%     run_name = int2str(run);  
%     fprintf('\t%s%s\n',base_name,run_name);
%     [lt,go,error] = try_catch_load(['lifetimes_' base_name(1:(end-1))],1,1);
%     if go,             
%       lifetimes = lt.lifetimes; clear lt
% %%
% % lifetimes = ngen;
% %       use = [use find(mutability==mu)];
% %% Start main algorithm
%       nlt = numel(lifetimes); %num data points
% %   num_bins = ceil(nlt/10);
% %   [freq,bin_loc] = hist(lifetimes,num_bins);
% %   csf = cumsum(freq);
%       slt = sort(lifetimes); %sorted lifetimes
%       csslt = ([length(slt):(-1):1]); %cumulative sum of sorted lifetimes
%       prob = ones(nlt +1,1); %initialize probability
%       prob(1) = 1; %initialize t = 0 probability
%       prob(2:end) = (csslt./nlt); %get probabilities
%       loc = ones(nlt +1,1); %initializes time
%       loc(1) = 0; %initialize time start
%       loc(2:end) = slt; %get rest of times
%       if do_survival_probability_plot,  
%         figure(825),  
%         hold on;
%         j = j +1;
%         plot(loc,prob,[colors(j) '-x']);
%       end
% %% End main algorithm
% %     end %exist lifetimes
% % %     hold off;
% %   end %for SIMS