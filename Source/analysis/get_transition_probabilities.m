%% get_transition_probabilities.m
%
global SIMOPTS;
% binsize = 50;
limit = SIMOPTS.limit;
n = 0;
N = numel(mutability)*numel(death_max)*numel(overpop);
tic;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  n = n +1;
  r = 0;
  ngens = zeros(numel(SIMS),1);
  mps = zeros(numel(SIMS),1);
  for run = SIMS, 
    r = r +1;
    [p,go] = try_catch_load(['population_' base_name int2str(run)],1,0);
    if go,  population = p.population;  clear p,  
    ngens(r) = numel(population);
    mps(r) = max(population);
    end
  end
  mp = max(mps);
  mngen = max(ngens);
  boost = mod(mp,binsize);
  bincen = ceil(binsize/2):binsize:(mp+boost+ceil(binsize/2));
  %lower bin values (-0.1 bc HeavisideMatrix checks strictly greater than)
  lowerbound = bincen -ceil(binsize/2) -0.1; 
  %upper bin values (-0.1 bc HeavisideMatrix checks strictly greater than)
  upperbound = bincen +floor(binsize/2) -0.1;
  w = zeros(numel(bincen));
  transition_probabilities = zeros(numel(bincen));
  r = 0;
  for run = SIMS, 
    r = r +1;
    run_name = int2str(run);
    [p,go] = try_catch_load(['population_' base_name run_name],1,1);
    if go==1, 
    population = p.population;  clear p
    %% Debug level start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     population = [300 400 200 320 410 360 600 500 230 400 340 440];
%     mps = max(population);
%     binsize = 50;
%     mp = max(mps);
%     boost = mod(mp,binsize);
%     bincen = ceil(binsize/2):binsize:(mp+boost+ceil(binsize/2));
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    this_p_in_bin = zeros(numel(bincen),mngen);
    norm_const = zeros(numel(bincen),1);
    for i = 1:numel(bincen), 
      lower_p = HeavisideMatrix(population,lowerbound(i)); %find values greater than lowerbound(i)
      upper_p = HeavisideMatrix(population,upperbound(i)); %find values less than upperbound(i)
      this_p_in_bin(i,1:ngens(r)) = mod(lower_p +upper_p,2); %common values are in bin
    end
    w = w +this_p_in_bin(:,(1:(ngens(r)-1)))*this_p_in_bin(:,(2:ngens(r)))';
    %% Debug level end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %             25  75  125 175 225 275 325 375 425 475 525 575 625
%     expect_w =  [[0   0   0   0   0   0   0   0   0   0   0   0   0];... %25
%                  [0   0   0   0   0   0   0   0   0   0   0   0   0];... %75
%                  [0   0   0   0   0   0   0   0   0   0   0   0   0];... %125
%                  [0   0   0   0   0   0   0   0   0   0   0   0   0];... %175
%                  [0   0   0   0   0   0   0.5 0   0.5 0   0   0   0];... %225
%                  [0   0   0   0   0   0   0   0   0   0   0   0   0];... %275
%                  [0   0   0   0   0   0   0   0   1   0   0   0   0];... %325
%                  [0   0   0   0   0   0   0   0   0   0   0   0   1];... %375
%                  [0   0   0   0   0.3 0   0.3 0.3 0   0   0   0   0];... %425
%                  [0   0   0   0   0   0   0   0   0   0   0   0   0];... %475
%                  [0   0   0   0   1   0   0   0   0   0   0   0   0];... %525
%                  [0   0   0   0   0   0   0   0   0   0   0   0   0];... %575
%                  [0   0   0   0   0   0   0   0   0   0   1   0   0]];   %625
%     help_bin = [bincen', this_p_in_bin; [0 population]];
%     help_w = [bincen', w; [0 bincen]];
%     help_ew = [bincen', expect_w; [0 bincen]];
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end %population
  end %for sims
  visit_probabilities = w/sum(sum(w));
  norm_consts = sum(w,2); %normalization constant for each binning row
  inorm = find(norm_consts);
  for j = inorm',  transition_probabilities(j,:) = w(j,:)/norm_consts(j); end
  skewed(n) = sum(sum(triu(visit_probabilities)))/sum(sum(triu(visit_probabilities')));
end %for mu
end %for dm
end %for op

if do_transition_mesh,  
  figure(63900000 +10000*mu +10*reproduction +limit);
  surf(transition_probabilities);
%   ae = [1 10:10:numel(bincen)]; k = 0;
%   label = char(numel(ae),1);
%   for lbl = ae,  
%     k = k +1;
%     label(k) = int2str(bincen(lbl));
%   end
%   set(gca,'XTickLabel',label,'YTickLabel',label);
  xlabel('N_{gen}'); ylabel('N_{gen+1}'); zlabel('transition probability');
  title([make_title_name(base_name,'') ' binsize = ' int2str(binsize)]);
end

if do_transition_ratio_plot,  
  figure(38900 +10*reproduction +limit);
  plot(mutability,skewed,'-x');
end