%% check_duplication.m
%
global SIMOPTS;
this_script = 'check_duplication';
n = 0;
N = numel(mutability)*numel(death_max)*numel(overpop);
R = numel(SIMS);
ngens = zeros(N,R);
duplicated = 0;
DUPLICATES = [];
tic;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  [base_name,dir_name] = NameAndCD(0,do_cd);
  populations = zeros(R,NGEN);
  go = ones(R,1);
  n = n +1;
  r = 0;
  fprintf([this_script ' ' base_name(1:(end-1)) '\n']);
  for run = SIMS, 
    r = r +1;
    run_name = int2str(run);
    [p,go(r)] = try_catch_load(['population_' base_name run_name],1,0);
    if go(r),  population = p.population;  clear p,  
    ngen = length(population);
    ngens(n,r) = ngen;
    populations(r,1:ngen) = population;
    end %population
  end %for SIMS
  un = numel(unique(ngens(n,:)));
%   old = duplicated;
  ds = [];  di = 0;
  if un~=R, 
    for i = 1:R,  
      for j = i+1:R,  
        if i~=j && go(i) && go(j),  
          p1 = populations(i,:);
          p2 = populations(j,:);
          dp = length(find(p1-p2));
          if ~dp, %if dp = 0, then p1==p2
          if ~length(find(j==ds)), %if run pair not yet deemed duplicated
            di = di +1;
            if di==1, ds(di) = i; di = di +1; end %include initial run duplicated
            ds(di) = j; %append runs duplicated
            duplicated = duplicated +1;
            if numel(mutability)>numel(death_max),  
              DUPLICATES = cat_row(DUPLICATES,[mu SIMS(i) SIMS(j)]);
            else, 
              DUPLICATES = cat_row(DUPLICATES,[dm SIMS(i) SIMS(j)]);
            end
            fprintf(['DUPLICATES: ' base_name(1:(end-1)) ' runs %1.0f and %1.0f \n'],...
                    SIMS(i),SIMS(j));
          end
          end
        end
      end
    end
  end
%   new = duplicated;
%   DUPLICATES(old+1:new,:)
end %for mu
end %for dm
end %for op
fprintf(['\n' 'Number of duplications: %1.0f \n'],duplicated);
if duplicated,  
  fprintf('Type DUPLICATES to see which sets were duplicated. \n');  
  fprintf('DUPLICATES is organized in rows of (value run run). \n');
end