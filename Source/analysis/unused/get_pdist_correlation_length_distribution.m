%% get_characteristic_lengths.m *******************************************************
global SIMOPTS;
this_script = 'get_characteristic_lengths';
tic;
%defaults
bms = SIMOPTS.basic_map_size;
X = (((bms(1)*2)-1)*2)-1; Y = (((bms(2)*2)-1)*2)-1;
if N_samples==[], N_samples = 10;  end
if epsilon==[], epsilon = SIMOPTS.overpop;  end
if de==[],  de = 1; end
if EPSILON==[], EPSILON = max([X Y])*4/5;  end
transient = 200;
%% initialize outputs
if do_characteristic_diameters==1,  CD = [];  CD_finite = [];  end
if do_correlation_lengths==1, 
  Neighborhoods = ([epsilon:de:EPSILON]);
  avg_N = zeros(length(Neighborhoods),length(SIMS));
  std_N = zeros(length(Neighborhoods),length(SIMS));
  AVG_N = zeros(length(Neighborhoods),length(mutability));
  STD_N = zeros(length(Neighborhoods),length(mutability));
end
%%
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  i = find(mu==mutability);
  for run = SIMS
    run_name = int2str(run);
    fprintf([this_script ' for ' base_name run_name '\n']);
    exp_name = [base_name run_name];
    if do_characteristic_diameters==1, 
      cd_name = ['cluster_diameters_' exp_name];
      [cd,go] = try_catch_load(cd_name,1,1);
      if go==1
        cluster_diameters = cd.cluster_diameters;  clear cd
        CD_finite = [CD_finite; cluster_diameters];
%         CD_finite = [CD_finite; cluster_diameters(cluster_diameters~=Inf)];
      end
    end
    
    if do_correlation_lengths==1, 
    go = 1; [p,go] = try_catch_load(['population_' exp_name],go,1);
    if go==1, [x,go] = try_catch_load(['trace_x_' exp_name],go,1);
    if go==1, [y,go] = try_catch_load(['trace_y_' exp_name],go,1);
    if go==1, 
      population = p.population;  clear p
      trace_x = x.trace_x;  clear x
      trace_y = y.trace_y;  clear y
      ngen = length(population);
      avg_n = zeros(length([epsilon:de:EPSILON]),ngen);
      u = 0;  v = 0;
      for gen=transient:ngen,
        script_gen_update(this_script,gen,base_name,run_name);
        u = v +1; v = sum(population(1:gen));
        n = zeros(length(Neighborhoods),N_samples);
        for j=1:N_samples, 
          x = X*rand(1) +0.5;  y = Y*rand(1) +0.5;
          e = 0;
          distances = sqrt((x-trace_x(u:v)).^2 +(y-trace_y(u:v)).^2);
          ri = 0;
          for radius=Neighborhoods, 
            ri = ri +1;
            n(ri,j) = length(find(distances<=radius));
          end
        end
        
%         avg_n(:,gen) = mean(n,2);
      end
%       avg_N(:,run) = mean(avg_n(:,1:ngen),2);
%       std_N(:,run) = std(avg_n')';
    end
    end
    end
%     AVG_N(:,i) = mean(avg_N,2);
%     STD_N(:,i) = std(avg_N')';
    end
    
  end
end
end
end
toc;

if do_characteristic_diameters==1, 
  [asdf,qwer] = hist(CD_finite,100);
  [m,b,sm,sb] = linear_fit(log10(Neighborhoods(25:75)'),log10((25:75)),log10((25:75)),1);
end
if do_correlation_lengths==1, 
  [mcl,bcl,smcl,sbcl] = linear_fit(log10(Neighborhoods(5:15)'),log10(AVG_N(5:15,end)),log10(STD_N(5:15)),1);
end
% xlswrite(['characteristic_lengths_' generalize_base_name(base_name) '__5sim'],...
%          [Neighborhoods, AVG_N]);