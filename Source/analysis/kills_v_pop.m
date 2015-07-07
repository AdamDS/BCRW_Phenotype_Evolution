%% All you need is base_name and run_name from either NamingScheme, NameAndCD, or manual input
% This also works with num_clusters replacing population (just do a word Find & Replace: ctrl+H)
%%
this_script = 'kills_v_pop';
fprintf([this_script '\n']);
global SIMOPTS;
zsn = zeros(NGEN,length(SIMS));
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  K1 = zsn;
  K2 = zsn;
  K3 = zsn;
  N = zsn;
  R1 = zsn;
  R2 = zsn;
  R3 = zsn;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  r = 0;
  for run = SIMS, 
    r = r +1;
    run_name = int2str(run);
    exp_name = [base_name run_name];
    go = 1; [p,go] = try_catch_load(['population_' exp_name],go,1);
    if go==1, 
      [k,go] = try_catch_load(['kills_' exp_name],go,1);
      if go==1, 
        population = p.population;  clear p
        ngen = length(population);
        kills = k.kills;  clear k
        K1(1:ngen,r) = kills(:,1);
        K2(1:ngen,r) = kills(:,2);
        K3(1:ngen,r) = kills(:,3);
        N(1:ngen,r) = population';
        R1(1:ngen,r) = N(1:ngen,r)*2;
        R2(1:ngen,r) = R1(1:ngen,r)-kills(1:ngen,1);
        R3(1:ngen,r) = R2(1:ngen,r)-kills(1:ngen,2);
        
        K1dN(1:ngen,r) = K1(1:ngen,r)./N(1:ngen,r);
        if do_kills_v_pop_plots==1, 
          %% overpopulation deaths
          if ~figure(74100000+mu*1000), figure(74100000+mu*1000);  end
            plot(R1(1:ngen,r),K1(1:ngen,r),'x'); hold on;
            xlabel('raw population');
            ylabel('coalescent deaths');
            title(make_title_name(base_name,''));
            figure(2), plot(R1(500:ngen,r),K1dN(500:ngen,r),'x'); hold on;
            xlabel('raw population');
            ylabel('coalescent deaths / raw population');
            title(make_title_name(base_name,''));

%           %% random deaths
%           figure(74200000+mu*1000+run)
%           plot(R2,kills(:,2),'x');
% 
%           %% wanderer deaths
%           figure(74300000+mu*1000+run)
%           plot(R3,kills(:,3),'x');
        end
        if record_kills_v_pop==1, 
          headers = cell(1,7);
          headers(1) = cellstr('Pop');
          headers(2) = cellstr('1st Raw Pop');
          headers(3) = cellstr('Overpop Death');
          headers(4) = cellstr('2nd Raw Pop');
          headers(5) = cellstr('Random Death');
          headers(6) = cellstr('3rd Raw Pop');
          headers(7) = cellstr('Wanderer Death');
          filename = ['kills_v_pop_' base_name(1:(length(base_name)-1))];
          % lbn = length(mutability); rnss = [3 3+lbn-1];
          %[SUCCESS,MESSAGE]=XLSWRITE(FILE,ARRAY,SHEET,RANGE)
          xlswrite(filename,headers(1:7), ...
            run,['A' int2str(1) ':G' int2str(1)]);
          xlswrite(filename,[N R1 kills(:,1) R2 kills(:,2) R3 kills(:,3)], ...
            run,['A' int2str(2) ':G' int2str(2001)]);
        end
      end
    end
  end
  
end
end
end