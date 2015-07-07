global SIMOPTS;
close all;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  for run = SIMS
    run_name = int2str(run);
    go = 1;
    [nd,go1,error] = try_catch_load(['num_descendants_' base_name run_name],go);
    if go1==0, 
      [par,go,error] = try_catch_load(['parents_' base_name run_name],go);
      if go==1, 
        [pop,go,error] = try_catch_load(['population_' base_name run_name],go);
        if go==1, 
          population = pop.population;  clear pop
          parents = par.parents;  clear par

          ngen = length(find(population>=limit));
          ipop = population(1);
          num_descendants = zeros(ipop,ngen);
          for indiv = 1:ipop
            working_on_indiv = indiv
            these_orgs = indiv;
            pars = these_orgs;
            g = 0;
            pu = 0; pv = 0;
            npu = 0;  npv = 0;
            num_pars = zeros(1,ngen);
            num_pars(1) = 1;
            for gen = 2:ngen
              g = g +1;
              pu = pv +1; pv = sum(population(2:gen));
              RO = parents(pu:pv,1);
      %         if limit>2, NN = parents(pu:pv,2);  end
              npu = npv +1; npv = length(pars);
              these_orgs = pars(npu:npv);
              num_pars(g) = length(these_orgs);
              for i = 1:num_pars(g)
                this_org = these_orgs(i);
                pars = [pars; find(RO==this_org)]; 
              end
            end
            num_descendants(indiv,:) = num_pars;
            rat(indiv,:) = num_pars./population(2:end);
          end
          save(make_data_name('num_descendants',base_name,run_name,0),...
                              'num_descendants');
        end
      end
    else
%       load(make_data_name('num_descendants',base_name,run_name,0));
      num_descendants = nd.num_descendants; clear nd
    end
    if do_descendants_mesh==1
      figure(mu*1000 +run); mesh(num_descendants);  title(make_title_name(base_name,run_name));
      xlabel('generation'); ylabel('ancestor'); zlabel('num\_descendants');
    end
  end
end
end
end