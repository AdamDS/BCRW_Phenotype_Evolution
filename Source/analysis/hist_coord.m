global SIMOPTS;
ocm = 0;
edge(1) = ((((basic_map_size(1)*2) -1)*2) -1) +0.5;
edge(2) = ((((basic_map_size(2)*2) -1)*2) -1) +0.5;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  ocm = ocm +1;
  if mu==plot_hc_mu(ocm)
    for run = SIMS
      if run==plot_hc_run
        run_name = int2str(run);
        load(['population_' base_name run_name]);
        load(['trace_x_' base_name run_name]);
        load(['trace_y_' base_name run_name]);
        if plot_hc_gen>length(population),  plot_hc_gen = population(end);  end
        u = 1 +sum(population(1:plot_hc_gen-1)); v = sum(population(1:plot_hc_gen));
        figure(3000000+1000*mu); 
          hist(trace_x(u:v),[0:0.25:edge(1)]); 
          xlim([0.5 edge(1)]);
        figure(3100000+1000*mu); 
          hist(trace_y(u:v),[0:0.25:edge(2)]); 
          xlim([0.5 edge(2)]);
        figure(3200000+1000*mu); 
          plot([0:0.25:edge(1)],log10(hist(trace_x(u:v),[0:0.25:edge(1)])),'x'); 
          xlim([0.5 edge(1)]);
        figure(3300000+1000*mu); 
          plot([0:0.25:edge(2)],log10(hist(trace_y(u:v),[0:0.25:edge(2)])),'x'); 
          xlim([0.5 edge(2)]);
      end
    end
  end
end
end
end