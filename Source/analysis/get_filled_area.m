global SIMOPTS;
areas = zeros(length(mutability),length(SIMS),NGEN);
m = 0;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  m = m +1;
  hash_length = mu;
  lsx = (((basic_map_size(2)*2) -1)*2) -1;
  lsy = (((basic_map_size(1)*2) -1)*2) -1;
  hlmx = mod(lsx,hash_length);
  hlmy = mod(lsy,hash_length);
  grid_horizontal = [0:hash_length:lsx+hlmx];
  grid_vertical = [0:hash_length:lsy+hlmy];
  for run = SIMS
    run_name = int2str(run);
    load(['population_' base_name run_name]);
    load(['trace_x_' base_name run_name]);
    load(['trace_y_' base_name run_name]);
    
    ngen = length(population);
    begin = ceil(0.9*ngen);
    u = 0;  v = 0;
    for gen = begin:ngen
      u = v +1; v = sum(population(1:gen));
      refs = u:v;
      area = 0;
      for h = 1:length(grid_horizontal)-1
        filled_x = length(find(trace_x(refs)>grid_horizontal(h) & ...
                               trace_x(refs)<grid_horizontal(h+1)));      
        if filled_x>1
          for v = 1:length(grid_vertical)-1

            filled_y = length(find(trace_y(refs)>grid_horizontal(v) & ...
                                   trace_y(refs)<grid_horizontal(v+1)));
            if filled_y>1
              area = area +1;
            end
          end
        end
      end
      area
      areas(m,run,gen) = area;
    end
  end
end
end
end