%% get_cluster_density_dimension.m
%
%
this_script = 'get_cluster_density_dimensions';
fprintf([this_script '\n']);
global SIMOPTS;
limit = SIMOPTS.limit;
samples = 1000;
init_z = zeros(length(overpop)*length(death_max)*length(mutability),length(SIMS));
M = init_z; B = init_z; EM = init_z;  EB = init_z;
SM = init_z; SB = init_z; SEM = init_z;  SEB = init_z;
LM = init_z; LB = init_z; LEM = init_z;  LEB = init_z;
D = zeros([size(init_z) samples]);
i = 0;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  r = 0;  i = i +1;
  for run = SIMS, 
    run_name = int2str(run);
    fprintf([this_script ' for ' base_name run_name '\n']);
%     [cdm,go,error] = try_catch_load(['cluster_diameters_' base_name run_name],1,1);
    [gr,go,error] = try_catch_load(['gyration_radii_' base_name run_name],1,1);
    if go, [onc,go,error] = try_catch_load(['orgsnclusters_' base_name run_name],go,1);
    if go,             
%       cluster_diameters = cdm.cluster_diameters;  clear cdm
      gyration_radii = gr.gyration_radii;  clear gr
      orgsnclusters = onc.orgsnclusters';  clear onc
%       if length(cluster_diameters)>length(orgsnclusters), 
%         fprintf('%1.0f > %1.0f \n',length(cluster_diameters),length(orgsnclusters));
%       elseif length(cluster_diameters)<length(orgsnclusters), 
%         fprintf('%1.0f < %1.0f \n',length(cluster_diameters),length(orgsnclusters));
%       else, 
%         fprintf('%1.0f = %1.0f \n',length(cluster_diameters),length(orgsnclusters));
%       end
%       finite = find(cluster_diameters~=inf);
      finite = 1:length(gyration_radii);
      lf = length(finite);
      picks = ceil(lf*rand(samples,1)); %random sampling
      [the_picks,ipicks] = unique(picks); %determine unique picks
      while length(ipicks)~=samples, %make sure the random samples aren't picked from same cluster
        need = samples-length(ipicks); %determine number of new picks needed
        new_picks = zeros(need,1);  %initialize the new picks
        new_picks = ceil(lf*rand(need,1)); %determine new picks
        [the_picks,ipicks] = unique([the_picks; new_picks]); %get unique picks
      end
%       density = orgsnclusters./cluster_diameters;
      lto = log10(orgsnclusters(finite(the_picks))/limit);
%       ltcd = log10(cluster_diameters(finite(the_picks))/op);
      ltgr = log10(gyration_radii(finite(the_picks))/op);
%       [m,b,em,eb] = linear_fit(log10(cluster_diameters(finite)),log10(orgsnclusters(finite)),[],0); %hold on;
      [m,b,em,eb] = linear_fit(log10(gyration_radii(finite)),log10(orgsnclusters(finite)),[],0); %hold on;
%       if run==3, 
%         figure(mu*100000 +dm*100 +SIMOPTS.limit +run);
%         [sm,sb,sem,seb] = linear_fit(ltcd,lto,[],1); %hold on;
%         xlabel('L = log_{10}(cluster\_diameter)'); ylabel('M = log_{10}(orgsnclusters)');
%       else
%         [sm,sb,sem,seb] = linear_fit(ltcd,lto,[],0); %hold on;
      [sm,sb,sem,seb] = linear_fit(ltgr,lto,[],0); %hold on;
%       end
      
%       ltd = log10(density(finite(the_picks)));
%       all_d = log10(orgsnclusters(finite)./density(finite))./log10(cluster_diameters(finite));
%       d_sample = all_d(the_picks);
%       [scd,iscd] = sort(cluster_diameters);
%       [scd,iscd] = sort(gyration_radii);
%       org_scd = orgsnclusters(iscd);
%       fin = find(scd~=inf);
%       cd_large = scd((end-length(fin)-samples):(end-length(fin)));
%       org_large = org_scd((end-length(fin)-samples):(end-length(fin)));
%       [lm,lb,lem,leb] = linear_fit(log10(cd_large),log10(org_large),[],1);  hold on;
    end
    r = r +1;
    M(i,r) = m;
    B(i,r) = b;
    EM(i,r) = em;
    EB(i,r) = eb;
    SM(i,r) = sm;
    SB(i,r) = sb;
    SEM(i,r) = sem;
    SEB(i,r) = seb;
%     LM(i,r) = lm;
%     LB(i,r) = lb;
%     LEM(i,r) = lem;
%     LEB(i,r) = leb;
%     D(i,r,:) = d_sample;
%     D(i,r)
%     fprintf(['density dimension %1.4f +- %1.4f \n'],m,em);
    end
  end
end
end
end