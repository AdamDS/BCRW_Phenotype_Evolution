%% check_data_set.m
%
%
now = datestr(now());
diary_file = ['check_data_set_' now(1:end-3) '.txt'];
change = find(diary_file==':' | diary_file=='-' | diary_file==' ');
diary_file(change) = '_';
diary(diary_file);
diary on;
this_script = 'check_data_set';
fprintf([this_script '\n']);
global SIMOPTS;
N = length(overpop)*length(death_max)*length(mutability)*length(SIMOPTS.SIMS);
totp = zeros(N,1);
bad_xy = [];
bad_sim = cell(N,1);
limit = SIMOPTS.limit;
min_land = 0.5;
max_land = 0.5 +2*((2*max(basic_map_size)) -1) -1;
smax_land = 0.5 +2*((2*min(basic_map_size)) -1) -1;
diag_land = sqrt(max_land^2 +smax_land^2);
i = 0;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  [base_name,dir_name] = NameAndCD(0,do_cd);
  for run = SIMS, 
    new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
    run_name = int2str(run);
    clus_name = cluster_name(base_name);
    
    fprintf(' %s for %s \n',this_script,[clus_name run_name]);
    i = i +1;
    
%     open_all_data;
    [p,gop,perror] = exist_load([new_dir_name 'population_' base_name run_name]);%,1,0);

    %Primary requirement, need population
    if gop, 
      population = p.population;  clear p,  
%       get_all_data;
      
      initialize_check_strings; 
      
      %% Check population values
      checks_p = [IPOP,limit,NGEN]; %check initial, min, & length
      comp_p = {equals,gte,lte};
      [cp,pout] = data_check_A(gop,str_p,population,checks_p,comp_p);
      
%% Check population dependencies
      min_p = min(population);
      max_p = max(population);
      len_p = length(population);
      sum_p = sum(population);  
      len_pmo = len_p-1;
      spmi = sum_p -IPOP;
      rep = 2;
      max_p2nc = ceil(max_p/limit);
      clear population

      [tx,gotx,txerror] = exist_load([new_dir_name 'trace_x_' base_name run_name]);%gop,0);
      if gotx, trace_x = tx.trace_x; clear tx, 
      else, trace_x = []; end
      checks_tx = [sum_p,min_land,max_land]; %check rows, min, & max
      comp_tx = {equals,gte,lte};
      [ctx,txout] = data_check_B(gotx,str_tx,trace_x,checks_tx,comp_tx);
      clear trace_x
      
      [ty,goty,tyerror] = exist_load([new_dir_name 'trace_y_' base_name run_name]);%gop,0);
      if goty, trace_y = ty.trace_y; clear ty, 
      else, trace_y = []; end
      checks_ty = [sum_p,min_land,max_land]; %check rows, min, & max
      comp_ty = {equals,gte,lte};
      [cty,tyout] = data_check_B(goty,str_ty,trace_y,checks_ty,comp_ty);
      clear trace_y      

      [tcs,gotcs,tcserror] = exist_load([new_dir_name 'trace_cluster_seed_' clus_name run_name]);%gop,0);
      if gotcs,  trace_cluster_seed = tcs.trace_cluster_seed;  clear tcs,  
      else, trace_cluster_seed = [];  end
      checks_tcs = [sum_p,(limit-1),1,max_p]; %check rows, cols, min, & max
      comp_tcs = {equals,equals,equals,equals};
      [ctcs,tcsout] = data_check_C(gotcs,str_tcs,trace_cluster_seed,checks_tcs,comp_tcs);
      clear trace_cluster_seed
      
      [sd,gosd,sderror] = exist_load([new_dir_name 'seed_distances_' clus_name run_name]);%gotcs,0);
      if gosd, seed_distances = sd.seed_distances; clear sd,  
      else, seed_distances = []; end
      checks_sd = [spmi,(limit-1),op,diag_land]; %check rows, cols, min, & max
      comp_sd = {equals,equals,gte,lte};
      [csd,sdout] = data_check_C(gosd,str_sd,seed_distances(IPOP+1:end,:),checks_sd,comp_sd);
      clear seed_distances

      [k,gok,kerror] = exist_load([new_dir_name 'kills_' base_name run_name]);%gop,0);
      if gok,  kills = k.kills;  max_kop = max(kills(:,1)); clear k,  
      else, kills = []; max_kop = []; end
      checks_k = [len_p,3,0,2*max_p]; %check rows, cols, min, & max
      comp_k = {equals,equals,gte,lte};
      [ck,kout] = data_check_C(gok,str_k,kills,checks_k,comp_k);
      clear kills
      
      [r,gor,rerror] = exist_load([new_dir_name 'rivalries_' base_name run_name]);%gop,0);
      if gor,  rivalries = r.rivalries;  clear r,  
      else, rivalries = []; end
      checks_r = [len_p,0,max_kop]; %check rows, min, & max
      comp_r = {equals,gte,lte};
      [cr,rout] = data_check_B(gor,str_r,rivalries,checks_r,comp_r);
      clear rivalries

%bacterial parents has two equal columns (parent id = 1 is stored twice,
%once in each column
      [par,gopar,parerror] = exist_load([new_dir_name 'parents_' base_name run_name]);%gop,0);
      if gopar,	parents = par.parents;  clear par,  
      else, parents = []; end
      checks_par = [spmi,rep,1,max_p]; %check rows, cols, min, & max
      comp_par = {equals,equals,gte,lte};
      [cpar,parout] = data_check_C(gopar,str_par,parents,checks_par,comp_par);
      clear parents
      
%% Check cluster dependencies
      [nc,gonc,ncerror] = exist_load([new_dir_name 'num_clusters_' clus_name run_name]);%gop,0);
      if gonc, 
        num_clusters = nc.num_clusters; clear nc,  
        sum_nc = sum(num_clusters);
        max_nc = max(num_clusters);
        ainc = num_clusters(1)+1;
        sum_ncmend = sum(num_clusters(1:(end-1)));
        sum_ncmfirst = sum(num_clusters(2:end));
      else, 
        num_clusters = []; 
        sum_nc = [];
        max_nc = [];
        ainc = [];
        sum_ncmend = [];
        sum_ncmfirst = [];
      end

      checks_nc = [len_p,1,max_p2nc]; %checks rows, min, & max
      comp_nc = {equals,gte,lte};
      [cnc,ncout] = data_check_B(gonc,str_nc,num_clusters,checks_nc,comp_nc);
      clear num_clusters
      
%       if gonc,  
        [onc,goonc,oncerror] = exist_load([new_dir_name 'orgsnclusters_' clus_name run_name]);%gonc,0);
        if goonc,  orgsnclusters = onc.orgsnclusters;  clear onc,  
        else, orgsnclusters = []; end
        checks_onc = [sum_nc,limit,max_p]; %check rows, min, & max
        comp_onc = {equals,gte,lte};
        [conc,oncout] = data_check_B(goonc,str_onc,orgsnclusters,checks_onc,comp_onc);
        clear orgsnclusters

        [tc,gotc,tcerror] = exist_load([new_dir_name 'trace_cluster_' clus_name run_name]);%gonc,0);
        if gotc, trace_cluster = tc.trace_cluster; clear tc,  
        else, trace_cluster = []; end
        checks_tc = [sum_p,1,max_nc]; %check rows, min, & max
        comp_tc = {equals,equals,equals};
        [ctc,tcout] = data_check_B(gotc,str_tc,trace_cluster,checks_tc,comp_tc);
        clear trace_cluster

        [cx,gocx,cxerror] = exist_load([new_dir_name 'centroid_x_' clus_name run_name]);%gonc,0);
        if gocx, centroid_x = cx.centroid_x; clear cx,  
        else, centroid_x = []; end
        checks_cx = [sum_nc,min_land,max_land]; %check rows, min, & max
        comp_cx = {equals,gte,lte};
        [ccx,cxout] = data_check_B(gocx,str_cx,centroid_x,checks_cx,comp_cx);
        clear centroid_x

        [cy,gocy,cyerror] = exist_load([new_dir_name 'centroid_y_' clus_name run_name]);%gonc,0);
        if gocy, centroid_y = cy.centroid_y; clear cy,  
        else, centroid_y = []; end
        checks_cy = [sum_nc,min_land,max_land]; %check rows, min, & max
        comp_cy = {equals,gte,lte};
        [ccy,cyout] = data_check_B(gocy,str_cy,centroid_y,checks_cy,comp_cy);
        clear centroid_y

        [gr,gogr,grerror] = exist_load([new_dir_name 'gyration_radii_' clus_name run_name]);%gocx,0);
        if gogr,  gyration_radii = gr.gyration_radii;  clear gr,  
        else, gyration_radii = []; end
        checks_gr = [sum_nc,op^2,inf]; %check rows, min, & max
        comp_gr = {equals,gte,lte};
        [cgr,grout] = data_check_D(gogr,str_gr,gyration_radii,checks_gr,comp_gr,ainc);
        clear gyration_radii

        [div,godiv,diverror] = exist_load([new_dir_name 'cluster_diversity_' clus_name run_name]);%gocx,0);
        if godiv,  cluster_diversity = div.cluster_diversity;  clear div,  
        else, cluster_diversity = []; end
        checks_div = [sum_nc,op,inf]; %check rows, min, & max
        comp_div = {equals,gte,lte};
        [cdiv,divout] = data_check_D(godiv,str_div,cluster_diversity,checks_div,comp_div,ainc);
        clear cluster_diversity

        [dia,godia,diaerror] = exist_load([new_dir_name 'cluster_diameters_' clus_name run_name]);%gocx,0);
        if godia,  cluster_diameters = dia.cluster_diameters;  clear dia,  
        else, cluster_diameters = []; end
        checks_dia = [sum_nc,op,inf]; %check rows, min, & max
        comp_dia = {equals,gte,lte};
        [cdia,diaout] = data_check_D(godia,str_dia,cluster_diameters,checks_dia,comp_dia,ainc);
        clear cluster_diameters

        [pl,gopl,plerror] = exist_load([new_dir_name 'path_length_' clus_name run_name]);%gocx,0);
        if gopl, path_length = pl.path_length; clear c,  
        else, path_length = []; end
        checks_pl = [sum_nc,op,inf]; %check rows, min, & max
        comp_pl = {equals,gte,lte};
        [cpl,plout] = data_check_B(gopl,str_pl,path_length,checks_pl,comp_pl);
        clear path_length
%       end
%% Check lineages dependencies
      [nd,gond,nderror] = exist_load([new_dir_name 'num_descendants_' base_name run_name]);%gop,0);
      if gond, num_descendants = nd.num_descendants; clear nd,  
      else, num_descendants = []; end
      checks_nd = [IPOP,len_pmo,0,max_p]; %check rows, cols, min, & max
      comp_nd = {equals,equals,gte,lte};
      [cnd,ndout] = data_check_C(gond,str_nd,num_descendants,checks_nd,comp_nd);
      clear num_descendants
      
%       if gonc,  
        [ndc,gondc,ndcerror] = exist_load([new_dir_name 'num_descendant_clusters_' clus_name run_name]);%gotc,0);
        if gondc,  
          num_descendant_clusters = ndc.num_descendant_clusters;  clear ndc,  
          max_ldc = max(sum(num_descendant_clusters,2));
        else, 
          num_descendant_clusters = []; 
          max_ldc = [];
        end
        checks_ndc = [IPOP,len_pmo,0,max_p]; %check rows, cols, min, & max
        comp_ndc = {equals,equals,gte,lte};
        [cndc,ndcout] = data_check_C(gondc,str_ndc,num_descendant_clusters,checks_ndc,comp_ndc);
        clear num_descendant_clusters

        [dc,godc,dcerror] = exist_load([new_dir_name 'descendant_clusters_' clus_name run_name]);%gondc,0);
        if godc, descendant_clusters = dc.descendant_clusters; clear dc,  
        else, descendant_clusters = []; end
        checks_dc = [IPOP,max_ldc,0,max_nc]; %check rows, cols, min, & max
        comp_dc = {equals,equals,gte,lte};
        [cdc,dcout] = data_check_C(godc,str_dc,descendant_clusters,checks_dc,comp_dc);
        clear descendant_clusters

        [ncp,goncp,ncperror] = exist_load([new_dir_name 'num_clusters_produced_' clus_name run_name]);%gotc,0);
        if goncp,  
          num_clusters_produced = ncp.num_clusters_produced;  clear ncp,  
          sum_ncp = sum(num_clusters_produced);
        else, 
          num_clusters_produced = []; 
          sum_ncp = [];
        end
        checks_ncp = [sum_ncmend,0,max_nc]; %check rows, min, & max
        comp_ncp = {equals,gte,lte};
        [cncp,ncpout] = data_check_B(goncp,str_ncp,num_clusters_produced,checks_ncp,comp_ncp);
        clear num_clusters_produced

        [cp,gocp,cperror] = exist_load([new_dir_name 'clusters_produced_' clus_name run_name]);%goncp,0);
        if goncp,  clusters_produced = cp.clusters_produced;  clear cp,  
        else, clusters_produced = []; end
        checks_cp = [sum_ncp,1,max_nc]; %check rows, min, & max
        comp_cp = {equals,equals,lte};
        [ccp,cpout] = data_check_B(gocp,str_cp,clusters_produced,checks_cp,comp_cp);
        clear clusters_produced

        [ncf,goncf,ncferror] = exist_load([new_dir_name 'num_clusters_fused_' clus_name run_name]);%gotc,0);
        if goncf,  
          num_clusters_fused = ncf.num_clusters_fused;  clear ncf,  
          sum_ncf = sum(num_clusters_fused);
        else, 
          num_clusters_fused = []; 
          sum_ncf = [];
        end
        checks_ncf = [sum_ncmfirst,0,max_nc]; %check rows, min, & max
        comp_ncf = {equals,gte,lte};
        [cncf,ncfout] = data_check_B(goncf,str_ncf,num_clusters_fused,checks_ncf,comp_ncf);
        clear num_clusters_fused

        [cf,gocf,cferror] = exist_load([new_dir_name 'clusters_fused_' clus_name run_name]);%goncf,0);
        if gocf, clusters_fused = cf.clusters_fused; clear cf,  
        else, clusters_fused = []; end
        checks_cf = [sum_ncf,0,max_nc]; %check rows, min, & max
        comp_cf = {equals,gte,lte};
        [ccf,cfout] = data_check_B(gocf,str_cf,clusters_fused,checks_cf,comp_cf);
        clear clusters_fused
%       end
    end %gop
    
    if do_passes && gop, 
      print_passes; 
      fprintf([base_name run_name]);% '\n']);
    end
    if do_fails && gop,  
      print_fails;
%       fprintf('\n');
    end
%     print_report(pout);
%     print_report(txout);
%     print_report(tyout);
%     print_report(tcsout);
%     print_report(sdout);
%     print_report(kout);
%     print_report(rout);
%     print_report(parout);
%     print_report(ncout);
%     print_report(oncout);
%     print_report(tcout);
%     print_report(cxout);
%     print_report(cyout);
%     print_report(divout);
%     print_report(diaout);
%     print_report(plout);
%     print_report(ndout);
%     print_report(ndcout);
%     print_report(dcout);
%     print_report(ncpout);
%     print_report(cpout);
%     print_report(ncfout);
%     print_report(cfout);
    %report out
    if do_pause_each,  pause;  end
  end
  fprintf('Done with parameter set! \n');
  if do_pause_set,  pause;  end
end
end
end
diary off;
%% Old code

%       [cp] = data_check_A(population,checks_p);
%       [p1] = print_check_result(str_p,cp(1,1),equals,checks_p(1),perror,cp(2,1));
%       [p2] = print_check_result(str_p,cp(1,2),gte,checks_p(2),perror,cp(2,2));
%       [p3] = print_check_result(str_p,cp(1,3),lte,checks_p(3),perror,cp(2,3));
%       pout = [p1,p2,p3];
%       ipop = population(1);
%       if IPOP==ipop, c_ipop = 1;  end

%       if min_p>=limit, c_minp = 1; end

%       find_p = find(population>=limit);
%       if len_p==length(find_p), c_ngen = 1; end

%       [c_lentx,~,c_mintx,c_maxtx] = check_variable(gotx,trace_x,[sum_p,0,min_land,max_land]);
%       if gotx==1, 
%         len_tx = length(trace_x); 
%         if sum_p==len_tx, c_lentx = 1;  end %length of trace_x
%         min_tx = min(trace_x);
%         if min_land<=min_tx,  c_mintx = 1;  end %minimum of trace_x
%         max_tx = max(trace_x);
%         if max_land>=max_tx,  c_maxtx = 1;  end %maximum of trace_x
%       end

%       [ctx] = data_check_B(gotx,trace_x,checks_tx);
%       [tx1] = print_check_result(str_tx,ctx(1),equals,checks_tx(1),txerror,ctx(2,1));
%       [tx2] = print_check_result(str_tx,ctx(2),lte,checks_tx(2),txerror,ctx(2,2));
%       [tx3] = print_check_result(str_tx,ctx(3),gte,checks_tx(3),txerror,ctx(2,3));
%       txout = [tx1,tx2,tx3];
%       [c_lenty,~,c_minty,c_maxty] = check_variable(goty,trace_y,[sum_p,0,min_land,max_land]);
%       if goty==1, 
%         len_ty = length(trace_y); 
%         if sum_p==len_ty, c_lenty = 1;  end  %length of trace_y
%         min_ty = min(trace_y);
%         if min_land<=min_ty,  c_minty = 1;  end %minimum of trace_y
%         max_ty = max(trace_y);
%         if max_land>=max_ty,  c_maxty = 1;  end %maximum of trace_y
%       end

%       [cty] = data_check_B(goty,trace_y,checks_ty);
%       [ty1] = print_check_result(str_ty,cty(1),equals,checks_ty(1),tyerror,cty(2,1));
%       [ty2] = print_check_result(str_ty,cty(2),lte,checks_ty(2),tyerror,cty(2,2));
%       [ty3] = print_check_result(str_ty,cty(3),gte,checks_ty(3),tyerror,cty(2,3));
%       tyout = [ty1,ty2,ty3];
%       [c_rowtcs,c_coltcs,c_mintcs,c_maxtcs] = check_variable(gotcs,trace_cluster_seed,...
%                                                              [sum_p,limit,1,max_p]);
%       if gotcs==1,  
%         [row_tcs,col_tcs] = size(trace_cluster_seed); 
%         if sum_p==row_tcs,  c_rowtcs = 1; end 
%         if (limit-1)==col_tcs,  c_coltcs = 1; end 
%         min_tcs = min(min(trace_cluster_seed));
%         if 1<=min_tcs,  c_mintcs = 1; end
%         max_tcs = max(max(trace_cluster_seed));
%         if max_p>=max_tcs,  c_maxtcs = 1; end
%       end
      
%       [ctcs] = data_check_C(gotcs,trace_cluster_seed,checks_tcs);
%       [tcs1] = print_check_result(str_tcs,ctcs(1,1),equals,checks_tcs(1),tcserror,ctcs(2,1));
%       [tcs2] = print_check_result(str_tcs,ctcs(1,2),equals,checks_tcs(2),tcserror,ctcs(2,2));
%       [tcs3] = print_check_result(str_tcs,ctcs(1,3),lte,1,tcserror,ctcs(2,3));
%       [tcs4] = print_check_result(str_tcs,ctcs(1,4),gte,ceil(max_p/limit),...
%                                   tcserror,ctcs(2,4));
%       tcsout = [tcs1,tcs2,tcs3,tcs4];
%       [c_lenty,~,c_minty,c_maxty] = check_variable(goty,trace_y,[sum_p,0,op,diag_land]);
%       if gosd==1, 
%         len_sd = size(seed_distances,1);  
%         if sum_p==len_sd, c_lensd = 1;  end 
%         min_sd = min(min(seed_distances));
%         if op<=min_sd,  c_minsd = 1;  end
%         max_sd = max(max(seed_distances));
%         if diag_land>=max_sd, c_maxsd = 1;  end
%       end
      
%       [csd] = data_check_B(gosd,seed_distances,checks_sd);
%       [sd1] = print_check_results(str_sd,csd(1,1),equals,checks_sd(1),sderror,csd(2,1));
%       [sd2] = print_check_results(str_sd,csd(1,2),lte,checks_sd(2),sderror,csd(2,2));
%       [sd3] = print_check_results(str_sd,csd(1,3),gte,checks_sd(3),sderror,csd(2,3));
%       sdout = [sd1,sd2,sd3];
% %       if gok==1,  
% %         len_k = size(kills,1);  
% %         if sum_p==len_k,  c_lenk = 1; end 
%         min_k = min(min(kills));  
% %         if 0<=min_k,  c_mink = 1; end
%         max_k = max(max(kills));
% %         if max_p>=max_k,  c_maxk = 1; end
      
