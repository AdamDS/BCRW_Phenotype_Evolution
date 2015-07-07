function [] = save_populations(base_name,dir_name)
global populations;
p_name = [dir_name 'populations_' base_name];
save(p_name,'populations');
end