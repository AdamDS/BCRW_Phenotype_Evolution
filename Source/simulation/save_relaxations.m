function [] = save_relaxations(base_name,dir_name)
global lifetimes;
relaxations = lifetimes;
r_name = [dir_name 'relaxations_' base_name];
save(r_name,'relaxations');
end