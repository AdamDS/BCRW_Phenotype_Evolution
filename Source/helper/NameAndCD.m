%% NameAndCD.m
% function [exp_name,dir_name] = NameAndCD(make_dir,do_cd)
% According to input options and settings (given in SIMOPTS), 
% the experiment name is generated.
%
% make_dir, if set to 1 will make a new directory if the directory indicated by dir_name
% does not exist. If set to 0, then if the directory indicated by dir_name does not exist, 
% the directory will not be created, and an error will be thrown to the command window. If 
% an error is thrown, then whatever is running will continue.
% do_cd, if set to 1 will change to the directory indicated by dir_name. If
% set to 0, then it will not change from the current directory.
%
% exp_name is the base name for simulation ordered as: biological model _ 
% mating scheme _ single mu/competition _ offspring distribution _ 
% landscape movement/fitness level _ landscape type _ landscape topology _ 
% initial population size _ mutability range _ 
% two mu competitors (mu1 population size _ mu1 _ v _ mu2 population size _ mu2)
% random death _ competition _ max generations
% dir_name is the corresponding directory name to exp_name
% 
% Note: if a particular option or setting isn't displayed in the file base name, then it
% is most likely that you have chosen the default setting, OR much less likely, the
% handling is not considered in the source code within this function. You will then have
% to update this function. See below.
%
% NOTE: Changing directories may be done, but it is preferential to use the
% dir_name to call data files by their full directory path+filename. See
% split_cd.m for details. A quick example: if split=10 then all data for 
% simulations 0-9 are in Data\0_9, and data for runs 10-19 are in 
% Data\10_19, etc.
% 
% NOTE: How to customize/use NameAndCD:
% You must plan to have your data organized according to an organizational
% tree of subfolders which specify the simulation options in the folder
% names. All data is then listed under it's specified simulation options
% path name in a folder labeled Data (if split~=0, then Data\a_b, where
% b-a-1=split). See previous NOTE for an example.
% This organization should help with opening data folders quicker while in
% Windows Explorer (may use Linux as well) as well (too much data file 
% storage in one folder can slow loading of data folders, so this helps 
% minimize exploring folders along with keeping the data nice and tidy). 
% For help organizing your folder tree, see NameAndCD.m for the dir_name 
% variable which tells Matlab where your data should be. Subsequently, 
% the source parameter in main_babies_#.m should contain everything up to 
% the folder where dir_name begins to accumulate directories. For example, 
% I have my data tree starting at
% 'R:\Babies_Root\' on my lab computers & 'C:\Babies_Root\' on my laptop.
%
% -ADS

function [exp_name,dir_name] = NameAndCD2(make_dir,do_cd),
global SIMOPTS;
if ~exist('do_cd','var'),  do_cd = 1;  end
exp_name = [];  dir_name = [];

%% Source directory
if SIMOPTS.source(1)==0 && SIMOPTS.OS==0, % Windows root
  dir_name = [dir_name 'C:\Babies_Root\'];
elseif SIMOPTS.source(1)==0 && SIMOPTS.OS==1, % Unix Bortas root
  dir_name = [dir_name '/export/home/scotta/Desktop/Babies_Root'];
elseif SIMOPTS.source(1)==1 && SIMOPTS.OS==1, % Unix Bortas root
  dir_name = [dir_name '/export/home/kingd/Desktop/Babies_Root'];
elseif SIMOPTS.source(1)==2 && SIMOPTS.OS==1, % Unix Bortas root
  dir_name = [dir_name '/export/home/bahars/Desktop/Babies_Root'];
elseif SIMOPTS.OS==0, % Custom Windows location
  dir_name = [dir_name SIMOPTS.source];
end

if SIMOPTS.loaded==1, exp_name = [exp_name SIMOPTS.load_name];  end

if SIMOPTS.no_bio,  exp_name = [exp_name 'Faux_'];  dir_name = [dir_name 'Faux\'];  end

%% Reproduction type
if SIMOPTS.reproduction==0, 
  exp_name = [exp_name]; dir_name = [dir_name 'Assortative_Mating\'];
elseif SIMOPTS.reproduction==1, 
  exp_name = 'Bacterial_'; dir_name = [dir_name 'Bacterial\'];
elseif SIMOPTS.reproduction==2, 
  exp_name = 'Random_Mating_'; dir_name = [dir_name 'Random_Mating\']; 
end

%% Mutability type
n2sb = proper_name(SIMOPTS.mu);
if SIMOPTS.exp_type==0, 
  exp_name = [exp_name 'Mu_' n2sb '_'];  dir_name = [dir_name 'Mu_x\'];
elseif SIMOPTS.exp_type==1, 
  exp_name = [exp_name 'Comp_']; dir_name = [dir_name 'Comp\']; 
elseif SIMOPTS.exp_type==3, 
  exp_name = 'Confined_';  dir_name = [dir_name 'Confined\'];  
end

%% Offspring distribution type
if SIMOPTS.distribution==0, 
  exp_name = [exp_name 'Uniform_']; dir_name = [dir_name 'Uniform\'];
else, 
  exp_name = [exp_name 'Normal_']; dir_name = [dir_name 'Normal\'];  
end

%% Landscape movement type
if SIMOPTS.landscape_movement>SIMOPTS.NGEN && ...
    SIMOPTS.landscape_heights(1)~=SIMOPTS.landscape_heights(2),
  exp_name = [exp_name 'Frozenscape_']; dir_name = [dir_name 'Frozenscape\'];
elseif SIMOPTS.landscape_heights(1)~=SIMOPTS.landscape_heights(2), 
    if SIMOPTS.shock==0, 
      exp_name = [exp_name 'Shifting_']; dir_name = [dir_name 'Shifting\'];
    elseif SIMOPTS.shock==1 || SIMOPTS.landscape_movement~=2, 
        exp_name = [exp_name int2str(SIMOPTS.landscape_movement) '_Shock_'];
        dir_name = [dir_name int2str(SIMOPTS.landscape_movement) '_Shock\'];
    end
elseif SIMOPTS.landscape_heights(1)==SIMOPTS.landscape_heights(2), 
  exp_name = [exp_name int2str(SIMOPTS.landscape_heights(1)) '_Flatscape_']; 
  dir_name = [dir_name int2str(SIMOPTS.landscape_heights(1)) '_Flatscape\']; 
end

%% Landscape size
if SIMOPTS.basic_map_size(1)~=12 || SIMOPTS.basic_map_size(2)~=12, 
  exp_name = [exp_name int2str(SIMOPTS.basic_map_size(1)) 'x' ...
    int2str(SIMOPTS.basic_map_size(2)) '_basic_map_size_'];
  dir_name = [dir_name int2str(SIMOPTS.basic_map_size(1)) 'x' ...
    int2str(SIMOPTS.basic_map_size(2)) '_basic_map_size\'];
end
if SIMOPTS.linear,  
  exp_name = [exp_name 'linear_'];  dir_name = [dir_name 'linear\'];
end

%% Landscape topology type
if SIMOPTS.periodic(1), 
  if SIMOPTS.periodic(2), 
    exp_name = [exp_name 'toroid_'];  dir_name = [dir_name 'toroid\'];
  else, 
    exp_name = [exp_name 'xcylinder_']; dir_name = [dir_name 'xcylinder\'];
  end
else, 
  if SIMOPTS.periodic(2), 
    exp_name = [exp_name 'ycylinder_']; dir_name = [dir_name 'ycylinder\'];
  end
end

%% Initial population size
if SIMOPTS.IPOP~=300,
  exp_name = [exp_name int2str(SIMOPTS.IPOP) '_IPOP_']; 
  dir_name = [dir_name int2str(SIMOPTS.IPOP) '_IPOP\'];
end

%% Mutability range
if SIMOPTS.range~=1 && SIMOPTS.exp_type==1, 
  exp_name = [exp_name int2str(SIMOPTS.range) '_range_']; 
  dir_name = [dir_name int2str(SIMOPTS.range) '_range\'];
end

%% Mutability dual
if SIMOPTS.exp_type==2, 
  exp_name = [exp_name int2str(SIMOPTS.bi(2,1)) '_' int2str(SIMOPTS.bi(1,1)*100)...
    '_v_' int2str(SIMOPTS.bi(2,2)) '_' int2str(SIMOPTS.bi(1,2)*100) '_'];
  dir_name = [dir_name int2str(SIMOPTS.bi(2,1)) '_' int2str(SIMOPTS.bi(1,1)*100)...
    '_v_' int2str(SIMOPTS.bi(2,2)) '_' int2str(SIMOPTS.bi(1,2)*100) '\'];
end

%% Random death size & type
if SIMOPTS.dm~=0.7, 
  n2sdm = proper_name(SIMOPTS.dm);
  if SIMOPTS.indiv_death==0, 
    exp_name = [exp_name n2sdm '_death_max_'];
    dir_name = [dir_name n2sdm '_death_max\'];
  elseif SIMOPTS.indiv_death==1, 
    exp_name = [exp_name n2sdm '_indiv_death_max_'];
    dir_name = [dir_name n2sdm '_indiv_death_max\'];
  end
end

%% Competition limit size
if SIMOPTS.op~=0.25, 
  n2sop = proper_name(SIMOPTS.op);
  exp_name = [exp_name n2sop '_overpop_'];
  dir_name = [dir_name n2sop '_overpop\'];
end

%% Competition type
if SIMOPTS.random_walk==1, 
  exp_name = [exp_name 'BARW_'];
  dir_name = [dir_name 'BARW\'];
end

%% Maximum number of generations
if SIMOPTS.NGEN~=2000, 
  exp_name = [exp_name int2str(SIMOPTS.NGEN) '_NGEN_'];
  dir_name = [dir_name int2str(SIMOPTS.NGEN) '_NGEN\'];
end

%% Lifetimes & Relaxations fix
if SIMOPTS.only_lt, 
  sims = SIMOPTS.SIMS;
  first = sims(1);
  last = sims(end);
  exp_name = [exp_name int2str(first) '_' int2str(last)];
end

dir_name = [dir_name 'Data\'];

if SIMOPTS.loaded,  dir_name = [dir_name 'Loaded\'];  end

if SIMOPTS.for_external==1, dir_name = [dir_name 'needs_sorting\']; end

%% Operating system slash fix
if SIMOPTS.OS(1)==1, 
  backslash = find(dir_name=='\');
  dir_name(backslash) = '/';
end

%% Create & change directory
if make_dir && ~exist(dir_name,'dir'),  %if want make_dir but dir DNE
  try,  mkdir(dir_name);
  catch error,  fprintf('Cannot mkdir: %s \n %s \n',dir_name,error);  end
end
if do_cd, 
  try,  cd(dir_name); 
  catch error,  fprintf('Cannot cd: %s \n %s \n',dir_name,error);  end
end
% if exist(dir_name,'dir')~=7,  %if dir_name does not exist
%   if make_dir,  mkdir(dir_name);  end
%   if do_cd, %wanted to cd
%     fprintf('ERROR: directory %s DNE',dir_name); %DNE = Does Not Exist
%   end
% else, %if dir_name does exist
%   if do_cd, cd(dir_name); end
% end
end
%% Old code
% if do_cd, 
%   try 
%     if make_dir,  
%       mkdir(dir_name);
%     end
%     if do_cd, 
%       cd(dir_name);
%     end
%     if SIMOPTS.OS==0, 
%       if ~SIMOPTS.loaded, 
%         cd([dir_name 'Data\']);
%       else, 
%         cd([dir_name 'Data\Loaded/']);
%       end
%     elseif SIMOPTS.OS==1, 
%       if ~SIMOPTS.loaded, 
%         cd([dir_name 'Data/']);
%       else, 
%         cd([dir_name 'Data/Loaded/']);
%       end
%     end
%   catch
%     if make_dir,
%       if SIMOPTS.OS==0, 
%         if ~SIMOPTS.loaded, 
%           mkdir([dir_name 'Data\']);
%           cd([dir_name 'Data\']);
%         else, 
%           mkdir([dir_name 'Data\Loaded\']);
%           cd([dir_name 'Data\Loaded\']);
%         end
%       elseif SIMOPTS.OS==1, 
%         if ~SIMOPTS.loaded, 
%           mkdir([dir_name 'Data/']);
%           cd([dir_name 'Data/']);
%         else, 
%           mkdir([dir_name 'Data/Loaded/']);
%           cd([dir_name 'Data/Loaded/']);
%         end
%       end
%     else, 
%       error = [dir_name 'directory does not exist']
%     end
%   end
% end