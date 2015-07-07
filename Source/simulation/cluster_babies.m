%always set these variables
matlab_ver      = 'R2012b';    % (MATLAB release supported by your license) R2009a R2009b R2010a
email           = 'adsbnc@mail.umsl.edu'; % your email address
email_opt       = 'a';       % qsub email options
h_rt            = '1:07:00'; % hard wall time
vf              = '1G';   % Amount of memory need per task
queue           = 'all.q'; % specify queue
min_cpu_slots   = 1;          % Min number of cpu slots needed for the job
max_cpu_slots   = 12;          % Max number of cpu slots needed for the job
disp('Please wait.. Sending job data to the Cluster.... ')
% Configure the scheduler - Do NOT modify these
%sge_options = ['-l vf=', vf, ',h_rt=', h_rt,' -m ', email_opt, ' -M ', email, ' -q ' , queue ];
sge_options = ['-l vf=', vf, ',h_rt=', h_rt,' -m ', email_opt, ' -M ', email ];
% SGEClusterInfo.setExtraParameter(sge_options);
sched = findResource();
% End of scheduler configuration
get(sched)
job2 = createJob(sched);
tic
% start of user specific commands
job2= batch('test_parfor', 'matlabpool', max_cpu_slots, 'FileDependencies', {'test_parfor.m'});
disp('Job submitted..')
datestr(clock)
% The following commands can be run once the job is submitted to view the results
disp ('Job sent to the cluster')
disp('USE >> waitForState(job)       to wait for job to be finished')
disp('USE >> job.State               to see job state')
disp('USE >> load(job,variable)      to load variables back in the workspace      OR') 
disp('USE >> results = getAllOutputArguments(job)      to load variables back in the workspace      AND') 
disp('USE >> results{:}               to see the results')

%
% %always set these variables
% matlab_ver      = 'R2012b';    % (MATLAB release supported by your license) R2009a R2009b R2010a
% email           = 'YOUREMAILID'; % your email address
% email_opt       = 'a';       % qsub email options
% h_rt            = '1:07:00'; % hard wall time
% vf              = '2G';   % Amount of memory need per task
% queue           = 'all.q' % specify queue
% min_cpu_slots   = 4;          % Min number of cpu slots needed for the job
% max_cpu_slots   = 32;          % Max number of cpu slots needed for the job
% disp('Please wait.. Sending job data to the Cluster.... ')
% % Configure the scheduler - Do NOT modify these
% %sge_options = ['-l vf=', vf, ',h_rt=', h_rt,' -m ', email_opt, ' -M ', email, ' -q ' , queue ];
% sge_options = ['-l vf=', vf, ',h_rt=', h_rt,' -m ', email_opt, ' -M ', email ];
% SGEClusterInfo.setExtraParameter(sge_options);
% sched = findResource();
% % End of scheduler configuration
% get(sched)
% job2 = createJob(sched);
% tic
% % start of user specific commands
% job2= batch('USERmFile', 'matlabpool', max_cpu_slots, 'FileDependencies', {'USERFUNCTIONS'});
% disp('Job submitted..')
% datestr(clock)
% % The following commands can be run once the job is submitted to view the results
% disp ('Job sent to the cluster')
% disp('USE >> waitForState(job)       to wait for job to be finished')
% disp('USE >> job.State               to see job state')
% disp('USE >> load(job,variable)      to load variables back in the workspace      OR') 
% disp('USE >> results = getAllOutputArguments(job)      to load variables back in the workspace      AND') 
% disp('USE >> results{:}               to see the results')

%
% % Example MATLAB submission script for running a parallel job on the 
% % SNS Aurora Cluster
% % 24-Feb-2012 Prentice Bisbal, prentice@ias.edu
% 
% % Do not modify these lines
% sched = findResource('scheduler', 'type', 'generic');
% set(sched,'ClusterMatlabRoot','/state/partition1/apps/MATLAB/R2012b');
% set(sched,'ClusterOsType','unix');
% set(sched,'ClusterSize',16);
% set(sched,'DestroyJobFcn',@destroyJobFcn);
% set(sched,'GetJobStateFcn',@getJobStateFcn);
% set(sched,'HasSharedFilesystem',false);
% 
% % Modify these lines to suit your job requirements.
% set(sched,'DataLocation','C:\Users\amviot\Research\Programming\Babies\Jobs');
% h_rt = '00:05:00';
% exclusive = 'false';
% 
% % Do not modify this line
% set(sched,'ParallelSubmitFcn',{@parallelSubmitFcn, h_rt, exclusive});
% 
% % Create parallel job with default settings. Which of these two functions 
% % you use depends on what functions you will be using in your parallel job. 
% %pjob = createParallelJob(sched);
% pjob = createMatlabPoolJob(sched);
% 
% % Specify the number of workers required for execution of your job
% set(pjob, 'MinimumNumberOfWorkers', 1);
% set(pjob, 'MaximumNumberOfWorkers', 8);
% 
% % Files needed by the job
% set(pjob, 'FileDependencies', {'test_parfor.m'});
% 
% % Add a task to the job. In this example, I'm calling MATLAB's Sum 
% % function to add 1+1. Replace @sum with a name your MATLAB program
% %createTask(pjob, @sum, 1, {[1 2 3]});
% createTask(pjob, @test_parfor, 1, {});
% 
% % Submit the job to the cluster
% submit(pjob);
% 
% % Wait for the job to finish running, and retrieve the results.
% % This is optional. Your program will block here until the parallel
% % job completes. If your program is writing it's results to file, you
% % many not want this, or you might want to move this further down in your
% % code, so you can do other stuff while pjob runs.
% waitForState(pjob, 'finished');
% results = getAllOutputArguments(pjob);
% 
% % This checks for errors from individual tasks and reports them.
% % very useful for debugging
% errmsgs = get(pjob.Tasks, {'ErrorMessage'});
% nonempty = ~cellfun(@isempty, errmsgs);
% celldisp(errmsgs(nonempty));
% 
% % Display the results
% disp(results);
% 
% % Destroy job
% % For parallel jobs, I recommend NOT using the destroy command, since it
% % causes the SGE jobs to exit with an Error due to a race condition. If you
% % insist on using it to clean up the 'Job' files and subdirectories in your
% % working directory, you must include the pause statement to avoid the job 
% % finishing in SGE with a Error. 
% %pause(16);
% %destroy(pjob);
% % 
% % 
% % % % % babiesFiles = {'Simulations','try_catch_load','setInitialMutabilities'
% % % % useProfile = 'Bortas';
% % % % useNCores = 8;
% % % % clus = parcluster(useProfile);
% % % % job = createCommunicatingJob(clus);
% % % % task = createTask(job, @main_babies_func);
% % % % communicatingSubmitFcn(job);
% % % % % job = batch('main_babies_0',...
% % % % %             'Profile',useProfile,...
% % % % %             'Matlabpool',useNCores);%,...
% % % % % %             'AttachedFiles',babiesFiles);
% % % % wait(job);
% % % % results = fetchOutputs(job);