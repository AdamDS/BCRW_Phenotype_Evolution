>> sched = findResource('scheduler', 'type', 'generic')
>> set(sched,'DataLocation','/home/user')
>> set(sched,'HasSharedFilesystem',true)
>> set(sched,'ClusterMatlabRoot','/usr/local/matlab-r14sp3')
>> set(sched,'SubmitFcn',@submitFunc)
>> j = createJob(sched);
>> createTask(j, @sum, 1, {[1 1]});
>> createTask(j, @sum, 1, {[2 2]});
>> createTask(j, @sum, 1, {[3 3]});
>> submit(j);

(source /usr/local/sge/default/common/settings.csh ; qsub -o "/home/user/Job2/Task1.log" /usr/local/matlab-r14sp3/toolbox/local/submit)
your job 53208 ("submit") has been submitted

(source /usr/local/sge/default/common/settings.csh ; qsub -o "/home/user/Job2/Task2.log" /usr/local/matlab-r14sp3/toolbox/local/submit)
your job 53209 ("submit") has been submitted

(source /usr/local/sge/default/common/settings.csh ; qsub -o "/home/user/Job2/Task3.log" /usr/local/matlab-r14sp3/toolbox/local/submit)
your job 53210 ("submit") has been submitted

>> waitForState(j)
>> results = getAllOutputArguments(j);
>> results

results = 

    [2]
    [4]
    [6]

>> 