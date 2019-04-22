function [] = RunAbaqus(file)
% RunAbaqus takes as input an .inp file and submits the job to Abaqus
% Input the file name with single quotes and no .inp

runfile = ['abaqus job=' file];
system(runfile);
% pause(20)
    a = 1;
    while (a ~=0)
        pause(1);
        addlck = [runfile '.lck'];
        a = exist(addlck);
    end
end