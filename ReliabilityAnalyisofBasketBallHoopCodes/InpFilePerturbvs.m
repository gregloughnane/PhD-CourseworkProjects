function [newvs] = InpFilePerturbvs(orig_input,mod_input,perturb)
% orig_input is the original input file
% mod_input is the new input file with perturbed values
% perturb is the perturbation parameter used to compute sensitivities via
% forward finite difference

fid = fopen(orig_input,'r'); %open data file
newfile = fopen(mod_input,'w'); %Need to name new file here

end_of_file = 0; %end of file identifier
ender = 0;
linecount = 0;

while end_of_file == 0 && ender == 0

lines = fgetl(fid);

    linecount = linecount+1; % Initiate counter
    if linecount < 63605
        ScanPrintLine(lines,newfile);
    elseif linecount == 63605
        ElasticModulusLine= sscanf(lines,'%s');
        A = str2num(ElasticModulusLine);
        Es = A(1);
        vs = A(2);
        sameEs = Es;
        newvs = vs*perturb + vs;
        PoissonRatioString = sprintf(', %1.3f',newvs);
        PrintFirst = sprintf(' %1.3e',sameEs);
        newEline = strcat(PrintFirst,PoissonRatioString);
        fprintf(newfile,'%s \n',newEline);
    elseif linecount > 63605
        ScanPrintLine(lines,newfile);
    end
    end_of_file = feof(fid);
    
end
    fclose(fid);
    fclose(newfile);
end