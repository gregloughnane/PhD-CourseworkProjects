function [newvb] = InpFilePerturbvb(orig_input,mod_input,perturb)
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
    if linecount < 63608
        ScanPrintLine(lines,newfile);
    elseif linecount == 63608
        ElasticModulusLine= sscanf(lines,'%s');
        A = str2num(ElasticModulusLine);
        Eb = A(1);
        vb = A(2);
        sameEb = Eb;
        newvb = vb*perturb + vb;
        PoissonRatioString = sprintf(', %1.3f',newvb);
        PrintFirst = sprintf(' %1.3e',sameEb);
        newEline = strcat(PrintFirst,PoissonRatioString);
        fprintf(newfile,'%s \n',newEline);
    elseif linecount > 63608
        ScanPrintLine(lines,newfile);
    end
    end_of_file = feof(fid);
    
end
    fclose(fid);
    fclose(newfile);
end