function [outputF] = InpFilePerturbF(orig_input,mod_input,perturb)
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
    if linecount < 63628
        ScanPrintLine(lines,newfile);
    elseif linecount == 63628
        ForceLine= sscanf(lines,'%s');
        A = char(ForceLine);
        if length(A) == 23
           force = str2num(A(1,17:23));
        elseif length(A) == 24
            force = str2num(A(1,17:24));
        end
        newF = -(force*perturb + force);
        outputF = abs(newF);
        PreString = '_PickedSet52, 2,';
        PrintSecond = sprintf(' %1.2e',newF);
        newFline = strcat(PreString,PrintSecond);
        fprintf(newfile,'%s \n',newFline);
    elseif linecount > 63628
        ScanPrintLine(lines,newfile);
    end
    end_of_file = feof(fid);
    
end
    fclose(fid);
    fclose(newfile);
end