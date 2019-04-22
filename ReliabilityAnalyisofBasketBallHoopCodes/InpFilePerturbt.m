function [newt] = InpFilePerturbt(orig_input,mod_input,perturb)
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
    if linecount < 29
        ScanPrintLine(lines,newfile);
    elseif linecount == 29
        Line= sscanf(lines,'%s');
        A = char(Line);
        if length(A) == 10;
            shiftt = A(1,6:10);
            B = A(1,7:10);
            C = str2num(B);
            D = abs(C)
            Origt = sprintf(',%1.1f',shiftt);
            Origtelement = sprintf(',%1.1f',D);
        elseif length(A) == 11
            shiftt = A(1,6:11);
            B = A(1,7:11);
            C = str2num(B);
            D = abs(C)
            Origt = sprintf(',%1.2f',shiftt);
            Origtelement = sprintf(',%1.1f',D);
        elseif length(A) == 12
            shiftt = A(1,6:12);
            B = A(1,7:12);
            C = str2num(B);
            D = abs(C)
            Origt = sprintf(',%1.3f',shiftt);
            Origtelement = sprintf(',%1.1f',D);
         elseif length(A) == 14
            shiftt = A(1,6:14);
            B = A(1,7:14);
            C = str2num(B);
            D = abs(C)
            Origt = sprintf(',%1.5f',shiftt);
            Origtelement = sprintf(',%1.1f',D);
        end
        newt = -(D*perturb + D)
        Replacementt = sprintf(',%1.5f',newt);
        newLine = strrep(A,Origt,Replacementt);
        fprintf(newfile,'%s \n',newLine);
    elseif linecount == 30
        ScanPrintLine(lines,newfile);
    elseif linecount >= 31 && linecount < 19810
        ScanPrintLine(lines,newfile);
    elseif linecount >= 19810 && linecount < 39589
        Line= sscanf(lines,'%s');
        A = char(Line);
        newt = D*perturb + D;
        Replacementt = sprintf(',%1.5f',newt);
        newLine = strrep(A,Origtelement,Replacementt);
        fprintf(newfile,'%s \n',newLine);
    elseif linecount >= 39589
        ScanPrintLine(lines,newfile);
    end
    end_of_file = feof(fid);
    
end
    fclose(fid);
    fclose(newfile);
end