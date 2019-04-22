function [] = InpFileCreator(orig_input,mod_input,newt,newF,newEb,newvb,newEs,newvs)
% orig_input is the original input file
% mod_input is the new input file with perturbed values
% newE is the new design variable E value
% newF is the new design variable E value

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
            D = abs(C);
            Origtelement = sprintf(',%1.1f',D);
        elseif length(A) == 11
            shiftt = A(1,6:11);
            B = A(1,7:11);
            C = str2num(B);
            D = abs(C);
            Origtelement = sprintf(',%1.1f',D);
        elseif length(A) == 12
            shiftt = A(1,6:12);
            B = A(1,7:12);
            C = str2num(B);
            D = abs(C);
            Origtelement = sprintf(',%1.1f',D);
        end
        Replacementt = sprintf(',%1.5f',-newt);
        newLine = strrep(A,shiftt,Replacementt);
        fprintf(newfile,'%s \n',newLine);
    elseif linecount == 30
        ScanPrintLine(lines,newfile);
    elseif linecount > 30 && linecount < 19810
        ScanPrintLine(lines,newfile);
    elseif linecount >= 19810 && linecount < 39589
        Line= sscanf(lines,'%s');
        A = char(Line);
        Replacementt = sprintf(',%1.5f',newt);
        newLine = strrep(A,Origtelement,Replacementt);
        fprintf(newfile,'%s \n',newLine);
    elseif linecount >= 19810 && linecount < 63605
        ScanPrintLine(lines,newfile);
    elseif linecount == 63605
        PoissonRatioString = sprintf(', %1.3f',newvs);
        PrintFirst = sprintf(' %1.3e',newEs);
        newEline = strcat(PrintFirst,PoissonRatioString);
        fprintf(newfile,'%s \n',newEline);
    elseif linecount > 63605 && linecount < 63608
        ScanPrintLine(lines,newfile);            
    elseif linecount == 63608
        PoissonRatioString = sprintf(', %1.3f',newvb);
        PrintFirst = sprintf(' %1.3e',newEb);
        newEline = strcat(PrintFirst,PoissonRatioString);
        fprintf(newfile,'%s \n',newEline);
    elseif linecount > 63608 && linecount < 63628
        ScanPrintLine(lines,newfile);
    elseif linecount == 63628
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