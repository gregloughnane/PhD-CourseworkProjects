function [] = PythonFileChanger(orig_input,mod_input,odbfilename,rptfilename)
% orig_input is the original input file
% mod_input is the new input file with perturbed values
% odbfilename is the filename for the perturbed Abaqus database file
% rptfilename is the .rpt file to be written by Python

fid = fopen(orig_input,'r'); %open data file
newfile = fopen(mod_input,'w'); %Need to name new file here

end_of_file = 0; %end of file identifier
ender = 0;
linecount = 0;

while end_of_file == 0 && ender == 0

lines = fgetl(fid);

    linecount = linecount+1; % Initiate counter
    
    if linecount < 23
        ScanPrintLinePYTHONONLY(lines,newfile)
    elseif linecount == 24
        Line= sscanf(lines,'%s');
        A = char(Line);
        newLine = strrep(A,'Job-30',odbfilename);
        fprintf(newfile,'%s \n',newLine);
    elseif linecount == 25
        ScanPrintLinePYTHONONLY(lines,newfile)
    elseif linecount == 26
        Line= sscanf(lines,'%s');
        A = char(Line);
        newLine = strrep(A,'Job-30',odbfilename);
        fprintf(newfile,'%s \n',newLine);
    elseif linecount == 27
        Line= sscanf(lines,'%s');
        A = char(Line);
        newLine = strrep(A,'stresses',rptfilename);
        fprintf(newfile,'%s \n',newLine);
    elseif linecount > 27
        ScanPrintLinePYTHONONLY(lines,newfile)
    end
    
        
    end_of_file = feof(fid);
    
end
    fclose(fid);
    fclose(newfile);
end