function [] = RemoveSpaces(orig_input,mod_input)

fid = fopen(orig_input,'r'); %open data file
newfile = fopen(mod_input,'w'); %Need to name new file here

end_of_file = 0; %end of file identifier
ender = 0;
linecount = 0;

while end_of_file == 0 && ender == 0

lines = fgetl(fid);

    linecount = linecount+1; % Initiate counter
    if linecount <= 39589
        ScanPrintLine(lines,newfile);
    end
    
    end_of_file = feof(fid);
    end
    fclose(fid);
    fclose(newfile);
end