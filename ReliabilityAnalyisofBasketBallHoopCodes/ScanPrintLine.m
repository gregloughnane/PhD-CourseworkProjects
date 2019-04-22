function[] = ScanPrintLine(lines,newfile)

    A = sscanf(lines,'%s'); 
    fprintf(newfile,'%s \n', A);
end