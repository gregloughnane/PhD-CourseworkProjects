function[] = ScanPrintLinePYTHONONLY(lines,newfile)

    [words,count] = textscan(lines,'%s','Whitespace');
    NumberOfWords = cellfun('length',words);
    
    if NumberOfWords == 1
        A = textscan(lines, '%s ','Delimiter',' \b\t');
        B = strcat(A{1});
    elseif NumberOfWords == 2
        A = textscan(lines, '%s %s','Delimiter',' \b\t');
        B = strcat(A{1},{' '},A{2});
    elseif NumberOfWords == 3
        A = textscan(lines, '%s %s %s','Delimiter',' \b\t');
        B = strcat(A{1},{' '},A{2},{' '},A{3});
    elseif NumberOfWords == 4
        A = textscan(lines, '%s %s %s %s','Delimiter',' \b\t');
        B = strcat(A{1},{' '},A{2},{' '},A{3},{' '},A{4});
    elseif NumberOfWords == 5
        A = textscan(lines, '%s %s %s %s %s','Delimiter',' \b\t');
        B = strcat(A{1},{' '},A{2},{' '},A{3},{' '},A{4},{' '},A{5});
    elseif NumberOfWords == 6
        A = textscan(lines, '%s %s %s %s %s %s','Delimiter',' \b\t');
        B = strcat(A{1},{' '},A{2},{' '},A{3},{' '},A{4},{' '},A{5},{' '},A{6}); 
    elseif NumberOfWords == 7
        A = textscan(lines, '%s %s %s %s %s %s %s','Delimiter',' \b\t');
        B = strcat(A{1},{' '},A{2},{' '},A{3},{' '},A{4},{' '},A{5}... 
                   ,A{6},{' '},A{7});
    elseif NumberOfWords == 8
        A = textscan(lines, '%s %s %s %s %s %s %s %s','Delimiter',' \b\t');
        B = strcat(A{1},{' '},A{2},{' '},A{3},{' '},A{4},{' '},A{5}... 
                   ,A{6},{' '},A{7},{' '},A{8});
   elseif NumberOfWords == 9
        A = textscan(lines, '%s %s %s %s %s %s %s %s %s','Delimiter',' \b\t');
        B = strcat(A{1},{' '},A{2},{' '},A{3},{' '},A{4},{' '},A{5}... 
                   ,A{6},{' '},A{7},{' '},A{8},{' '},A{9});
   elseif NumberOfWords == 14
        A = textscan(lines, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',' \b\t');
        B = strcat(A{1},{' '},A{2},{' '},A{3},{' '},A{4},{' '},A{5}... 
                   ,A{6},{' '},A{7},{' '},A{8},{' '},A{9},{' '},A{10},{' '} ...
                   ,A{11},{' '},A{12},{' '},A{13},{' '},A{14});
    end
    
    fwrite(newfile,B{1,:});
    fprintf(newfile, '\n');
    
end