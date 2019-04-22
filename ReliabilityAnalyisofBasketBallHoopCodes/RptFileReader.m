function [maxstress] = RptFileReader(rpt_input)
% rpt_input is the rpt file created after an Abaqus job has ran

fid = fopen(rpt_input,'r'); %open data file

end_of_file = 0; %end of file identifier
ender = 0;
linecount = 0;

while end_of_file == 0 && ender == 0
lines = fgetl(fid);
linecount = linecount+1; % Initiate counter

    if linecount == 19552
        Maxstress= sscanf(lines,'%s');
        A = char(Maxstress);
        if length(A) == 18
            s = str2double(A(1,8:18));
        elseif length(A) == 17
            s = str2double(A(1,8:17));
        end
    end
%       if linecount == 6637
%         Maxstress= sscanf(lines,'%s');
%         A = char(Maxstress);
%         if length(A) == 16
%             s = str2double(A(1,6:16));
%         end
%     end

    
    end_of_file = feof(fid);
end
    
    stress = sprintf('%1.6f',s);
    maxstress =str2num(stress);
    format SHORTE

    fclose(fid);
end