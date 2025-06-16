%% Create FASTA file from Excel

data = shRNAref;
ind = length(data(1,:));

% Removes duplicates
uniqueIDs = unique(data(:,1));
uniqueRef = {};
for i = [1:length(uniqueIDs(:,1))]
    uniqueID = uniqueIDs{i,1};
    loc = find(strcmp(lower(data(:,1)),lower(uniqueID)));
    uniqueRef = [uniqueRef; uniqueID data(loc(1),2)];
end

txtFile = '';
for i = [1:length(uniqueRef(:,1))]
    val1 = uniqueRef{i,1};
    val2 = uniqueRef{i,2};
    line = ['>',val1,'?',val2];
    
    if i == 1
        txtFile = [txtFile,line];
        continue
    end
    
    txtFile = [txtFile,'???',line];
end
txtFile = regexprep(txtFile,'?','\n');
dlmwrite('shRNAReference.fasta',txtFile,'delimiter','');

clear data i ind line txtFile val1 val2 uniqueIDs uniqueRef loc uniqueID

%% Trims down a fastq file

filename = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\samlpe4_S4_L001_R1_001.fastq';
cutLead = 50;
cutLag = 10;

fid = fopen(filename,'r');
fid2 = fopen('output.fastq','w');
i = 1;
tline = fgetl(fid);
while ischar(tline)
    if i == 1
        if sum(ismember(tline,'@'))
            fprintf(fid2,'%s\n',tline);
        else
            tline = tline([46:end-6]);
            fprintf(fid2,'%s\n',tline);
        end
        i = i+1;
        continue
    end
    
    tline = fgetl(fid);
    if sum(ismember(tline,'@'))
        fprintf(fid2,'%s\n',tline);
    else
        tline = tline([46:end-6]);
        fprintf(fid2,'%s\n',tline);
    end
    i = i+1;
    
    if rem(i,10000)==0
        i
    end
    
end
fclose(fid);
fclose(fid2);












