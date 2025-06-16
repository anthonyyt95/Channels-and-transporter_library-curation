%% Acquires all the ICTs for the 2nd library + positive & negative controls
%
% Input:
%   fullICTs                X-by-2 cell list of gene names (col 1) and ncbi
%                           refseq ID (col 2)
%   excludRef               X-by-2 cell list of gene names (col 1) and ncbi
%                           refseq ID (col 2). These genes will be removed
%   ctrls                   Control genes (same format as above)

fullICTs = input1;
excludeRef = input2;

% Removes ICTs already included in the 1st sublibrary ('excludeRef')
exclude = [];
fullICTsIDs = cell2mat(fullICTs(:,2));
for i = [1:length(excludeRef(:,1))]
    gene = excludeRef{i,1};
    id = excludeRef{i,2};
    
    loc = find(fullICTsIDs == id);
    if not(isempty(loc))
        exclude = [exclude; loc];
    end
end
fullICTs(exclude,:) = [];
output = fullICTs;

clear fullICTs excludeRef exclude fullICTsIDs i gene id loc

%% Acquires gRNA sequences given NCBI IDs
%
% Input:
%   filename                Char of file location
%   geneList                X-by-2 ell list of gene names (col1) and ncbi
%                           IDs (col2)
%

filename = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\Data\2020.02.07_ICTlistForgRNADesign\2020.02.07_gRNAsequences\2020.02.14_FinalSequences\ICTgRNAs_773.xlsx';
geneList = list;

ref = readtable(filename,'ReadVariableNames',false);
ref = table2cell(ref);

if isnumeric(ref{1,1})
    for i = [1:length(ref(:,1))]
        ref{i,1} = char(string(ref{i,1}));
    end
end

output = {};
missing = {};
for i = [1:length(geneList(:,1))]
    gene = geneList{i,2};
    gene = char(string(gene));
    loc = find(strcmp(lower(ref(:,1)),lower(gene)));
    if not(isempty(loc))
        output = [output; ref(loc,:)];
    else
        missing = [missing; geneList(i,:)]; 
    end
end

clear filename geneList ref i gene loc a

%% Appends adaptors sequences to list of gRNA sequences

filename = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\C2_Datasheet.xlsx';
prime5 = 'CAATTGGAGAAAAGCCTTGTTTG';
prime3 = 'GTTTTAGAGCTAGGATCCTAGCAAGTT';
colName = 'sgRNA Sequence';

data = importdata(filename);

% Finds appropriate column
header = data.textdata.Sequences(1,:);
col = find(strcmp(colName,header));
gSeq = data.textdata.Sequences([2:end],col);

% Appends 5' & 3' adapters
prime5List = repelem({prime5},length(gSeq))';
prime3List = repelem({prime3},length(gSeq))';
appendList = [prime5List,gSeq,prime3List];
output = {};
for i = [1:length(appendList(:,1))]
    append = appendList(i,:);
    append = strjoin(append,'');
    output{i,1} = append;
end

delete('output.xlsx');
xlswrite('output.xlsx',output);

clear append appendList col colName data filename gSeq header i prime3
clear prime3List prime5 prime5List




