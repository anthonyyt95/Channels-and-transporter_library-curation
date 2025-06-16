%% Extracts RefSeq & NCBI IDs for each gene (human - first pass)

load('2019.11.06_workspace_ictList.mat', 'ICT797')

refFile = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Pan-Experiment Data\Gene ID conversions\mm_RefSeq_NCBI.txt';
list = ICT797;

tic
refStruc = importdata(refFile);
ref = refStruc.textdata(:,[1:7]);
refAppend = [refStruc.textdata(1,8); num2cell(refStruc.data(:,1))];
ref = [ref refAppend];

output = {};
missing = {};
multipleNCBIs = {};
multipleENS = {};
step = 1;
for i = [1:length(list(:,1))]
    aliases = strsplit(list{i,1},' ');
    gene = aliases{1};
    
    loc = find(strcmp(lower(ref(:,6)),lower(gene)));
    if not(isempty(loc))
        geneIDs = ref(loc,:);
        
        % Isolates the RefSeq IDs
        pre_nmIDs = unique(geneIDs(:,5));
        ind = cellfun(@isempty,pre_nmIDs);
        ind = find(ind==1);
        pre_nmIDs(ind,:) = [];
        nmIDs = '';
        for j = [1:length(pre_nmIDs(:,1))]
            if j == 1
                nmIDs = pre_nmIDs{j,1};
            else
                nmIDs = [nmIDs, ' ', pre_nmIDs{j,1}];
            end
        end
        
        % Of the NM_IDs, isolates the ncbi ID
        ncbiID = unique(cell2mat(geneIDs(:,8)));
        
        % Of the NM_IDs, isolates the ensemble IDs
        pre_ensID = unique(geneIDs([2:end],1));
        ensID = '';
        if length(pre_ensID) > 1
            for j = [1:length(pre_ensID)]
                if j == 1
                    ensID = pre_ensID{j};
                else
                    ensID = [ensID,' ',pre_ensID{j}];
                end
            end
        end
                    
        % Creates output
        if not(length(ncbiID) == 1);
            eval(['multipleNCBIs.',gene,' = geneIDs']);
            continue
        else
            output{step,1} = gene;
            output{step,2} = list{i,1};
            output{step,3} = nmIDs;
            output{step,4} = ncbiID;
            output{step,5} = ensID;
            step = step + 1;
        end
    else
        missing = [missing; list{i,1}];
        continue
    end
end
       

toc
clear aliases ans gene geneIDs i ICT797 ind j list loc ncbiID nmIDs
clear pre_nmIDs ref refAppend refFile refStruc step

%% Obtains gene homologs in mouse 

list = input;
refFile = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Pan-Experiment Data\Gene ID conversions\hs_mm_homologs.csv';

tic
refStruc = readtable(refFile);
heading = refStruc.Properties.VariableNames;
refStruc = table2cell(refStruc);
refStruc = [heading; refStruc];

% Separates by taxa
mouseInd = find(cell2mat(refStruc([2:end],3))==10090)+1;
humanInd = find(cell2mat(refStruc([2:end],3))==9606)+1;
mouse = refStruc(mouseInd,:);
human = refStruc(humanInd,:);

% Acquires homologs by gene name
output = {};
step = 1;
missingInHuman = {};
noHomolog = {};
for i = [1:length(input(:,1))]
    hsGene = list{i,1};
    loc = find(strcmp(lower(human(:,4)), lower(hsGene)));
    
    if not(isempty(loc))
        hlogID = human{loc,1};
        
        loc2 = find(cell2mat(mouse(:,1))==hlogID);
        if not(isempty(loc2))
            output{step,1} = hsGene;
            output{step,2} = human{loc,5};
            output{step,3} = human{loc,11};
            output{step,4} = human{loc,13};
            
            output{step,5} = mouse{loc2,4};
            output{step,6} = mouse{loc2,5};
            output{step,7} = mouse{loc2,11};
            output{step,8} = mouse{loc2,13};
            step = step + 1;
        else
            noHomolog = [noHomolog; hsGene];
            continue
        end
    else
        missingInHuman = [missingInHuman; hsGene];
        continue
    end
end

toc
clear heading hlogID hsGene human humanInd i list loc loc2 
clear mouse mouseInd refFile refStruc step


%% Extracts RefSeq & NCBI IDs for each gene (mouse - first pass)

load('2019.11.06_workspace_ictList.mat', 'ICT797')

refFile = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Pan-Experiment Data\Gene ID conversions\mm_RefSeq_NCBI.txt';
list = input;


tic
refStruc = importdata(refFile);
ref = refStruc.textdata(:,[1:7]);
refAppend = [refStruc.textdata(1,8); num2cell(refStruc.data(:,1))];
ref = [ref refAppend];

output = {};
missing = {};
step = 1;
for i = [1:length(list(:,1))]
    aliases = strsplit(list{i,1},' ');
    gene = aliases{1};
    
    loc = find(strcmp(lower(ref(:,6)),lower(gene)));
    if not(isempty(loc))
        geneIDs = ref(loc,:);
        
        % Isolates the RefSeq IDs
        pre_nmIDs = unique(geneIDs(:,5));
        ind = cellfun(@isempty,pre_nmIDs);
        ind = find(ind==1);
        pre_nmIDs(ind,:) = [];
        nmIDs = '';
        for j = [1:length(pre_nmIDs(:,1))]
            if j == 1
                nmIDs = pre_nmIDs{j,1};
            else
                nmIDs = [nmIDs, ' ', pre_nmIDs{j,1}];
            end
        end
        
        % Of the NM_IDs, isolates the ncbi ID
        ncbiID = unique(cell2mat(geneIDs(:,8)));
        if not(length(ncbiID) == 1);
            print('Multiple NCBI IDs');
        else
            output{step,1} = gene;
            output{step,2} = list{i,1};
            output{step,3} = nmIDs;
            output{step,4} = ncbiID;
            step = step + 1;
        end
    else
        missing = [missing; list{i,1}];
        continue
    end
end
       

toc
clear aliases ans gene geneIDs i ICT797 ind j list loc ncbiID nmIDs
clear pre_nmIDs ref refAppend refFile refStruc step

%% Extracts alias information based on NCBI IDs

list = icts;

urlParent = 'https://www.ncbi.nlm.nih.gov/gene/';

output = {};
for i = [1:length(list(:,1))]
    i
    ncbiID = char(string(list{i,1}));
    url = [urlParent, ncbiID];
    source = urlread(url);
    
    loc = strfind(source,'summaryDiv');
    source = source([loc:end]);
    loc = strfind(source,'Also known as');
    if not(isempty(loc))
        source = source([loc:end]);
        loc2 = strfind(source,'</dt>');
        source = source([1:loc2(2)]);
        
        loc = strfind(source,'<dd>');
        source = source([loc+4:end]);
        loc = strfind(source,'</dd>');
        source = source([1:loc-1]);
        
        aliases = regexprep(source,'; ',',');
        output{i,1} = aliases;
    else
        output{i,1} = [];
    end
end

clear aliases i list loc loc2 ncbiID source url urlParent

%% Acquire top-4-pick sequences from gRNA Designer results

dirName = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\Data\2020.02.07_ICTlistForgRNADesign\2020.02.07_gRNAsequences\NegControls.txt';
fileNames = dir(dirName);

output = {};
step = 1;
pickVal = [1 2 3 4];
for i = [1:length(fileNames)]
    fileName = fileNames(i).name;
    loc = strfind(lower(fileName),'.txt');
    if isempty(loc)
        continue
    end
    
    if step == 1
        dirParent = fileNames(i).folder;
    end    
    fileName = [dirParent,'\',fileName];
    
    data = readtable(fileName);
    data = table2cell(data);
    
    % Acquires picked genes
    ind = [];
    for j = [1:length(data(:,1))]
        pick = data{j,51};
        loc = find(pickVal==pick);
        if not(isempty(loc))
            ind = [ind; j];
        end
    end
    outputVal = data(ind,:);
    
    output = [output; outputVal];
end

%% Acquires gRNA sequences given NCBI IDs

filename = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\Data\2020.02.07_ICTlistForgRNADesign\2020.02.07_gRNAsequences\2020.02.14_FinalSequences\NegControls.csv';
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

filename = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\2020.02.24_FinalList.xlsx';
%prime5 = 'CAATTGGAGAAAAGCCTTGTTTG';
%prime3 = 'GTTTTAGAGCTAGGATCCTAGCAAGTT';
prime5 = 'GGAGAAAAGCCTTGTTTG';
prime3 = 'GTTTTAGAGCTAGGATCCTAGC';
colName = 'sgRNA Sequence';

data = importdata(filename);


% Finds appropriate column
header = data.textdata.Sequences(1,:);
col = find(strcmp(colName,header));
gSeq = data.textdata.Sequences([2:end],col);

%% Appends 5' & 3' adapters
prime5List = repelem({prime5},length(gSeq))';
prime3List = repelem({prime3},length(gSeq))';
appendList = [prime5List,gSeq,prime3List];
output = {};
for i = [1:length(appendList(:,1))]
    append = appendList(i,:);
    append = strjoin(append,'');
    output{i,1} = append;
end
%%
delete('output.xlsx');
xlswrite('output.xlsx',output);

clear append appendList col colName data filename gSeq header i prime3
clear prime3List prime5 prime5List








