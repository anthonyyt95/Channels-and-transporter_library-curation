%% Orders ICT scores (from Google Sheets)

data = info;

include = [];
maxScore = max(cell2mat(data(:,13)));
scoreList = cell2mat(data(:,13));
output = {};
for i = [1:maxScore+1]
    metric = (maxScore+1) - i;
    ind = find(scoreList==metric);
    if not(isempty(ind))
        output = [output; data(ind,:)];
    else
        continue
    end
end
xlswrite('output.xlsx',output);

clear data include maxScore scoreList output i metric ind

%% Negative control genes: searches to see if query is expressed in any immune cells

load('E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\06.12.18_B Cell Experiments\Matlab\2019.11.10_MouseRNAseqDatabases.mat', 'haemopediaMouse_full', 'immgenMouseRNAseq_full');

dataSets = {'haemopediaMouse_full';
    'immgenMouseRNAseq_full'};

query = 'SMARCA1';
prctleThreshold = 0.00001;

output = {};
missingDatasets = {};
for i = [1:length(dataSets(:,1))]
    dataSet = dataSets{i,1};
    eval(['data = ',dataSet,'.values;']);
    eval(['genes = ',dataSet,'.genes;']);
    eval(['cells = ',dataSet,'.cells;']);
    
    loc = find(strcmp(lower(genes), lower(query)));
    if isempty(loc)
        missingDatasets = [missingDatasets; dataSet];
        continue
    end
    
    for j = [1:length(cells(1,:))]
        val = data(:,j);
        val = rescale(val,0,100);
        val = val(loc,1);
        val = num2cell(val);
        cell = cells{1,j};
        
        output = [output; cell, val];
    end
end

outputUnder = {};
outputOver = {};
sumNum = 0;
for i = [1:length(output(:,1))]
    val = output{i,2};
    cell = output{i,1};
    if val < prctleThreshold
        sumNum = sumNum + 1;
        outputUnder = [outputUnder; output(i,:)];
    else
        outputOver = [outputOver; output(i,:)];
    end
end
metric1 = (sumNum/length(output(:,1)))*100;
metric2 = (length(outputUnder)/length(outputOver));

clear cell cells data dataSet dataSets genes i j haemopediaMouse_full
clear immgenMouseRNAseq_full loc output prctleThreshold query sumNum val


%% Negative control genes: searches to see if query is expressed in any immune cells

load('E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\06.12.18_B Cell Experiments\Matlab\2019.11.10_MouseRNAseqDatabases.mat', 'haemopediaMouse_full', 'immgenMouseRNAseq_full');

data = immgenMouseRNAseq_full;
queries = input;

genes = immgenMouseRNAseq_full.genes;
values = immgenMouseRNAseq_full.values;

bplotData = [];
missing = {};
for i = [1:length(queries(:,1))]
    gene = queries{i,1};
    
    loc = find(strcmp(lower(genes),lower(gene)));
    if not(isempty(loc))
        bplotData = [bplotData, values(loc,:)'];
    else
        missing = [missing; gene];
    end
end
bplotData = log(bplotData+1);
b = boxplot(bplotData,'PlotStyle','compact');

clear haemopediaMouse_full immgenMouseRNAseq_full data queries genes values 
clear bplotData i gene loc














    
    
    
    
    
    
    



