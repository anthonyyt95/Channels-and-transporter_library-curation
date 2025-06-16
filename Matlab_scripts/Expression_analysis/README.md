%% Acquires GeneRIF information with NCBI gene IDs

ncbiIDs = list;

parentURL1 = 'https://www.ncbi.nlm.nih.gov/gene/?db=gene&term=';
parentURL2 = '&report=generif&page=';
output = {};
step = 1;
for i = [1:length(ncbiIDs(:,1))]
    ncbiID = char(string(ncbiIDs{i,1}));
    url = [parentURL1,ncbiID,parentURL2,char(string(1))];
    i
    source = urlread(url);
    
    % Determines if Gene RIFs exist for the ID
    loc = strfind(source,'No GeneRIF available for this Gene ID');
    if not(isempty(loc))
        output{i,1} = ncbiID;
        output{i,2} = 'N/A';
        output{i,3} = 'N/A';
        continue
    end
    
    % Acquires total page number
    loc = strfind(source,'<table class="generif jig-ncbigrid');
    source = source([loc:end]);
    loc = strfind(source,'</tfoot>');
    source = source([1:loc(1)]);
    loc = strfind(source,'last</a>');
    if isempty(loc)
        pageTot = '1';
    else
        source = source([1:loc(1)]);
        source = reverse(source);

        loc = strfind(source,'"');
        pageTot = source([loc(1)+1:loc(2)-1]);
        pageTot = reverse(pageTot);
    end
    
    % Collects information through the pages
    rifOut = {};
    paperOut = {};
    for j = [1:str2double(pageTot)]
        url = [parentURL1,ncbiID,parentURL2,char(string(j))];
        source = urlread(url);
        loc = strfind(source,'<table class="generif jig-ncbigrid');
        source = source([loc:end]);
        loc = strfind(source,'<div id="messagearea_bottom"');
        source = source([1:loc]);
        
        loc1 = strfind(source,'<tbody>');
        loc2 = strfind(source,'</tbody>');
        source = source([loc1(1):loc2(1)]);
        
        rowLocStart = strfind(source,'<tr');
        rowLocEnd = strfind(source,'</tr>');
        
        % Acquires table information
        for k = [1:length(rowLocStart)]
            rowInfo = source([rowLocStart(k):rowLocEnd(k)]);
            loc = strfind(rowInfo,'</td>');
            annot = rowInfo([9:loc(1)-1]);
            assocPaper = extractHTMLText(source([rowLocStart(k):rowLocEnd(k)]));
            
            rifOut = [rifOut, annot];
            paperOut = [paperOut, assocPaper];
        end
    end
    
    rifOut = strjoin(rifOut,' ---- ');
    paperOut = strjoin(paperOut,' ---- ');
    
    output{i,1} = ncbiID;
    output{i,2} = rifOut;
    output{i,3} = paperOut;
end

% Create 'output' format
fileOut = {};
for i = [1:length(output(:,1))]
    val = output{i,3};
    val = regexprep(val,'\n',': ');
    fileOut{i,1} = string(output{i,1});
    fileOut{i,2} = string(output{i,2});
    fileOut{i,3} = string(val);
end
output = fileOut;

clear annot ans assocPaper i j k loc loc1 loc2 ncbiID ncbiIDs pageTot 
clear paperOut parentURL1 parentURL2 rifOut rowInfo rowLocEnd rowLocStart
clear source url
    
%% Acquires genes with associations to immune function

keyTerms = {' t cell';
    't-cell';
    ' b cell';
    'b-cell';
    'lymphocyt';
    'macrophag';
    'immun';
    'monocyt';
    'neutrophil';
    'granulocyt';
    'antibod';
    'cytokin';
    'inflamm'};
data = input;

output = {};
step = 1;
for i = [1:length(data(:,1))]
    rifAssoc = strsplit(data{i,3},' ---- ');
    
    rifOut = {};
    for j = [1:length(rifAssoc)]
        rif = rifAssoc{j};
        
        for k = [1:length(keyTerms(:,1))]
            keyTerm = keyTerms{k};
            loc = strfind(lower(rif),keyTerm);
            if not(isempty(loc))
                rifOut = [rifOut, rif];
            end
        end
    end
    rifRelevant = length(rifOut);
    rifTotal = length(rifAssoc);
    
    rifOut = strjoin(rifOut,' ---- ');
    output{step,1} = data{i,1};
    output{step,2} = str2double(data{i,2});
    output{step,3} = rifOut;
    output{step,4} = rifRelevant;
    output{step,5} = rifTotal;
    output{step,6} = (rifRelevant/rifTotal)*100;
    step = step + 1;
end

%% Determines if a list of genes has negligible expression

load('2020.02.08_NewExpressionAnalysis.mat', 'ThSubsetsMouse')
load('2020.02.08_NewExpressionAnalysis.mat', 'immgenMouse')

structure1 = immgenMouse;
structure2 = ThSubsetsMouse;
botPercent = 0.5;
geneList = list;
cellList = clist;

tic
datasets = {'structure1'; 'structure2'};

% Equalizes genes between the 2 datasets
genes1 = structure1.genes;
genes2 = structure2.genes;
include2 = [];
include1 = [];
missing = {};
for i = [1:length(genes1(:,1))]
    gene = genes1{i,1};
    
    loc = find(strcmp(lower(genes2),lower(gene)));
    if not(isempty(loc))
        include2 = [include2; loc(1)];
        include1 = [include1; i];
    else
        missing = [missing; gene];
    end
end
genes2 = genes2(include2,:);
values2 = structure2.values(include2,:);
genes1 = genes1(include1,:);
values1 = structure1.values(include1,:);
genesOut = genes1;

valuesOut = [values1, values2];
cellsOut = [structure1.cells, structure2.cells];
annotOut = [structure1.annot; structure2.annot];

% Specifies cells to be queried
include = [];
for i = [1:length(cellList(:,1))]
    cell = cellList{i,1};
    loc = find(strcmp(lower(cellsOut),lower(cell)));
    if not(isempty(loc))
        include = [include loc];
    end
end
valuesOut = valuesOut(:,include);
cellsOut = cellsOut(1,include);
annotOut = annotOut(include,:);

% Determines if a given gene has some expression
botGenes = {};
botNum = ceil(botPercent*length(genesOut));
for i = [1:length(valuesOut(1,:))]
    sortCell = [genesOut, num2cell(valuesOut(:,i))];
    sortCell = sortrows(sortCell,2,'asc');
    botGenes = [botGenes; sortCell([1:botNum],1)];
end
output = {};
for i = [1:length(geneList(:,1))]
    gene = geneList{i,1};
    loc = find(strcmp(lower(botGenes),lower(gene)));
    numExpressed = length(cellsOut) - length(loc);
    output = [output; gene num2cell(numExpressed)];
end
    
%%

data1 = list1;
data2 = list2;

for i = [1:length(data1(:,1))]
    data1{i,1} = strtrim(data1{i,1});
end
for i = [1:length(data2(:,1))]
    data2{i,1} = strtrim(data2{i,1});
end

data1 = unique(data1);
data2 = unique(data2);

exclude = [];
for i = [1:length(data1(:,1))]
    gene = data1{i,1};
    loc = find(strcmp(lower(data2),lower(gene)));
    if not(isempty(loc))
        exclude = [exclude; i];
    end
end
data1(exclude,:) = [];
output = data1;


    