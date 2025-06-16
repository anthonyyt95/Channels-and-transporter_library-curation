%% Acquires positive controls based on CERES 

load('C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\04.10.19_CancerDatabases&Analysis\Matlab\20190426_workspace_CRISPRi.mat');

data = ceres;

values = ceres.values;
genes = ceres.genes;
avgValues = mean(values,2);
sortOut = [genes, num2cell(avgValues)];
output = sortrows(sortOut,2,'asc');
    
clear ceres data values genes avgValues sortOut 

%% Acquires most depleting CERES genes associated with the most depleting genes from the REGNASE-1 paper

load('C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\04.10.19_CancerDatabases&Analysis\Matlab\20190426_workspace_CRISPRi.mat');
fileName = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\Data\20200223_REGNASE1TcellWholeGenomeScreen\41586_2019_1821_MOESM5_ESM.xlsx';

regnaseData = importdata(fileName);
ceresData = ceres;

% Sorts CERES data
values = ceresData.values;
genes = ceresData.genes;
avgValues = mean(values,2);
sortOut = [genes, num2cell(avgValues)];
ceresSorted = sortrows(sortOut,2,'asc');

% Sorts REGNASE-1 Brie screen data
values = regnaseData.data.Sheet1(:,3);
genes = regnaseData.textdata.Sheet1([2:end],2);
sortOut = [genes, num2cell(values)];
regnaseSorted = sortrows(sortOut,2,'asc');

% Acquires genes with depleting effects in both
output = [];
step = 1;
for i = [1:length(ceresSorted([1:100],1))]
    gene = ceresSorted{i,1};
    loc = find(strcmp(lower(regnaseSorted(:,1)),lower(gene)));
    if not(isempty(loc))
        output{step,1} = gene;
        output{step,2} = ceresSorted{i,2};
        output{step,3} = regnaseSorted{loc,2};
        step = step + 1;
    end
end

%% Boxplot graph

load('C:\Users\antho\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\04.10.19_CancerDatabases&Analysis\Matlab\20190426_workspace_CRISPRi.mat');

geneQuery = input;
data = ceres;

values = ceres.values;
genes = ceres.genes;
plotValues = [];
plotGenes = {};
for i = [1:length(geneQuery(:,1))]
    gene = geneQuery{i,1};
    loc = find(strcmp(lower(genes),lower(gene)));
    if not(isempty(loc))
        plotValues = [plotValues, values(loc,:)'];
        plotGenes = [plotGenes; gene];
    end
end
figure
hold
b = boxplot(plotValues,'PlotStyle','compact');
xRange = get(gca,'xlim');
line = plot(xRange,[-0.5,-0.5]);
line.Color = [0.5 0.5 0.5];
set(gca,'FontSize',15);
ylim([-3,1]);
    

    