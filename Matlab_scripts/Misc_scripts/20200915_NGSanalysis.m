%% Create fasta file

%refSeqName = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\C1_Datasheet.xlsx';
%refSeqName = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\C2_Datasheet.xlsx'
%refSeqName = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\06.26.18_Misc\2020.04.06_MetabolicTransporterList\C3_datasheet.xlsx';
refSeqName = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\Data\2020.10.01_AdjunctLibraryCloning\gRNASeq_LargeLibraries\2020.10.05_C3Large.xlsx';
fileOutName = 'C3.2reference.fasta';

refSeq = importdata(refSeqName);

geneNCBI = refSeq.data.Sequences(:,1);
geneNames = refSeq.textdata.Sequences([2:end],5);
pickOrder = refSeq.data.Sequences(:,51);
genegRNA = refSeq.textdata.Sequences([2:end],20);
output = {};
emptyLine = {''};
for i = [1:length(geneNCBI(:,1))]
    geneName = regexprep(geneNames{i},' ','-');
    line1 = ['>',geneName,'_',num2str(pickOrder(i)),' [gene=',geneNames{i},'] [ncbi=',num2str(geneNCBI(i)),']'];
    line2 = genegRNA{i};
    
    if i == 1
        output = [output; line1; line2];
    else
        output = [output; emptyLine; line1; line2];
    end
end

file = fopen(fileOutName,'w');
for i = [1:length(output(:,1))]
    fprintf(file,'%s\n',output{i});
end
fclose(file);

clearvars emptyLine fileOutName genegRNA geneName geneNames geneNCBI i line1
clearvars line2 output pickOrder refSeq refSeqName file

%% Acquires counts & merges lanes

parent = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\03.31.19_ICT Analysis\Data\2020.07.03_LibraryCloning\NGSdata\2020.09.15_NGS\C3library\counts';
sampleKeyWords = {'C3-M';
    'C3-P'};

% Combines files
files = dir(parent);
countOut = [];
gRNAOut = {};
sampOut = {};
for i = [1:length(files(:,1))-4]
    fileName = files(i).name;
    
    loc = strfind(fileName,'.txt');
    if isempty(loc)
        continue
    end
    
    address = [parent,'\',fileName];
    data = readtable(address);
    data = table2cell(data);
    data(1,:) = [];
    
    counts = cell2mat(data(:,1));
    gRNA = data(:,2);
    fileName = strsplit(fileName,'.');
    fileName = fileName{1};
    
    countOut = [countOut, counts];
    gRNAOut = [gRNAOut, gRNA];
    sampOut = [sampOut, fileName];
end

% Merges lanes
countOut2 = [];
for i = [1:length(sampleKeyWords)]
    keyWord = sampleKeyWords{i};
    
    ind = [];
    for j = [1:length(sampOut)]
        loc = strfind(sampOut{j},keyWord);
        if not(isempty(loc))
            ind = [ind, j];
        end
    end
    
    % Adds lanes
    val = sum(countOut(:,ind),2);
    countOut2 = [countOut2, val];
end

clear address countOut counts data fileName files gRNA i
clear ind j keyWord loc parent sampleKeyWords val

%% Creates histogram

data1 = input1;
data2 = input2;
data3 = input3;

data1 = log10(data1);
data2 = log10(data2);
data3 = log10(data3);

figure
hold
h = histogram(data1,20);
h.FaceColor = [0.5 1 0.5];
set(gca,'LineWidth',3);
set(gca,'TickDir','out');
set(gca,'FontSize',25);
xlim([0,5]);

figure
hold
h = histogram(data2,20);
h.FaceColor = [0.5 0.7 1];
set(gca,'LineWidth',3);
set(gca,'TickDir','out');
set(gca,'FontSize',25);
xlim([0,5]);

figure
hold
h = histogram(data3,35);
h.FaceColor = [0.5 0 0.2];
set(gca,'LineWidth',3);
set(gca,'TickDir','out');
set(gca,'FontSize',25);
xlim([0,5]);

clear data1 data2 data3 h

%% Correlation plot

data1 = input1;
data2 = input2;
sColor = [0.5 0 0.2];

data1 = log10(data1);
data2 = log10(data2);
[p ~] = corr(data1,data2);

figure
hold
s = scatter(data1,data2,'o');
s.MarkerFaceColor = sColor;
s.MarkerFaceAlpha = 0.3;
s.MarkerEdgeColor = [0 0 0];
s.MarkerEdgeAlpha = 0.3;
s.SizeData = 30;
set(gca,'LineWidth',3);
set(gca,'TickDir','out');
set(gca,'FontSize',25);

clear data1 data2 sColor s

%% Acquire range of 95th & 5th percentile gRNAs

data = input;

upper = prctile(input,95);
lower = prctile(input,5);
range = [upper,lower];
range = log10(range);
range

clear data upper lower







