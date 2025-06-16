%% Acquires proportion of newICTs represented in each FirstTier outputs

struc = list;
dlist = newICTs;

% Acquires the 1st genes in 'dlist'
for i = [1:length(dlist(:,1))]
    aliases = strsplit(dlist{i,1},' ');
    dlist{i,1} = aliases{1};
end

% Acquires 'dlist' representation for each FirstTier cutoff
listnum = fields(struc);
output = {};
for i = [1:length(listnum(:,1))]
    eval(['data = struc.',listnum{i,1},'.FirstTier(:,1);']);
    
    genesRepresented = {};
    for j = [1:length(data(:,1))]
        gene = data{j,1};
        loc = find(strcmp(lower(dlist), lower(gene)));
        if not(isempty(loc))
            genesRepresented = [genesRepresented; gene];
        end
    end
    
    eval(['output.',listnum{i,1},' = genesRepresented;']);
end

% Plots percentage of representation by cutoff
val = [];
for i = [1:length(listnum(:,1))]
    eval(['data = output.',listnum{i,1},';']);
    val = [val; length(data)];
end
yval = val;
xval = 1:length(yval);
p = plot(xval,yval)
p.LineWidth = 3;
p.Color = [0 0 0];
set(gca,'LineWidth',3)
set(gca,'box','off')
set(gca,'TickDir','out')
set(gca,'FontSize',14)
set(gca,'XGrid','on')
set(gca,'YGrid','on')

clear aliases data dlist gene genesRepresented i j listnum
clear loc p struc val xval yval


%% Plots FirstTier number against percentile cutoff

struc = list;

listnum = fields(struc);
val = [];
for i = [1:length(listnum(:,1))]
    eval(['data = struc.',listnum{i,1},'.FirstTier;']);
    val = [val; length(data)];
end
yval = val;
xval = 1:length(yval);
p = plot(xval,yval)
p.LineWidth = 3;
p.Color = [0 0 0];
set(gca,'LineWidth',3)
set(gca,'box','off')
set(gca,'TickDir','out')
set(gca,'FontSize',14)
set(gca,'XGrid','on')
set(gca,'YGrid','on')

%% Acquires metrics for input fields

struc = list;
newICTlist = newICTs;
v1List = v1;
inputFields = {'list25',
    'list30',
    'list35',
    'list40'};

output1 = [];
output2 = [];
for i = [1:length(inputFields(:,1))]
    fieldName = inputFields{i,1};
    eval(['data = struc.',fieldName,'.FirstTier(:,1);']);
    
    % Acquires proportion of FirstTier genes in newICTlist
    newGenesInList = {};
    for j = [1:length(newICTlist(:,1))]
        aliases = strsplit(newICTlist{j,1},' ');
        
        for k = [1:length(aliases)]
            gene = aliases{k};
            loc = find(strcmp(lower(data), lower(gene)));
            if not(isempty(loc))
                newGenesInList = [newGenesInList; aliases{1}];
                break
            end
        end
    end
    output1 = [output1; length(newGenesInList)/length(newICTlist)];
    
    % Acquires proportion of FirstTier genes in v1List
    V1GenesInList = {};
    for j = [1:length(v1List(:,1))]
        gene = v1List{j,1};
        loc = find(strcmp(lower(data), lower(gene)));
        if not(isempty(loc))
            V1GenesInList = [V1GenesInList; gene];
        end
    end
    output2 = [output2; length(V1GenesInList)/length(v1List)];
end

clear aliases data fieldName gene i inputFields j k loc newGenesInList
clear newICTlist struc V1GenesInList v1List







